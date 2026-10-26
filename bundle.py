"""Bundle-ready detail-map data for the app (written to out/bundle; nothing is written outside work/ and out/).

usage: python -I bundle.py
reads  work/osm_features.pkl (osm_extract.py), work/terrain/*.json (terrain.py), work/land.pkl and
       work/korea_crop.npz (land.py: everything is cropped to Korea; the boundary rivers stay whole),
       the general-map scene MapBase.json (read only; $DETAILMAP_MAPBASE, else ../../Daedong/Resources/Map/MapBase.json,
       else data/MapBase.json: see common.mapbase_path)
writes out/bundle/peaks_{A,B,C}.bin   glyph lists per zoom band (A ~4x, B ~8x, C ~16x)
       out/bundle/lines_{A,B}.bin     coastline, rivers (A) and streams (B), lakes and wide rivers
       out/bundle/ridges.bin          major divides for the faint connecting ridge wash
       out/bundle/land_{A,B}.bin      Korean land polygons per tile (land.py), and detail_area.png from land.py
       out/bundle/labels.bin          peaks, streams and rivers, 시·군·읍 names (Korean)
       out/general/glyph_moves.json   the general map's hand-drawn mountains moved onto real summits
       out/general/MapBase_moved.json the same, applied to a copy of the general scene (for previews)
       out/bundle/chars.txt, chars_missing.txt  characters the labels need / not in the current font subset
       out/sizes.json
       work/bands.pkl, work/peaks.json, work/ridges.json   (check data for a preview renderer; not used here)

File format (all little endian): see geom.pack. Tile payloads:
  peaks: count, then per glyph: dx, dy (zigzag varint, 1/64 pt, sorted by y), kind (0 sharp, 1 round, 2 mound,
         3 spire of the previous head), height in 1/16 map pt at the band's level zoom, zoom-in x8, zoom-full x8,
         flags (1 famous, 2 rough rock: add flanking spires, 4 has name)
  lines: count, then per feature: class, parts, points (zigzag varint deltas; A 1/64 pt, B 1/256 pt)
  ridges: like lines plus one weight byte per point
  labels (one deflated stream): count, then per label: kind (1 peak, 2 river, 3 stream, 4 city, 5 town/county seat),
         min zoom x8, x, y (varint, 1/64 pt), [peaks: elevation varint] [streams: angle byte], name length, UTF-8
"""
import glob, json, math, os, pickle, struct, sys, zlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import WORK, OUT, HERE, meters_per_point, unproj, mapbase_path  # noqa: E402
from geom import rdp, clip_poly, varint, zz, encode_tile, pack, area_km2, deflate  # noqa: E402
import glyphs as G  # noqa: E402

BAND_LEVEL = {'A': 4.0, 'B': 8.0, 'C': 16.0}     # the detail layer's tile levels (zoom 2~4, 4~8, 8~16)
LINE_TILE = {'A': 128.0, 'B': 32.0}
LINE_Q = {'A': 64.0, 'B': 256.0}
TOL_M = {'A': 25.0, 'B': 6.0}
MPP = meters_per_point(38.0)


def split_tiles(pts, T):
    """Cut a polyline at tile borders: {(i, j): [piece, ...]}."""
    out = {}
    cur, ct = [pts[0]], (int(pts[0][0] // T), int(pts[0][1] // T))
    for a, b in zip(pts, pts[1:]):
        tb = (int(b[0] // T), int(b[1] // T))
        if tb == ct:
            cur.append(b)
            continue
        ts = []
        for axis in (0, 1):
            lo, hi = sorted((a[axis], b[axis]))
            k = math.floor(lo / T) + 1
            while k * T <= hi:
                if b[axis] != a[axis]:
                    ts.append((k * T - a[axis]) / (b[axis] - a[axis]))
                k += 1
        for t in sorted(ts):
            q = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
            cur.append(q)
            if len(cur) > 1:
                out.setdefault(ct, []).append(cur)
            nxt = min([u for u in ts if u > t] + [1.0])
            mq = (a[0] + (t + nxt) / 2 * (b[0] - a[0]), a[1] + (t + nxt) / 2 * (b[1] - a[1]))
            ct = (int(mq[0] // T), int(mq[1] // T))
            cur = [q]
        cur.append(b)
        ct = tb
    if len(cur) > 1:
        out.setdefault(ct, []).append(cur)
    return out


def simp(xy, tol_m):
    pts = [tuple(p) for p in xy.tolist()]
    return rdp(pts, tol_m / MPP) if len(pts) > 2 else pts


def main():
    os.makedirs(os.path.join(OUT, 'bundle'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'general'), exist_ok=True)
    feats = pickle.load(open(os.path.join(WORK, 'osm_features.pkl'), 'rb'))
    crop = np.load(os.path.join(WORK, 'korea_crop.npz'))
    CM, CP = crop['mask'], float(crop['p'])
    BORDER = ('압록강', '두만강')

    def inside(x, y):
        r, c = int(y * CP), int(x * CP)
        return 0 <= r < CM.shape[0] and 0 <= c < CM.shape[1] and bool(CM[r, c])

    def crop_line(xy):
        """Pieces of a polyline inside Korea."""
        r = np.clip((xy[:, 1] * CP).astype(int), 0, CM.shape[0] - 1)
        c = np.clip((xy[:, 0] * CP).astype(int), 0, CM.shape[1] - 1)
        ins = CM[r, c]
        out, cur = [], []
        for p, k in zip(xy, ins):
            if k:
                cur.append(p)
            elif cur:
                if len(cur) > 1:
                    out.append(np.array(cur))
                cur = []
        if len(cur) > 1:
            out.append(np.array(cur))
        return out
    n0 = (len(feats['coast']), len(feats['water']), len(feats['areas']), len(feats['pois']))
    feats['coast'] = [p for xy in feats['coast'] for p in crop_line(xy)]
    feats['water'] = [dict(w, xy=w['xy']) if w['name'] in BORDER else None for w in feats['water']] and \
        [w for w in feats['water'] if w['name'] in BORDER] + \
        [dict(w, xy=p) for w in feats['water'] if w['name'] not in BORDER for p in crop_line(w['xy'])]
    feats['areas'] = [a for a in feats['areas'] if a['name'] in BORDER or
                      inside(*np.asarray(a['rings'][0]).mean(axis=0))]
    feats['pois'] = [q for q in feats['pois'] if inside(q['x'], q['y'])]
    print('cropped to Korea: coast/water/areas/pois', n0, '->',
          (len(feats['coast']), len(feats['water']), len(feats['areas']), len(feats['pois'])))
    # ---------------- summits
    peaks, ridges = [], []
    for f in sorted(glob.glob(os.path.join(WORK, 'terrain', '*.json'))):
        d = json.load(open(f))
        peaks += d['peaks']
        ridges += d['ridges']
    peaks = [p for p in peaks if inside(p['x'], p['y'])]
    ridges = [r for r in ridges if inside(r['pts'][len(r['pts']) // 2 // 2 * 2], r['pts'][len(r['pts']) // 2 // 2 * 2 + 1])]
    grid = G.Grid([(p['x'], p['y']) for p in peaks])
    named = [q for q in feats['pois'] if q['kind'] == 'peak']
    taken = set()
    for l in sorted(named, key=lambda l: -(l.get('ele') or 0)):
        r = 300 / MPP
        best, bp = None, -1
        for k in grid(l['x'] - r, l['y'] - r, l['x'] + r, l['y'] + r):
            p = peaks[k]
            if k in taken or math.hypot(p['x'] - l['x'], p['y'] - l['y']) > r:
                continue
            if p['prom'] > bp:
                best, bp = k, p['prom']
        if best is not None:
            taken.add(best)
            peaks[best]['name'] = l['name']
            peaks[best]['wd'] = l['wd']
    G.clusters(peaks)
    json.dump(peaks, open(os.path.join(WORK, 'peaks.json'), 'w'), ensure_ascii=False)
    json.dump(ridges, open(os.path.join(WORK, 'ridges.json'), 'w'))
    sizes = {'peaks_total': len(peaks), 'named_peaks': len(taken)}
    counts = {}
    for band, level in BAND_LEVEL.items():
        zm = level / math.sqrt(2)          # a level's tiles are seen from level/2 to level; sizes are set for the middle
        order = sorted((i for i, p in enumerate(peaks) if 'head' not in p and G.zoom_in(G.eff_prom(p)) < zm),
                       key=lambda i: -G.eff_prom(peaks[i]))
        boxes = {}
        kept = []
        for i in order:
            p = peaks[i]
            k = G.kind(p)
            s_ = G.smoothstep(G.zoom_in(G.eff_prom(p)), G.zoom_full(G.eff_prom(p)), zm)
            Hs = G.screen_height(p, zm, k) * (0.55 + 0.45 * s_)
            W = Hs * G.TEMPLATES[k]['aspect']
            cx, cy = int(p['x'] // 4), int(p['y'] // 4)
            hidden = False
            for gx in (cx - 1, cx, cx + 1):
                for gy in (cy - 1, cy, cy + 1):
                    for bx0, by0, bx1, by1 in boxes.get((gx, gy), ()):
                        if bx0 < p['x'] < bx1 and by0 < p['y'] < by1:
                            hidden = True
            if hidden:
                continue
            box = (p['x'] - W * 0.35 / zm, p['y'] - 2 / zm, p['x'] + W * 0.35 / zm, p['y'] + Hs * 0.55 / zm)
            for gx in range(int(box[0] // 4), int(box[2] // 4) + 1):
                for gy in range(int(box[1] // 4), int(box[3] // 4) + 1):
                    boxes.setdefault((gx, gy), []).append(box)
            kept.append(i)
        tiles = {}
        n = 0
        for i in kept:
            p = peaks[i]
            k = G.kind(p)
            Hfull = G.screen_height(p, zm, k)
            W = Hfull * G.TEMPLATES[k]['aspect']
            entries = [(p, k, Hfull)]
            if k == 'sharp':
                for j in p.get('members', []):
                    q = peaks[j]
                    if abs(q['x'] - p['x']) * zm > 1.3 * W or abs(q['y'] - p['y']) * zm > 0.8 * Hfull:
                        continue
                    entries.append((q, 'spire', Hfull * min(0.85, 0.5 + 0.3 * math.sqrt(q['prom'] / max(p['prom'], 1)))))
            for q, kk, H in entries:
                tiles.setdefault((int(p['x'] // 64), int(p['y'] // 64)), []).append((q, kk, H, p))
                n += 1
        payload = {}
        for (i, j), items in tiles.items():
            b = bytearray()
            varint(len(items), b)
            px = py = 0
            for q, kk, H, head in items:
                ix, iy = int(round((q['x'] - i * 64) * 64)), int(round((q['y'] - j * 64) * 64))
                varint(zz(ix - px), b)
                varint(zz(iy - py), b)
                px, py = ix, iy
                b.append(G.KIND_ID[kk])
                varint(int(round(H / zm * 16)), b)            # full height in 1/16 map pt
                b.append(min(255, int(round(G.zoom_in(G.eff_prom(head)) * 8))))
                b.append(min(255, int(round(G.zoom_full(G.eff_prom(head)) * 8))))
                fam = (head['ele'] >= 1200 and head['prom'] >= 250) or (head.get('wd') and head['prom'] >= 150)
                nm = sum(1 for j2 in head.get('members', [])) if kk != 'spire' else 0
                b.append((1 if fam and kk != 'spire' else 0) | (2 if kk == 'sharp' and head['rough'] >= 8 and nm < 2 else 0)
                         | (4 if q.get('name') and kk != 'spire' else 0))
            payload[(i, j)] = bytes(b)
        data = pack(payload, b'DDP1')
        open(os.path.join(OUT, 'bundle', f'peaks_{band}.bin'), 'wb').write(data)
        sizes[f'peaks_{band}.bin'] = len(data)
        counts[band] = n
    sizes['glyphs_per_band'] = counts
    # ---------------- lines and lakes per band
    bands = {}
    for band in ('A', 'B'):
        T, q, tol = LINE_TILE[band], LINE_Q[band], TOL_M[band]
        tiles = {}
        bf = {'coast': [], 'river': [], 'stream': [], 'water': []}
        for xy in feats['coast']:
            if len(xy) > 1:
                bf['coast'].append(simp(xy, tol))
        for w in feats['water']:
            if len(w['xy']) < 2 or (band == 'A' and w['kind'] != 'river'):
                continue
            bf['river' if w['kind'] == 'river' else 'stream'].append(simp(w['xy'], tol))
        for a in feats['areas']:
            rings = [r for r in a['rings'] if len(r) >= 4]
            if not rings:
                continue
            ar = area_km2([tuple(p) for p in rings[0].tolist()])
            if band == 'A' and ar < 0.05:
                continue
            rs = [simp(r, tol) for r in rings]
            rs = [r for r in rs if len(r) >= 3]
            if rs:
                bf['water'].append({'rings': rs, 'river': a['river'], 'area': ar, 'name': a['name']})
        bands[band] = bf
        for cls in ('coast', 'river', 'stream'):
            for line in bf[cls]:
                for key, pieces in split_tiles(line, T).items():
                    for pc in pieces:
                        tiles.setdefault(key, []).append((cls, [pc]))
        for wa in bf['water']:
            xs = [p[0] for p in wa['rings'][0]]
            ys = [p[1] for p in wa['rings'][0]]
            cls = 'riverarea' if wa['river'] else 'water'
            for j in range(int(min(ys) // T), int(max(ys) // T) + 1):
                for i in range(int(min(xs) // T), int(max(xs) // T) + 1):
                    if int(min(xs) // T) == int(max(xs) // T) and int(min(ys) // T) == int(max(ys) // T):
                        parts = wa['rings']
                    else:
                        parts = [clip_poly(r, i * T, j * T, i * T + T, j * T + T) for r in wa['rings']]
                        parts = [pp for pp in parts if len(pp) >= 3]
                    if parts:
                        tiles.setdefault((i, j), []).append((cls, parts))
        payload = {k: encode_tile(v, k[0] * T, k[1] * T, q) for k, v in tiles.items()}
        data = pack(payload, b'DDL1')
        open(os.path.join(OUT, 'bundle', f'lines_{band}.bin'), 'wb').write(data)
        sizes[f'lines_{band}.bin'] = len(data)
        sizes[f'lines_{band}_features'] = {k: len(v) for k, v in bf.items()}
    pickle.dump(bands, open(os.path.join(WORK, 'bands.pkl'), 'wb'))
    # ---------------- ridge wash lines
    tiles = {}
    for r in ridges:
        if max(r['w']) < 0.45:
            continue
        pts = [(r['pts'][i], r['pts'][i + 1]) for i in range(0, len(r['pts']) - 1, 2)]
        key = (int(pts[len(pts) // 2][0] // 128), int(pts[len(pts) // 2][1] // 128))
        tiles.setdefault(key, []).append((pts, r['w']))
    payload = {}
    for (i, j), ls in tiles.items():
        b = bytearray()
        varint(len(ls), b)
        for pts, w in ls:
            varint(len(pts), b)
            px = py = 0
            for (x, y), wi in zip(pts, w):
                ix, iy = int(round((x - i * 128) * 64)), int(round((y - j * 128) * 64))
                varint(zz(ix - px), b)
                varint(zz(iy - py), b)
                px, py = ix, iy
                b.append(min(255, int(round(wi * 255))))
        payload[(i, j)] = bytes(b)
    data = pack(payload, b'DDR1')
    open(os.path.join(OUT, 'bundle', 'ridges.bin'), 'wb').write(data)
    sizes['ridges.bin'] = len(data)
    # ---------------- Korean land polygons per tile
    land = pickle.load(open(os.path.join(WORK, 'land.pkl'), 'rb'))
    for band in ('A', 'B'):
        T, q, tol = LINE_TILE[band], LINE_Q[band], TOL_M[band]
        tiles = {}

        def qclip(poly, x0, y0, size):
            if size <= T + 1e-6:
                tiles.setdefault((int(round(x0 / T)), int(round(y0 / T))), []).append(('land', [poly]))
                return
            h = size / 2
            for qx, qy in ((x0, y0), (x0 + h, y0), (x0, y0 + h), (x0 + h, y0 + h)):
                c = clip_poly(poly, qx, qy, qx + h, qy + h)
                if len(c) >= 3:
                    qclip(c, qx, qy, h)
        rings = [land['mainland']] + [r for r in land['islands'] if len(r) >= 4]
        for r in rings:
            pr = simp(r, tol)
            if len(pr) < 3:
                continue
            xs = [p[0] for p in pr]
            ys = [p[1] for p in pr]
            i0, i1, j0, j1 = int(min(xs) // T), int(max(xs) // T), int(min(ys) // T), int(max(ys) // T)
            if i0 == i1 and j0 == j1:
                tiles.setdefault((i0, j0), []).append(('land', [pr]))
            elif (i1 - i0 + 1) * (j1 - j0 + 1) <= 16:
                for j in range(j0, j1 + 1):
                    for i in range(i0, i1 + 1):
                        c = clip_poly(pr, i * T, j * T, i * T + T, j * T + T)
                        if len(c) >= 3:
                            tiles.setdefault((i, j), []).append(('land', [c]))
            else:
                size = T
                while size < max(1360, 2400):
                    size *= 2
                qclip(pr, 0.0, 0.0, size)
        payload = {k: encode_tile(v, k[0] * T, k[1] * T, q) for k, v in tiles.items()}
        data = pack(payload, b'DDG1')
        open(os.path.join(OUT, 'bundle', f'land_{band}.bin'), 'wb').write(data)
        sizes[f'land_{band}.bin'] = len(data)
        sizes[f'land_{band}_tiles'] = len(payload)
    # detail area (land.py) as 64-pt tiles of 256 x 256 bytes: missing = nothing, empty payload = all in
    from PIL import Image
    area = np.array(Image.open(os.path.join(WORK, 'detail_area.png')))
    payload = {}
    for j in range(area.shape[0] // 256):
        for i in range(area.shape[1] // 256):
            sub = area[j * 256:(j + 1) * 256, i * 256:(i + 1) * 256]
            if sub.max() == 0:
                continue
            payload[(i, j)] = b'' if sub.min() == 255 else sub.tobytes()
    data = pack(payload, b'DDA1')
    open(os.path.join(OUT, 'bundle', 'area.bin'), 'wb').write(data)
    sizes['area.bin'] = len(data)
    sizes['area_tiles'] = {'partial': sum(1 for v in payload.values() if v), 'full': sum(1 for v in payload.values() if not v)}
    # ---------------- labels
    labels = []
    for p in peaks:
        if p.get('name') and 'head' not in p and p['name'] != '삼각점' and not any(c.isdigit() for c in p['name']):
            labels.append((1, max(3.0, G.zoom_in(G.eff_prom(p))), p['x'], p['y'], p['name'], int(p['ele']), None))
    best = {}
    for w in feats['water']:
        if not w['name'] or len(w['xy']) < 2:
            continue
        xy = w['xy']
        L = float(((xy[1:] - xy[:-1]) ** 2).sum(axis=1).__pow__(0.5).sum())
        key = (w['name'], w['kind'] == 'river')
        if key not in best or L > best[key][0]:
            best[key] = (L, xy, w['kind'])
    for (name, is_river), (L, xy, kind) in best.items():
        m = len(xy) // 2
        a, b = xy[max(0, m - 1)], xy[min(len(xy) - 1, m + 1)]
        ang = math.degrees(math.atan2(b[1] - a[1], b[0] - a[0])) % 180
        labels.append((2 if is_river else 3, 6.0 if is_river else 12.0, float(xy[m][0]), float(xy[m][1]), name, None, ang))
    for q in feats['pois']:
        if q['kind'] != 'place':
            continue
        n = q['name']
        if q['place'] == 'city' or n.endswith('시'):
            labels.append((4, 2.5, q['x'], q['y'], n, None, None))
        elif q['place'] == 'town' or n.endswith('군') or n.endswith('읍'):
            labels.append((5, 4.0, q['x'], q['y'], n, None, None))
    # Korean names only: features across the Yalu and Tumen come with Chinese or Russian names
    def korean(n):
        return any(0xAC00 <= ord(c) <= 0xD7A3 for c in n) and not any(
            0x0400 <= ord(c) <= 0x04FF or (0x4E00 <= ord(c) <= 0x9FFF and not any(0xAC00 <= ord(h) <= 0xD7A3 for h in n.split('(')[0]))
            for c in n)
    labels = [l for l in labels if korean(l[4])]
    b = bytearray()
    varint(len(labels), b)
    chars = set()
    for kind, mz, x, y, name, ele, ang in labels:
        b.append(kind)
        b.append(min(255, int(round(mz * 8))))
        varint(int(round(x * 64)), b)
        varint(int(round(y * 64)), b)
        if kind == 1:
            varint(max(0, ele), b)
        if kind in (2, 3):
            b.append(int(ang / 180 * 255) & 255)
        nb = name.encode('utf-8')
        varint(len(nb), b)
        b += nb
        chars |= set(name)
    data = deflate(bytes(b))
    open(os.path.join(OUT, 'bundle', 'labels.bin'), 'wb').write(data)
    sizes['labels.bin'] = len(data)
    sizes['labels_by_kind'] = {k: sum(1 for l in labels if l[0] == v) for k, v in
                               (('peak', 1), ('river', 2), ('stream', 3), ('city', 4), ('county/town', 5))}
    with open(os.path.join(WORK, 'labels.tsv'), 'w') as f:
        for l in labels:
            f.write('\t'.join(str(v) for v in l) + '\n')
    # ---------------- characters for the font subset
    ks = set()
    for c in range(0xAC00, 0xD7A4):
        try:
            if len(chr(c).encode('euc-kr')) == 2:
                ks.add(chr(c))
        except UnicodeEncodeError:
            pass
    base = ks | {chr(c) for c in range(0x20, 0x7F)} | set('·~…‘’“”「」『』〈〉《》–年月日％')
    app_src = os.path.join(HERE, '..', '..', 'Daedong')
    if os.path.isdir(app_src):
        for f in glob.glob(os.path.join(app_src, '**', '*.swift'), recursive=True):
            base |= set(open(f, encoding='utf-8').read())
    else:
        print('app sources (../../Daedong) not found: chars_missing is checked against Hangul (KS X 1001) and ASCII only')
    need = ''.join(sorted(c for c in chars if ord(c) >= 0x20))
    missing = ''.join(sorted(c for c in chars if c not in base and ord(c) >= 0x20))
    open(os.path.join(OUT, 'bundle', 'chars.txt'), 'w').write(need)
    open(os.path.join(OUT, 'bundle', 'chars_missing.txt'), 'w').write(missing)
    sizes['chars'] = len(need)
    sizes['chars_missing'] = missing
    # ---------------- general map: hand-drawn mountains onto real summits
    base = json.load(open(mapbase_path()))
    gl = next(o for o in base['ops'] if o['t'] == 'glyphs' and o.get('items'))
    orig_keys = [tuple(g['p'][:2]) for g in gl['items']]
    if base.get('detailApplied'):
        print('general map already carries the moves (detailApplied is set in the scene); keeping glyph_moves.json')
        moves = [tuple(m) for m in json.load(open(os.path.join(OUT, 'general', 'glyph_moves.json')))['moves']]
    else:
        moves = G.move_general(base, peaks, index=grid)
    json.dump({'note': 'index = position of the glyph item in the original MapBase.json glyph list; apply dx, dy to '
                       'p, x and d (x, y pairs; d is x, y, r), then re-sort the items by base line '
                       '(mean of p[1] and p[13]). Done by apply in glyphs.move_general.',
               'moves': [[k, dx, dy] for k, dx, dy in moves]},
              open(os.path.join(OUT, 'general', 'glyph_moves.json'), 'w'))
    if not base.get('detailApplied'):
        json.dump(base, open(os.path.join(OUT, 'general', 'MapBase_moved.json'), 'w'), separators=(',', ':'))
    sizes['general_glyphs_moved'] = f'{len(moves)} of {len(orig_keys)}'
    total = sum(v for k, v in sizes.items() if k.endswith('.bin') or k.endswith('.png'))
    sizes['bundle_total_bytes'] = total
    json.dump(sizes, open(os.path.join(OUT, 'sizes.json'), 'w'), indent=1, ensure_ascii=False)
    print(json.dumps(sizes, indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
