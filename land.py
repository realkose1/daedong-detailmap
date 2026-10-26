"""Korean land polygons for the detail layer, and the masks that keep everything on the Korean side.

usage: python -I land.py
reads  work/osm_features.pkl, work/territory.pkl, the general-map scene MapBase.json (general coastline; $DETAILMAP_MAPBASE,
       else ../../Daedong/Resources/Map/MapBase.json, else data/MapBase.json: see common.mapbase_path)
writes work/land.pkl            {'mainland': ring, 'islands': [rings]} in map points
       work/korea_crop.npz      Korean land (+300 m) at 8 cells per point, for cropping lines, labels, summits
       work/detail_area.png     where the detail layer draws (bundle.py tiles it into out/bundle/area.bin): Korean land and the sea within 25 km of it
                                (fading out from 15 km), 4 px per map point, 8-bit; nothing over foreign land
       out/land_check.json      how far the general map's coastline sits from the OSM coastline

The mainland ring: OSM coastline from the 압록강 mouth round the peninsula to the 두만강 mouth, then up the
두만강 centre line, across the 백두산 watershed, and down the 압록강 centre line (the land border follows the
rivers). Islands: closed coastline rings on the Korean side (inside South Korea's territory, or near the
North Korean coast and not across the border).
"""
import json, math, os, pickle, sys
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import WORK, OUT, HERE, unproj, proj, meters_per_point, mapbase_path  # noqa: E402

KEY = lambda p: (round(float(p[0]), 5), round(float(p[1]), 5))


def join(lines):
    """Join polylines end to end (shared end points). Returns (closed rings, open chains)."""
    pool = [np.asarray(l, dtype=np.float64) for l in lines if len(l) >= 2]
    ends = {}
    for i, l in enumerate(pool):
        ends.setdefault(KEY(l[0]), set()).add(i)
        ends.setdefault(KEY(l[-1]), set()).add(i)
    used = [False] * len(pool)
    rings, chains = [], []
    for i in range(len(pool)):
        if used[i]:
            continue
        used[i] = True
        cur = pool[i]
        for direction in (1, -1):
            while True:
                tail = KEY(cur[-1]) if direction == 1 else KEY(cur[0])
                if KEY(cur[0]) == KEY(cur[-1]) and len(cur) > 3:
                    break
                nxt = [j for j in ends.get(tail, ()) if not used[j]]
                if not nxt:
                    break
                j = nxt[0]
                used[j] = True
                l = pool[j]
                if direction == 1:
                    cur = np.vstack([cur, l[1:]]) if KEY(l[0]) == tail else np.vstack([cur, l[::-1][1:]])
                else:
                    cur = np.vstack([l[:-1], cur]) if KEY(l[-1]) == tail else np.vstack([l[::-1][:-1], cur])
        (rings if KEY(cur[0]) == KEY(cur[-1]) and len(cur) > 3 else chains).append(cur)
    return rings, chains


def area(r):
    x, y = r[:, 0], r[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def nearest(chain, p):
    d = ((chain - p) ** 2).sum(axis=1)
    k = int(np.argmin(d))
    return k, math.sqrt(d[k])


def main():
    feats = pickle.load(open(os.path.join(WORK, 'osm_features.pkl'), 'rb'))
    terr = pickle.load(open(os.path.join(WORK, 'territory.pkl'), 'rb'))
    rings, chains = join(feats['coast'])
    chains.sort(key=lambda c: -len(c))
    print('coast rings', len(rings), 'open chains', len(chains), [len(c) for c in chains[:6]])
    # boundary rivers: centre lines from the mouth up to the source (ways may not share end nodes at dams)
    main_chain = chains[0]

    def river_path(name):
        rs, cs = join([w['xy'] for w in feats['water'] if w['name'] == name and w['kind'] == 'river'])
        pool = [c for c in cs + rs if len(c) >= 5]
        ends = [(nearest(main_chain, c[0])[1], k, 0) for k, c in enumerate(pool)] + \
               [(nearest(main_chain, c[-1])[1], k, 1) for k, c in enumerate(pool)]
        _, k, e = min(ends)
        path = pool.pop(k)
        if e == 1:
            path = path[::-1]                      # mouth first
        while pool:
            best = min(((math.hypot(*(c[0] - path[-1])), k, 0) for k, c in enumerate(pool)), default=None)
            best2 = min(((math.hypot(*(c[-1] - path[-1])), k, 1) for k, c in enumerate(pool)), default=None)
            d, k, e = min(best, best2)
            if d > 12:                             # ~6 km
                break
            c = pool.pop(k)
            path = np.vstack([path, c if e == 0 else c[::-1]])
        return path
    yalu, tumen = river_path('압록강'), river_path('두만강')     # both mouth first
    ky, dy = nearest(main_chain, yalu[0])
    kt, dt = nearest(main_chain, tumen[0])
    print('mouth gaps (pt): 압록강', round(dy, 2), '두만강', round(dt, 2), 'sources',
          [round(v, 3) for v in unproj(*yalu[-1])], [round(v, 3) for v in unproj(*tumen[-1])])
    a, b = sorted((ky, kt))
    coast_part = main_chain[a:b + 1]
    if coast_part[:, 1].max() < main_chain[:, 1].max() - 1:   # must go round the south coast
        coast_part = np.vstack([main_chain[b:], main_chain[:a + 1]])
    if math.hypot(*(coast_part[0] - yalu[0])) > math.hypot(*(coast_part[-1] - yalu[0])):
        coast_part = coast_part[::-1]                          # start at the 압록강 mouth
    ring = np.vstack([coast_part, tumen, yalu[::-1]])
    ring = np.vstack([ring, ring[:1]])
    yalu = yalu[::-1]
    tumen = tumen[::-1]
    print('mainland ring points', len(ring), 'area pt2', round(area(ring)))
    # islands on the Korean side
    kr = max(terr.get('KR', []), key=len, default=None)
    P = 2
    def rasterize(rs, scale, w=1360, h=2400):
        im = Image.new('L', (int(w * scale), int(h * scale)), 0)
        dr = ImageDraw.Draw(im)
        for r in rs:
            dr.polygon([(x * scale, y * scale) for x, y in r.tolist()], fill=255)
        return np.array(im) > 0
    krmask = rasterize(terr.get('KR', []), P)
    main_mask = rasterize([ring], P)
    dist_main = ndi.distance_transform_edt(~main_mask) / P          # map points
    yalu_mouth, tumen_mouth = yalu[-1], tumen[-1]
    islands = []
    for r in rings:
        c = r.mean(axis=0)
        ci, cj = int(c[1] * P), int(c[0] * P)
        if not (0 <= ci < krmask.shape[0] and 0 <= cj < krmask.shape[1]):
            continue
        lo, la = unproj(*c)
        if krmask[ci, cj]:
            islands.append(r)
            continue
        near = dist_main[ci, cj] * meters_per_point(la) < 40000
        across = (c[0] < yalu_mouth[0] - 2 and c[1] < yalu_mouth[1] + 30) or (c[0] > tumen_mouth[0] + 1)
        japan = la < 35.2 and lo > 128.9
        if near and not across and not japan:
            islands.append(r)
    print('islands', len(islands), 'of', len(rings))
    pickle.dump({'mainland': ring, 'islands': islands, 'yalu': yalu, 'tumen': tumen, 'coast_part': coast_part},
                open(os.path.join(WORK, 'land.pkl'), 'wb'))
    # crop mask (Korean land + 300 m), 8 per point
    C = 8
    land8 = rasterize([ring] + islands, C)
    buf = int(round(300 / meters_per_point(38) * C))
    crop = ndi.binary_dilation(land8, iterations=buf)
    np.savez_compressed(os.path.join(WORK, 'korea_crop.npz'), mask=crop, p=C)
    # detail area: Korean land + sea within 25 km (fade from 15 km), never over the land beyond the border.
    # "Beyond" = north of the boundary rivers, closed off seaward by short lines from the two river mouths
    # (west-south-west of the 압록강 mouth, east-south-east of the 두만강 mouth) up to the top of the map.
    P = 4
    land2 = rasterize([ring] + islands, P)
    d = ndi.distance_transform_edt(~land2) / P * meters_per_point(37.5) / 1000      # km
    a = np.clip((25 - d) / 10, 0, 1)
    a = a * a * (3 - 2 * a)
    ym, tm = yalu[-1], tumen[-1]            # mouths (yalu and tumen are source -> mouth here)
    wedge = np.vstack([[ym[0] - 30, ym[1] + 8], [ym[0] - 30, -50], [tm[0] + 40, -50], [tm[0] + 40, tm[1] + 10],
                       tm[None, :], tumen[::-1], [[(yalu[0][0] + tumen[0][0]) / 2, (yalu[0][1] + tumen[0][1]) / 2]],
                       yalu, ym[None, :]])
    beyond = rasterize([wedge], P)
    # soft seaward edge of the wedge (no straight seam in the sea), hard on land
    db = ndi.distance_transform_edt(~beyond) / P * meters_per_point(40) / 1000
    soft = np.clip(db / 8, 0, 1)
    a = a * soft * soft * (3 - 2 * soft)
    a[beyond & ~land2] = 0
    a[land2] = 1
    a = np.round(a * 15) / 15                     # 16 levels: smooth enough, compresses well
    Image.fromarray(np.round(a * 255).astype(np.uint8)).save(os.path.join(WORK, 'detail_area.png'), optimize=True)
    # validation against the general map's coastline (the transition zoom is 2.6-3.8x)
    base = json.load(open(mapbase_path()))
    gen =[np.array(r).reshape(-1, 2) for r in base['geoms']['korea']]
    dist_osm = ndi.distance_transform_edt(~(np.array(Image.fromarray(np.zeros((1, 1), np.uint8))) > 0)) if False else None
    coastline = Image.new('L', (1360 * 4, 2400 * 4), 0)
    dr = ImageDraw.Draw(coastline)
    for r in [ring] + islands:
        dr.line([(x * 4, y * 4) for x, y in r.tolist()], fill=255, width=1)
    dco = ndi.distance_transform_edt(~(np.array(coastline) > 0)) / 4                  # map points to OSM coast
    bl = Image.new('L', (1360 * 4, 2400 * 4), 0)
    ImageDraw.Draw(bl).line([(x * 4, y * 4) for x, y in np.vstack([yalu, tumen]).tolist()], fill=255, width=1)
    dbl = ndi.distance_transform_edt(~(np.array(bl) > 0)) / 4                        # map points to the land border
    ds = []
    for r in gen:
        for x, y in r:
            i, j = int(y * 4), int(x * 4)
            if 0 <= i < dco.shape[0] and 0 <= j < dco.shape[1] and dbl[i, j] > 6:   # coast only, not the land border
                ds.append(dco[i, j])
    ds = np.array(ds)
    # a matching coastline for the general map (applied to the general-map scene by the app's own map build): the OSM coast from the
    # 압록강 mouth round to the 두만강 mouth, closed along the general map's own land border (Natural Earth), so the
    # border still meets its neighbouring-land polygon; plus the OSM islands. Simplified to 0.5 map pt.
    sys.path.insert(0, HERE)
    from geom import rdp
    ne = max((np.array(r).reshape(-1, 2) for r in base['geoms']['korea2']), key=len)
    iy = int(np.argmin(((ne - coast_part[0]) ** 2).sum(axis=1)))
    it = int(np.argmin(((ne - coast_part[-1]) ** 2).sum(axis=1)))
    arc1 = ne[it:iy + 1] if it <= iy else np.vstack([ne[it:], ne[:iy + 1]])
    arc2 = (ne[iy:it + 1] if iy <= it else np.vstack([ne[iy:], ne[:it + 1]]))[::-1]
    arc = arc1 if arc1[:, 1].max() < arc2[:, 1].max() else arc2       # the northern one (the land border)
    cp = rdp([tuple(q) for q in coast_part.tolist()], 0.5)
    gen_main = np.vstack([np.array(cp), arc])
    gen_rings = [[round(v, 1) for q in gen_main.tolist() for v in q]]
    for r in islands:
        if area(r) < 0.5:
            continue
        q = rdp([tuple(t) for t in r.tolist()], 0.5)
        if len(q) >= 4:
            gen_rings.append([round(v, 1) for t in q for v in t])
    os.makedirs(os.path.join(OUT, 'general'), exist_ok=True)
    json.dump({'note': 'Korean land for the general map: OSM coastline (simplified 0.5 pt) closed along the Natural '
                       'Earth land border, plus OSM islands. Rings are flat x, y lists in map points.',
               'rings': gen_rings},
              open(os.path.join(OUT, 'general', 'korea_land_general.json'), 'w'), separators=(',', ':'))
    gdist = []
    for r in gen_rings:
        a2 = np.array(r).reshape(-1, 2)
        for x, y in a2:
            i2, j2 = int(y * 4), int(x * 4)
            if 0 <= i2 < dco.shape[0] and 0 <= j2 < dco.shape[1] and dbl[i2, j2] > 6:
                gdist.append(dco[i2, j2])
    gdist = np.array(gdist)
    chk = {'general_coast_vertices': int(ds.size),
           'offset_map_pt': {'median': round(float(np.median(ds)), 2), 'p90': round(float(np.percentile(ds, 90)), 2),
                             'p99': round(float(np.percentile(ds, 99)), 2)},
           'offset_screen_pt_at_3x': {'median': round(float(np.median(ds)) * 3, 1), 'p90': round(float(np.percentile(ds, 90)) * 3, 1)},
           'islands': len(islands), 'mainland_points': int(len(ring)),
           'new_general_offset_map_pt': {'median': round(float(np.median(gdist)), 2), 'p90': round(float(np.percentile(gdist, 90)), 2),
                                         'p99': round(float(np.percentile(gdist, 99)), 2)}}
    json.dump(chk, open(os.path.join(OUT, 'land_check.json'), 'w'), indent=1)
    print(chk)


def chains_as_rings(chains):
    return []


if __name__ == '__main__':
    main()
