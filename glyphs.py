"""Mountain glyphs on real summits, in the general map's brush vocabulary (option B).

Library used by bundle.py (data for the app).

Glyph items use the app's scene format: p (7 points, two cubic curves, the 4th point is the summit), f (famous),
x (준법 strokes, quadratic a-c-b), d (태점 x, y, r), and o (opacity while fading in; MapScene needs this one field).

Choices (all tunable here):
  * type from terrain around the summit: steep or rough -> sharp spire (금강·설악), middling -> rounded 토산,
    low or gentle -> low mound;
  * rocky massifs become clusters: nearby lower sharp summits join the most prominent one as extra spires;
  * a summit is fully in at zoom_full(prominence) and fades in (opacity and size) from zoom_full / 1.8;
  * screen height grows with zoom ** 0.42, with a cap per zoom, and is enlarged a little below 6x so the
    ink mass matches the general map's large hand-drawn mountains during the 2.5x-4x hand-over.
"""
import math

TEMPLATES = {
    'sharp': dict(aspect=0.72, p=[[-0.59, 1.0], [-0.39, 0.75], [-0.1, 0.18], [0, 0], [0.08, 0.22], [0.19, 0.7], [0.41, 1.0]],
                  x=[[(-0.025, 0.12), (-0.14, 0.55), (-0.17, 0.92)], [(-0.06, 0.22), (-0.25, 0.55), (-0.28, 0.92)],
                     [(0.03, 0.16), (0.08, 0.5), (0.1, 0.86)]],
                  d=[(0.04, 0.32, 1.7), (0.04, 0.17, 1.4), (0.17, 0.31, 1.4)]),
    'spire': dict(aspect=0.5, p=[[-0.55, 1.0], [-0.36, 0.72], [-0.09, 0.16], [0, 0], [0.08, 0.2], [0.2, 0.7], [0.45, 1.0]],
                  x=[[(-0.03, 0.15), (-0.13, 0.55), (-0.15, 0.9)]],
                  d=[(0.05, 0.3, 1.3)]),
    'round': dict(aspect=1.25, p=[[-0.5, 1.0], [-0.33, 0.6], [-0.18, 0.05], [0, 0], [0.18, 0.08], [0.33, 0.62], [0.5, 1.0]],
                  x=[[(-0.08, 0.18), (-0.16, 0.5), (-0.2, 0.85)]],
                  d=[(0.16, 0.30, 1.8)]),
    'mound': dict(aspect=2.0, p=[[-0.56, 1.0], [-0.4, 0.6], [-0.18, 0.05], [0, 0], [0.18, 0.08], [0.26, 0.62], [0.44, 1.0]],
                  x=[], d=[]),
}
KIND_ID = {'sharp': 0, 'round': 1, 'mound': 2, 'spire': 3}
REF = 1.5          # band unit -> screen point (the general map's glyphs read best around 1.5x)


def hash1(i, s):
    h = (i * 374761393 + s * 668265263) & 0xFFFFFFFF
    h = ((h ^ (h >> 13)) * 1274126177) & 0xFFFFFFFF
    h ^= h >> 16
    return h / 4294967295.0


def smoothstep(e0, e1, x):
    t = min(max((x - e0) / (e1 - e0), 0), 1)
    return t * t * (3 - 2 * t)


def kind(p):
    if (p['prom'] < 90 and p['ele'] < 600) or p['slope'] < 14:
        return 'mound'
    if p['slope'] >= 27 or p['rough'] >= 8.0:
        return 'sharp'
    return 'round'


def eff_prom(p):
    pr = p['prom']
    if p.get('name') and p.get('wd') and p['ele'] >= 1000:
        pr = max(pr, 160)
    return pr


def zoom_full(pr):
    """Fully in at: prominence 300 m -> 3.2x, 120 m -> 6.4x, 50 m -> 12.8x."""
    return 3.2 * (300.0 / max(pr, 1)) ** 0.774


def zoom_in(pr):
    return zoom_full(pr) / 1.8


def cap(zoom):
    return 58 if zoom < 6 else (64 if zoom < 12 else 72)


def screen_height(p, zoom, k=None):
    k = k or kind(p)
    pr = max(p['prom'], 30)
    h = min(56, max(11, 13 + 9.5 * math.log2(pr / 50)))
    h *= 1 + 0.18 * smoothstep(1000, 1600, p['ele'])
    h *= {'sharp': 1.1, 'round': 1.0, 'mound': 0.75, 'spire': 1.0}[k]
    h *= (zoom / 4) ** 0.42
    h *= 1 + 0.45 * (1 - smoothstep(3, 6, zoom))      # match the general map's ink mass at 2.5x-4x
    return min(h, cap(zoom))


def clusters(peaks):
    """Group nearby sharp summits into clusters around the most prominent one (head). Sets p['head'] / p['members']."""
    from bisect import bisect_left, bisect_right
    sharp = sorted([i for i, p in enumerate(peaks) if kind(p) == 'sharp'], key=lambda i: -peaks[i]['prom'])
    by_x = sorted(sharp, key=lambda i: peaks[i]['x'])
    xs = [peaks[i]['x'] for i in by_x]
    for i in sharp:
        p = peaks[i]
        if 'head' in p or p['prom'] < 60:
            continue
        r_km = min(1.6, 0.7 + p['prom'] / 1500)
        r = r_km * 1000 / 540.0                       # map points (close enough for grouping)
        lo, hi = bisect_left(xs, p['x'] - r), bisect_right(xs, p['x'] + r)
        mem = []
        for j in by_x[lo:hi]:
            q = peaks[j]
            if j == i or 'head' in q or 'members' in q or q['prom'] >= p['prom']:
                continue
            if math.hypot(q['x'] - p['x'], q['y'] - p['y']) <= r:
                mem.append(j)
        mem.sort(key=lambda j: math.hypot(peaks[j]['x'] - p['x'], peaks[j]['y'] - p['y']))
        mem = mem[:6]
        if mem:
            p['members'] = mem
            for j in mem:
                peaks[j]['head'] = i
    return peaks


def glyph_item(k, x, y, H, seed, famous, opacity):
    t = TEMPLATES[k]
    W = H * t['aspect'] * (0.88 + 0.24 * hash1(seed, 1))
    Hh = H * (0.92 + 0.16 * hash1(seed, 2))
    pts = []
    for n, (rx, ry) in enumerate(t['p']):
        jx = (hash1(seed, 10 + n) - 0.5) * 0.06 if n not in (0, 3, 6) else 0
        jy = (hash1(seed, 20 + n) - 0.5) * 0.06 if n not in (0, 3, 6) else 0
        pts += [round(x + (rx + jx) * W, 2), round(y + (ry + jy) * Hh, 2)]
    item = {'p': pts, 'f': 1 if famous else 0}
    if t['x'] and H >= 10:
        n = len(t['x']) if H >= 18 else 1
        item['x'] = [round(v, 2) for s in t['x'][:n] for rx, ry in s for v in (x + rx * W, y + ry * Hh)]
    if t['d'] and H >= 12:
        n = len(t['d']) if H >= 22 else 1
        item['d'] = [round(v, 2) for rx, ry, r in t['d'][:n] for v in (x + rx * W, y + ry * Hh, r)]
    if opacity < 0.999:
        item['o'] = round(opacity, 3)
    return item


def seed_of(p):
    return (int(p['x'] * 1000) * 31) ^ int(p['y'] * 1000)


def frame_glyphs(peaks, zoom, rect, margin=60, index=None):
    """Glyphs for one zoom inside rect (map points): (items in band units, labels, stats, K)."""
    K = zoom / REF
    x0, y0, x1, y1 = rect
    m = margin / zoom
    cand = index(x0 - m, y0 - m, x1 + m, y1 + 3 * m) if index else range(len(peaks))
    placed, chosen = [], []
    for i in sorted(cand, key=lambda i: -eff_prom(peaks[i])):
        p = peaks[i]
        if 'head' in p:
            continue                                   # drawn as part of its cluster
        if not (x0 - m <= p['x'] <= x1 + m and y0 - m <= p['y'] <= y1 + 3 * m):
            continue
        s = smoothstep(zoom_in(eff_prom(p)), zoom_full(eff_prom(p)), zoom)
        if s <= 0:
            continue
        k = kind(p)
        Hs = screen_height(p, zoom, k) * (0.55 + 0.45 * s)
        sx, sy = (p['x'] - x0) * zoom, (p['y'] - y0) * zoom
        if any(bx0 < sx < bx1 and by0 < sy < by1 for bx0, by0, bx1, by1 in placed):
            continue
        W = Hs * TEMPLATES[k]['aspect']
        placed.append((sx - W * 0.35, sy - 2, sx + W * 0.35, sy + Hs * 0.55))
        chosen.append((p, k, Hs, s, p['x'], p['y']))
        if k == 'sharp':
            # the cluster: real lower spires around the head, or two flanking ones on rough rock
            mem = [peaks[j] for j in p.get('members', [])]
            n = 0
            for q in mem:
                dxs = (q['x'] - p['x']) * zoom
                if abs(dxs) > 1.3 * W or abs((q['y'] - p['y']) * zoom) > 0.8 * Hs:
                    continue
                hq = Hs * min(0.85, 0.5 + 0.3 * math.sqrt(q['prom'] / max(p['prom'], 1)))
                chosen.append((q, 'spire', hq, s, q['x'], q['y']))
                n += 1
            if n < 2 and p['rough'] >= 8.0 and Hs >= 18:
                sd = seed_of(p)
                for side, f in ((-1, 0.62), (1, 0.52)):
                    ox = side * (0.28 + 0.08 * hash1(sd, 40 + side)) * W / zoom
                    oy = (0.10 + 0.06 * hash1(sd, 50 + side)) * Hs / zoom
                    chosen.append(({'x': p['x'] + ox, 'y': p['y'] + oy, 'prom': p['prom'] * f, 'ele': p['ele']},
                                   'spire', Hs * f, s, p['x'] + ox, p['y'] + oy))
    items, labels = [], []
    for p, k, Hs, s, x, y in sorted(chosen, key=lambda c: c[5] + c[2] / zoom):    # north first
        sd = (int(x * 1000) * 31) ^ int(y * 1000)
        famous = k != 'spire' and ((p['ele'] >= 1200 and p['prom'] >= 250) or (p.get('wd') and p['prom'] >= 150))
        items.append(glyph_item(k, x * K, y * K, Hs / REF, sd, famous, s))
        if p.get('name') and k != 'spire':
            labels.append({'kind': 'gpeak', 'name': p['name'], 'x': x, 'y': y, 'ele': p['ele'], 'wd': bool(p.get('wd')),
                           'area': Hs, 'prom': p['prom']})
    stats = {'glyphs': len(items), 'kinds': {kk: sum(1 for c in chosen if c[1] == kk) for kk in TEMPLATES}}
    return items, labels, stats, K


def ridge_wash(ridges, K, rect, zoom):
    """Faint connecting wash along the major divides: tapered filled ribbons (no line caps to show)."""
    x0, y0, x1, y1 = rect
    rings = []
    for l in ridges:
        pts = l['pts']
        xs, ys = pts[0::2], pts[1::2]
        if max(xs) < x0 - 10 or min(xs) > x1 + 10 or max(ys) < y0 - 10 or min(ys) > y1 + 10:
            continue
        w = l['w']
        if max(w) < 0.45:
            continue
        P = [(xs[i] * K, ys[i] * K) for i in range(len(xs))]
        n = len(P)
        if n < 3:
            continue
        s = [0.0]
        for i in range(1, n):
            s.append(s[-1] + math.hypot(P[i][0] - P[i - 1][0], P[i][1] - P[i - 1][1]))
        tot = s[-1] or 1
        left, right = [], []
        for i in range(n):
            a, b = P[max(0, i - 1)], P[min(n - 1, i + 1)]
            dx, dy = b[0] - a[0], b[1] - a[1]
            L = math.hypot(dx, dy) or 1e-6
            t = s[i] / tot
            env = smoothstep(0, 0.18, t) * smoothstep(0, 0.18, 1 - t)
            half = (8 + 10 * w[i]) * env / 2
            left.append((P[i][0] - dy / L * half, P[i][1] + dx / L * half))
            right.append((P[i][0] + dy / L * half, P[i][1] - dx / L * half))
        ring = left + right[::-1]
        rings.append([round(v, 2) for q in ring for v in q])
    return rings


def scene(items, ridges, K, rect, zoom):
    """Option-B glyph passes, with a finer wash buffer than the general map (no blocky bases on large glyphs)."""
    geoms = {'ridgewash': ridge_wash(ridges or [], K, rect, zoom)}
    ops = []
    if geoms['ridgewash']:
        ops += [{'t': 'wash', 'r': 4, 'cap': 6}, {'t': 'shift', 'dy': -2},
                {'t': 'fill', 'g': 'ridgewash', 'c': '#151515', 'a': 0.035},
                {'t': 'endShift'}, {'t': 'endWash'}]
    ops += [{'t': 'wash', 'r': 1.8, 'cap': 3.5},
            {'t': 'glyphs', 'items': items, 'm': 'wash', 'a': 0.48, 'depth': 0.6, 'tex': 0.14}, {'t': 'endWash'},
            {'t': 'glyphs', 'same': True, 'm': 'ink', 'occ': 0.92, 'top': 0.38, 'a': 0.82, 'tex': 0.3,
             'dot': 0.6, 'ds': 0.8, 'depth': 0.5}]
    return {'w': 1360 * K, 'h': 2400 * K, 'geoms': geoms, 'ops': ops}


def move_general(base, peaks, index=None, radius=16):
    """Move the general map's hand-drawn mountains onto real summits. Returns [(item index, dx, dy)]."""
    gl = next(o for o in base['ops'] if o['t'] == 'glyphs' and o.get('items'))
    items = gl['items']
    big = [i for i, p in enumerate(peaks) if p['prom'] >= 120]
    used = set()
    moves = []
    order = sorted(range(len(items)), key=lambda k: -(((items[k]['p'][1] + items[k]['p'][13]) / 2) - items[k]['p'][7]))
    for k in order:
        g = items[k]
        sx, sy = g['p'][6], g['p'][7]
        H = (g['p'][1] + g['p'][13]) / 2 - sy
        cand = index(sx - radius, sy - radius, sx + radius, sy + radius, big) if index else big
        best, bs = None, -1e9
        for j in cand:
            if j in used:
                continue
            p = peaks[j]
            d = math.hypot(p['x'] - sx, p['y'] - sy)
            if d > radius:
                continue
            score = math.log(p['prom']) * (1 + H / 30) - d / 4
            if score > bs:
                best, bs = j, score
        if best is None:
            continue
        used.add(best)
        dx, dy = peaks[best]['x'] - sx, peaks[best]['y'] - sy
        g['p'] = [round(v + (dx if i % 2 == 0 else dy), 2) for i, v in enumerate(g['p'])]
        if 'x' in g:
            g['x'] = [round(v + (dx if i % 2 == 0 else dy), 2) for i, v in enumerate(g['x'])]
        if 'd' in g:
            g['d'] = [round(v + (dx if i % 3 == 0 else dy if i % 3 == 1 else 0), 2) for i, v in enumerate(g['d'])]
        moves.append((k, round(dx, 2), round(dy, 2)))
    order2 = sorted(range(len(items)), key=lambda i: (items[i]['p'][1] + items[i]['p'][13]) / 2)
    gl['items'] = [items[i] for i in order2]
    return moves


class Grid:
    """Bucket index of points (map points) for fast rectangle queries."""
    def __init__(self, pts, cell=16.0):
        self.cell = cell
        self.b = {}
        for i, (x, y) in enumerate(pts):
            self.b.setdefault((int(x // cell), int(y // cell)), []).append(i)

    def __call__(self, x0, y0, x1, y1, subset=None):
        out = []
        c = self.cell
        for gx in range(int(x0 // c), int(x1 // c) + 1):
            for gy in range(int(y0 // c), int(y1 // c) + 1):
                out += self.b.get((gx, gy), [])
        if subset is not None:
            s = set(subset)
            out = [i for i in out if i in s]
        return out
