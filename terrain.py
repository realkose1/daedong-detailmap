"""Terrain from the DEM for one 1-degree tile (with a margin): summits with prominence, and the major
drainage divides (ridges) that carry the faint connecting wash between mountain glyphs.

usage: python -I terrain.py <lat> <lon> <korea_mask.npz> <out.json>
Works on a 67 m grid (8 cells per map point). Summits and ridges are kept only inside the tile proper, so
neighbouring tiles join without duplicates; the 0.2 degree margin gives the prominence and the catchments
room to reach past the tile edge.
"""
import json, math, os, sys
import numpy as np
from scipy import ndimage as ndi
from rasterio.warp import Resampling
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import proj, unproj, meters_per_point, dem_paths_for  # noqa: E402
from dem import load_dem, smoothstep  # noqa: E402

D = 8          # cells per map point
F = 1
MARGIN = 0.2   # degrees
MIN_PROM = 20.0


def d8(s, area):
    """Receivers (flat index) by steepest descent of surface s, and accumulation of `area` along them."""
    h, w = s.shape
    best = np.zeros_like(s)
    rec = np.arange(h * w, dtype=np.int64).reshape(h, w)
    idx = rec.copy()
    pad = np.pad(s, 1, mode='edge')
    ipad = np.pad(idx, 1, mode='edge')
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            nb = pad[1 + dy:1 + dy + h, 1 + dx:1 + dx + w]
            drop = (s - nb) / math.hypot(dx, dy)
            m = drop > best
            best[m] = drop[m]
            rec[m] = ipad[1 + dy:1 + dy + h, 1 + dx:1 + dx + w][m]
    rec = rec.ravel()
    order = np.argsort(-s.ravel(), kind='stable')
    acc = area.ravel().astype(np.float64).tolist()
    recl = rec.tolist()
    for i in order.tolist():
        r = recl[i]
        if r != i:
            acc[r] += acc[i]
    return rec, np.array(acc).reshape(h, w)


def chaikin(pts, n):
    for _ in range(n):
        if len(pts) < 3:
            return pts
        out = [pts[0]]
        for a, b in zip(pts, pts[1:]):
            out += [(0.75 * a[0] + 0.25 * b[0], 0.75 * a[1] + 0.25 * b[1]), (0.25 * a[0] + 0.75 * b[0], 0.25 * a[1] + 0.75 * b[1])]
        out.append(pts[-1])
        pts = out
    return pts


LEVELS = [(3000, 1.0), (900, 0.62), (250, 0.36), (70, 0.18)]   # catchments in 67 m cells (~13, 4, 1.1, 0.3 km2)
LEVEL_A = 0.62


def jump(ptr, n=24):
    for _ in range(n):
        ptr = ptr[ptr]
    return ptr


def smooth_chain(pts, k):
    """Moving average along the chain (end points fixed): removes the staircase of the cell-edge grid."""
    if k <= 0 or len(pts) < 2 * k + 2:
        return pts
    xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
    ker = np.ones(2 * k + 1) / (2 * k + 1)
    xp = np.convolve(np.pad(xs, k, mode='edge'), ker, mode='valid')
    yp = np.convolve(np.pad(ys, k, mode='edge'), ker, mode='valid')
    xp[0], yp[0], xp[-1], yp[-1] = xs[0], ys[0], xs[-1], ys[-1]
    return list(zip(xp.tolist(), yp.tolist()))


def rdp_idx(pts, eps):
    if len(pts) < 3:
        return list(range(len(pts)))
    keep = [False] * len(pts)
    keep[0] = keep[-1] = True
    st = [(0, len(pts) - 1)]
    while st:
        a, b = st.pop()
        (ax, ay), (bx, by) = pts[a], pts[b]
        dx, dy = bx - ax, by - ay
        L = math.hypot(dx, dy)
        md, mi = 0.0, -1
        for i in range(a + 1, b):
            d = math.hypot(pts[i][0] - ax, pts[i][1] - ay) if L < 1e-9 else abs(dy * (pts[i][0] - ax) - dx * (pts[i][1] - ay)) / L
            if d > md:
                md, mi = d, i
        if md > eps:
            keep[mi] = True
            st += [(a, mi), (mi, b)]
    return [i for i, k in enumerate(keep) if k]


def prominence(z, valid):
    h, w = z.shape
    zf = z.ravel()
    order = np.argsort(-zf, kind='stable')
    order = order[valid.ravel()[order]]
    parent = np.full(h * w, -1, dtype=np.int64)
    peak = np.full(h * w, -1, dtype=np.int64)
    prom = {}
    par = parent.tolist()
    pk = peak.tolist()
    zl = zf.tolist()

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    offs = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    for i in order.tolist():
        r, c = divmod(i, w)
        roots = set()
        for dr, dc in offs:
            rr, cc = r + dr, c + dc
            if 0 <= rr < h and 0 <= cc < w:
                j = rr * w + cc
                if par[j] != -1:
                    roots.add(find(j))
        if not roots:
            par[i] = i
            pk[i] = i
            continue
        roots = sorted(roots, key=lambda q: -zl[pk[q]])
        main = roots[0]
        par[i] = main
        for q in roots[1:]:
            prom[pk[q]] = zl[pk[q]] - zl[i]
            par[q] = main
    for i in order.tolist():
        if par[i] == i and pk[i] not in prom:
            prom[pk[i]] = zl[pk[i]] - min(zl)   # highest in its island / region
    return prom



def fill_depressions(z, land, eps=0.01):
    """Priority-flood: raise every pit to its spill height (plus a tiny slope) so water always reaches the sea
    or the window edge. Flat paddies and the DSM's flattened reservoirs otherwise break the drainage."""
    import heapq
    h, w = z.shape
    zf = z.astype(np.float64).ravel().copy()
    lf = land.ravel()
    done = ~lf.copy()
    edge = np.zeros((h, w), bool)
    edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
    edge |= land & ndi.binary_dilation(~land)
    heap = [(zf[i], i) for i in np.flatnonzero(edge.ravel() & lf).tolist()]
    heapq.heapify(heap)
    for _, i in heap:
        done[i] = True
    zl = zf.tolist()
    dn = done.tolist()
    offs = (-w - 1, -w, -w + 1, -1, 1, w - 1, w, w + 1)
    push, pop = heapq.heappush, heapq.heappop
    while heap:
        zc, i = pop(heap)
        c = i % w
        for o in offs:
            j = i + o
            if j < 0 or j >= h * w or dn[j]:
                continue
            cj = j % w
            if abs(cj - c) > 1:
                continue
            dn[j] = True
            zj = zl[j]
            if zj <= zc:
                zj = zc + eps
                zl[j] = zj
            push(heap, (zj, j))
    return np.array(zl, dtype=np.float64).reshape(h, w)


def trace_edges(ev, eh, minw):
    h = ev.shape[0]
    w = eh.shape[1]
    W1 = w + 1
    adj = {}

    def add(a, b, s):
        adj.setdefault(a, []).append((b, s))
        adj.setdefault(b, []).append((a, s))
    rr, cc = np.nonzero(ev >= minw)
    for r, c, sv in zip(rr.tolist(), cc.tolist(), ev[rr, cc].tolist()):
        add(r * W1 + c + 1, (r + 1) * W1 + c + 1, sv)
    rr, cc = np.nonzero(eh >= minw)
    for r, c, sv in zip(rr.tolist(), cc.tolist(), eh[rr, cc].tolist()):
        add((r + 1) * W1 + c, (r + 1) * W1 + c + 1, sv)
    used = set()
    chains = []

    def walk(a, b, s):
        nodes, ss = [a, b], [s]
        used.add((min(a, b), max(a, b)))
        prev, cur = a, b
        while len(adj[cur]) == 2:
            nb = [x for x in adj[cur] if x[0] != prev]
            if not nb:
                break
            nx, sv = nb[0]
            key = (min(cur, nx), max(cur, nx))
            if key in used:
                break
            used.add(key)
            nodes.append(nx)
            ss.append(sv)
            prev, cur = cur, nx
        return nodes, ss
    for node, nbs in adj.items():
        if len(nbs) != 2:
            for nb, sv in nbs:
                if (min(node, nb), max(node, nb)) not in used:
                    chains.append(walk(node, nb, sv))
    for node, nbs in adj.items():
        for nb, sv in nbs:
            if (min(node, nb), max(node, nb)) not in used:
                chains.append(walk(node, nb, sv))
    return chains, W1


def divides(rec, acc, landc):
    h, w = landc.shape
    n = h * w
    idx = np.arange(n)
    lf = landc.ravel()
    ev = np.zeros((h, w - 1), np.float32)
    eh = np.zeros((h - 1, w), np.float32)
    for T, wt in LEVELS[:2]:
        ch = (acc.ravel() >= T) & lf
        recc = np.where(ch[rec], rec, idx)
        nxt = np.where(ch, recc, idx)
        donors = np.bincount(nxt[ch & (nxt != idx)], minlength=n)
        stop = ch & ((donors >= 2) | (nxt == idx))
        ptr = np.where(stop | ~ch, idx, nxt)
        g = jump(ptr)
        key = np.where(ch, np.where(nxt != idx, g[nxt], idx), -1)
        ptr2 = np.where(ch | (rec == idx), idx, rec)
        first = jump(ptr2)
        lab = np.where(lf & ch[first], key[first], -1).reshape(h, w)
        near = ndi.binary_dilation(ch.reshape(h, w), iterations=1)
        a, b = lab[:, :-1], lab[:, 1:]
        mv = (a != b) & (a >= 0) & (b >= 0) & ~near[:, :-1] & ~near[:, 1:]
        a, b = lab[:-1, :], lab[1:, :]
        mh = (a != b) & (a >= 0) & (b >= 0) & ~near[:-1, :] & ~near[1:, :]
        ev = np.maximum(ev, mv * wt)
        eh = np.maximum(eh, mh * wt)
    return ev, eh


def main():
    lat, lon = int(sys.argv[1]), int(sys.argv[2])
    kmask = np.load(sys.argv[3])
    out = sys.argv[4]
    x0, y1 = proj(lon - MARGIN, lat - MARGIN)
    x1, y0 = proj(lon + 1 + MARGIN, lat + 1 + MARGIN)
    w, h = int((x1 - x0) * D), int((y1 - y0) * D)
    paths = dem_paths_for(lon - MARGIN, lat - MARGIN, lon + 1 + MARGIN, lat + 1 + MARGIN)
    z = load_dem(x0, y0, w, h, paths, d=D, resampling=Resampling.max)
    valid = ~np.isnan(z)
    z = np.where(valid, z, 0).astype(np.float32)
    land = valid & (z > 0)
    if land.sum() == 0:
        json.dump({'peaks': [], 'ridges': []}, open(out, 'w'))
        return
    mpp = meters_per_point(lat + 0.5)
    px_m = mpp / D
    G = lambda a, m: ndi.gaussian_filter(a, m / px_m)
    # summits
    prom = prominence(ndi.gaussian_filter(z, 0.7), land)
    gy, gx = np.gradient(G(z, 50), px_m)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    rad = max(1, int(round(400 / px_m)))
    sm = ndi.uniform_filter(slope, 2 * rad + 1)
    rough = np.sqrt(np.clip(ndi.uniform_filter(slope ** 2, 2 * rad + 1) - sm ** 2, 0, None))
    km = kmask['mask']
    kp = float(kmask['p'])

    def korean(x, y):
        r, c = int(y * kp), int(x * kp)
        return 0 <= r < km.shape[0] and 0 <= c < km.shape[1] and km[r, c]
    peaks = []
    for i, p in prom.items():
        if p < MIN_PROM:
            continue
        r, c = divmod(i, w)
        x, y = x0 + (c + 0.5) / D, y0 + (r + 0.5) / D
        lo, la = unproj(x, y)
        if not (lon <= lo < lon + 1 and lat <= la < lat + 1) or not korean(x, y):
            continue
        peaks.append({'x': round(x, 3), 'y': round(y, 3), 'ele': round(float(z[r, c])), 'prom': round(float(p)),
                      'slope': round(float(sm[r, c]), 1), 'rough': round(float(rough[r, c]), 1)})
    # major divides
    landc = land
    zfill = fill_depressions(ndi.gaussian_filter(z, 1.0), landc)
    rec, acc = d8(np.where(landc, zfill, -1.0), landc.astype(np.float64))
    ev, eh = divides(rec, acc, landc)
    mass = smoothstep(60, 1000, G(z, 700) - 0.6 * G(z, 8000))
    chains, W1 = trace_edges(ev, eh, LEVEL_A - 1e-3)
    ridges = []
    for nodes, ss in chains:
        if len(nodes) < 30:                     # ~2 km
            continue
        mid = nodes[len(nodes) // 2]
        mx, my = x0 + (mid % W1) / D, y0 + (mid // W1) / D
        lo, la = unproj(mx, my)
        if not (lon <= lo < lon + 1 and lat <= la < lat + 1) or not korean(mx, my):
            continue
        pts = [(x0 + (nd % W1) / D, y0 + (nd // W1) / D) for nd in nodes]
        wt = []
        for nd, sv in zip(nodes, [ss[0]] + ss):
            r, c = min(h - 1, nd // W1), min(w - 1, nd % W1)
            wt.append(sv * (0.3 + 0.7 * float(mass[r, c])))
        sp = chaikin(smooth_chain(pts, 4), 2)
        n0 = len(pts)
        ws = [wt[min(n0 - 1, int(round(k * (n0 - 1) / max(1, len(sp) - 1))))] for k in range(len(sp))]
        keep = rdp_idx(sp, 25 / mpp)
        ridges.append({'pts': [round(float(v), 3) for k in keep for v in sp[k]], 'w': [round(float(ws[k]), 2) for k in keep]})
    json.dump({'peaks': peaks, 'ridges': ridges}, open(out, 'w'))
    print(lat, lon, 'peaks', len(peaks), 'ridges', len(ridges))


if __name__ == '__main__':
    main()
