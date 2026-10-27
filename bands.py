"""Mountain ranges as layered horizontal bands (oblique view): the data for the app's detail layer.

usage: python -I bands.py            whole country -> out/bundle/bands_{4,8,16}.bin  (windows cached in work/bands_full.pkl)
       python -I bands.py --pack     only pack again from the cache
       (window() computes the bands of a single window, handy for previews)

The land is cut into horizontal slabs on a fixed map grid (spacing per level: 16, 8, 4 map pt; a finer level keeps
the coarser level's lines and adds layers in between). For each slab the crest profile is the highest relief (height
above the valley floor, 8 km neighbourhood) met along each column of the slab, so ridges that cross the slab rise as
natural, irregular undulations. Per sample (on a fixed x grid, `STEP`):
  h  crest height above the line, map pt (exaggerated; tapers out on low hills)
  b  how far the body hangs below its line (map pt), so a crest stands near its real ridge
  w  ink weight 0..1: relief over the level's threshold x how much this stretch is the real crest here rather than
     the slope climbing to the next layer (dominance over the neighbour slabs); short scraps are removed
  k  crag share 0..1 (cliff cells over 42 degrees at the DEM's own 30 m within ~600 m, near the crest)
Heights in map pt shrink only slowly with the level (zm^-0.35), so a level change does not make the ranges jump.

File (geom.pack, magic DDB1), one per level; tiles are horizontal strips (i = 0, j = strip of STRIP map pt).
Payload: varint line count; per line: varint n (the line sits at y = n * SPACING), varint run count;
  per run: varint start (global sample index, x = index * STEP), varint count,
           count bytes h (x HQ), count bytes b (x HQ), count bytes w (x 255), count bytes k (x 255);
  then varint accent count, varint accents (global sample indices of notable summits on this line).
"""
import json, math, os, sys
from multiprocessing import Pool
import numpy as np
from scipy import ndimage as ndi
from rasterio.warp import Resampling
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import unproj, meters_per_point, dem_paths_for, WORK, OUT  # noqa: E402
from dem import load_dem, smoothstep  # noqa: E402
from geom import pack  # noqa: E402

D = 8                                                   # DEM cells per map point (~67 m)
LEVELS = (4, 8, 16)
SPACING = {4: 16.0, 8: 8.0, 16: 4.0}                    # map pt between layers
STEP = {4: 0.5, 8: 0.25, 16: 0.25}                      # map pt between samples along a line
SIG = {4: (380.0, 780.0), 8: (330.0, 700.0), 16: (220.0, 520.0)}   # relief (m): faint .. full weight
ACCENT = {4: 600, 8: 380, 16: 250}                      # prominence (m) of summits that get an accent
PROM = {4: (300.0, 800.0), 8: (220.0, 650.0), 16: (160.0, 500.0)}   # prominence (m): weight boost around a summit
MINLEN = {4: 6.0, 8: 3.5, 16: 2.0}                      # map pt: shortest stretch worth a brush stroke
HQ = {4: 8.0, 8: 12.0, 16: 16.0}                        # h and b quantum (units per map pt)
STRIP = {4: 128, 8: 64, 16: 32}                         # map pt per strip (file tile)
FLOOR = 8000.0                                          # m: valley-floor neighbourhood
BASE = 0.4                                              # share of the (smoothed) height the body hangs down
A = 52.0                                                # screen pt for a 1000 m rise at 5.66x
CORE = 128.0                                            # map pt per window (national build)
PAD = 56.0                                              # map pt of DEM around a window
XPAD = 30.0                                             # map pt of samples computed past the window's sides


def varint(v):
    out = bytearray()
    v = int(v)
    while True:
        b = v & 0x7F
        v >>= 7
        if v:
            out.append(b | 0x80)
        else:
            out.append(b)
            return bytes(out)


def href(L):
    """Height of a 1200 m rise at level L (map pt): the fixed yardstick for line weight in the app."""
    zm = L / math.sqrt(2)
    return A * 1.2 ** 0.9 * (zm / 5.66) ** 0.65 / zm


def vnoise(u, seed):
    """1-D value noise in 0..1 on integer lattice points of u (deterministic)."""
    i = np.floor(u).astype(np.int64)
    f = u - i

    def hsh(k):
        v = (k * 374761393 + seed * 668265263) & 0xFFFFFFFF
        v = ((v ^ (v >> 13)) * 1274126177) & 0xFFFFFFFF
        return (v & 0xFFFF) / 65535.0
    t = f * f * (3 - 2 * f)
    return hsh(i) * (1 - t) + hsh(i + 1) * t


def window(x0, y0, x1, y1, levels=LEVELS, peaks=None, korea=None):
    """Lines whose y lies in [y0, y1), samples whose x lies in [x0, x1) on the global grid.
    Returns {L: [(n, i0, h, b, w, k, acc), ...]} with numpy arrays starting at global sample index i0."""
    X0, Y0 = x0 - PAD, y0 - PAD
    w, h = int(round((x1 - x0 + 2 * PAD) * D)), int(round((y1 - y0 + 2 * PAD) * D))
    lo0, la1 = unproj(X0, Y0)
    lo1, la0 = unproj(X0 + w / D, Y0 + h / D)
    z = load_dem(X0, Y0, w, h, dem_paths_for(lo0, la0, lo1, la1), d=D, resampling=Resampling.average)
    valid = ~np.isnan(z)
    z = np.where(valid, z, 0).astype(np.float32)
    land = valid & (z > 0)
    if korea is not None:
        land &= korea(X0, Y0, w, h)
    if not land.any():
        return {L: [] for L in levels}
    px_m = meters_per_point((la0 + la1) / 2) / D
    floor = ndi.gaussian_filter(ndi.minimum_filter(z, size=int(FLOOR / px_m)), FLOOR / 2 / px_m)
    relief = np.clip(ndi.gaussian_filter(z, 0.7) - floor, 0, None)
    relief[~land] = 0
    # crags, at the DEM's own resolution, only near the core
    cx0, cy0 = x0 - XPAD, y0 - 26.0
    cw, ch = int(round((x1 - x0 + 2 * XPAD) * 2 * D)), int(round((y1 - y0 + 52.0) * 2 * D))
    c0, ca1 = unproj(cx0, cy0)
    c1, ca0 = unproj(cx0 + cw / (2 * D), cy0 + ch / (2 * D))
    z2 = np.nan_to_num(load_dem(cx0, cy0, cw, ch, dem_paths_for(c0, ca0, c1, ca1), d=2 * D,
                                resampling=Resampling.bilinear)).astype(np.float32)
    gy, gx = np.gradient(z2, px_m / 2)
    cliff = (np.degrees(np.arctan(np.hypot(gx, gy))) > 42).astype(np.float32)
    crag2 = smoothstep(0.08, 0.24, ndi.uniform_filter(cliff, int(600 / (px_m / 2))))
    crag = np.zeros_like(relief)
    r0c, c0c = int(round((cy0 - Y0) * D)), int(round((cx0 - X0) * D))
    small = crag2[:ch // 2 * 2, :cw // 2 * 2].reshape(ch // 2, 2, cw // 2, 2).max(axis=(1, 3))
    crag[r0c:r0c + small.shape[0], c0c:c0c + small.shape[1]] = small[:h - r0c, :w - c0c]
    out = {}
    for L in levels:
        S, step = SPACING[L], STEP[L]
        lo, hi = SIG[L]
        zm = L / math.sqrt(2)
        i_lo, i_hi = math.ceil(x0 / step), math.ceil(x1 / step)            # core samples [i_lo, i_hi)
        e_lo, e_hi = math.ceil((x0 - XPAD) / step), math.ceil((x1 + XPAD) / step)
        xs = np.arange(e_lo, e_hi) * step
        cols = np.clip(((xs - X0) * D).astype(int), 0, w - 1)
        band = max(1, int(step * D))

        def slab(n):
            r0, r1 = int(round((n * S - S / 2 - Y0) * D)), int(round((n * S + S / 2 - Y0) * D))
            r0, r1 = max(0, r0), min(h, r1)
            if r1 <= r0:
                return None
            rr = ndi.maximum_filter1d(relief[r0:r1, :], band, axis=1)
            top = rr.max(axis=0)
            high = rr >= 0.6 * top[None, :]          # crags count only up near the crest
            rk = ndi.maximum_filter1d(np.where(high, crag[r0:r1, :], 0), band, axis=1).max(axis=0)
            return top[cols], rk[cols]
        n0, n1 = math.ceil(y0 / S), math.ceil(y1 / S)                       # lines with y in [y0, y1)
        prof = {n: slab(n) for n in range(n0 - 1, n1 + 1)}
        lines = []
        for n in range(n0, n1):
            if prof[n] is None:
                continue
            p, rk = prof[n]
            if p.max() < lo * 0.55:
                continue
            nb = np.maximum(*(prof[m][0] if prof[m] is not None else np.zeros_like(p) for m in (n - 1, n + 1)))
            dom = smoothstep(0.62, 0.95, ndi.gaussian_filter1d(p / np.maximum(nb, 1.0), 3))
            # a rocky summit is narrow: widen it a little before smoothing so a dome like 백운대 keeps its crags
            rock = ndi.gaussian_filter1d(ndi.maximum_filter1d(rk, max(1, int(round(1.5 / step)))), 0.75 / step)
            ps = ndi.gaussian_filter1d(p, 1.4)
            # crags: pinnacles from ridged noise on the global x grid, only where it is rocky
            lam = 11.0 / L
            r1 = (1 - np.abs(2 * vnoise(xs / lam, 3) - 1)) ** 3
            r2 = (1 - np.abs(2 * vnoise(xs / (lam * 0.45), 7) - 1)) ** 3
            pr = np.maximum(ps * (1 + rock * (0.5 * (r1 + 0.25 * r2) - 0.2)), 0)
            sig = smoothstep(lo, hi, ps)
            shape = smoothstep(lo * 0.55, lo * 1.05, ps)
            wt = sig * (0.12 + 0.88 * dom ** 1.2)
            # prominent summits (북한산, 관악산, 무등산 ...) stand out from their own plain even when the massif is small
            for q in (peaks or ()):
                if q['prom'] >= PROM[L][0] and abs(q['y'] - n * S) <= S * 0.75:
                    boost = smoothstep(PROM[L][0], PROM[L][1], q['prom']) * 0.85
                    fall = np.exp(-0.5 * ((xs - q['x']) / 1.6) ** 2) * (1 - abs(q['y'] - n * S) / (S * 0.75)) ** 0.5
                    wt = np.maximum(wt, boost * fall * smoothstep(lo * 0.5, lo, ps))
            # no short scraps: weak stretches have to run for a while to be drawn; strong ones stay even when narrow
            kmin = max(3, int(round(MINLEN[L] / step)))
            on = wt > 0.05
            lab, nlab = ndi.label(on)
            strong = np.zeros(len(wt), bool)
            if nlab:
                mx = np.asarray(ndi.maximum(wt, lab, index=np.arange(1, nlab + 1)))
                strong[on] = mx[lab[on] - 1] >= 0.45
            wt = ndi.gaussian_filter1d(np.where(strong, wt, ndi.grey_opening(wt, size=kmin)), kmin / 5)
            ws = ndi.gaussian_filter1d(wt, kmin / 3)
            shape = shape * smoothstep(0.03, 0.3, ws)
            hs = A * (pr / 1000) ** 0.9 * (zm / 5.66) ** 0.65 * shape
            hm = hs / zm
            bs = BASE * ndi.gaussian_filter1d(hm, 6.0 / step)
            a, b = i_lo - e_lo, i_hi - e_lo
            hc = hm[a:b]
            if hc.max() <= 0.02:
                continue
            acc = []
            if peaks is not None:
                for q in peaks:
                    if q.get('name') and q['prom'] >= ACCENT[L] and abs(q['y'] - n * S) <= S / 2 and x0 <= q['x'] < x1:
                        gi = int(round(q['x'] / step))
                        if 0 <= gi - e_lo < len(wt) and wt[gi - e_lo] > 0.4:
                            acc.append(gi)
            lines.append((n, i_lo, hc.astype(np.float32), bs[a:b].astype(np.float32), wt[a:b].astype(np.float32),
                          rock[a:b].astype(np.float32), acc))
        out[L] = lines
    return out


def _korea():
    crop = np.load(os.path.join(WORK, 'korea_crop.npz'))
    m, p = crop['mask'], float(crop['p'])
    # a little room past the coast and the border rivers (the detail layer's area mask decides what shows)
    m = ndi.binary_dilation(m[::2, ::2], iterations=8)
    p = p / 2

    def korea(X0, Y0, w, h):
        rr = ((Y0 + (np.arange(h) + 0.5) / D) * p).astype(int)
        cc = ((X0 + (np.arange(w) + 0.5) / D) * p).astype(int)
        ok_r = (rr >= 0) & (rr < m.shape[0])
        ok_c = (cc >= 0) & (cc < m.shape[1])
        out = np.zeros((h, w), bool)
        sub = m[np.clip(rr, 0, m.shape[0] - 1)][:, np.clip(cc, 0, m.shape[1] - 1)]
        out[:] = sub & ok_r[:, None] & ok_c[None, :]
        return out
    return korea, m, p


_PEAKS = None
_KOREA = None


def _job(args):
    global _PEAKS, _KOREA
    if _PEAKS is None:
        allp = json.load(open(os.path.join(WORK, 'peaks.json')))
        _PEAKS = [q for q in allp if 'head' not in q and q['prom'] >= min(min(ACCENT.values()), min(v[0] for v in PROM.values()))]
        _KOREA = _korea()[0]
    x0, y0 = args
    peaks = [q for q in _PEAKS if x0 - XPAD <= q['x'] < x0 + CORE + XPAD and y0 - 20 <= q['y'] < y0 + CORE + 20]
    try:
        return args, window(x0, y0, x0 + CORE, y0 + CORE, peaks=peaks, korea=_KOREA)
    except Exception as e:                      # a window outside every DEM tile
        print('skip', args, e, flush=True)
        return args, {L: [] for L in LEVELS}


def main():
    korea, m, p = _korea()
    jobs = []
    for j in range(int(2400 // CORE) + 1):
        for i in range(int(1360 // CORE) + 1):
            x0, y0 = i * CORE, j * CORE
            sub = m[int(y0 * p):int((y0 + CORE) * p), int(x0 * p):int((x0 + CORE) * p)]
            if sub.size and sub.any():
                jobs.append((x0, y0))
    print('windows', len(jobs), flush=True)
    import pickle
    cache = os.path.join(WORK, 'bands_full.pkl')
    if '--pack' in sys.argv and os.path.exists(cache):
        full = pickle.load(open(cache, 'rb'))
        jobs = []
    else:
        full = {L: {} for L in LEVELS}
    with Pool(6) as pool:
        for k, (args, res) in enumerate(pool.imap_unordered(_job, jobs)):
            for L, lines in res.items():
                for (n, i0, hh, bb, ww, kk, acc) in lines:
                    d = full[L].setdefault(n, {'h': {}, 'acc': set()})
                    d['h'][i0] = (hh, bb, ww, kk)
                    d['acc'].update(acc)
            if k % 10 == 0:
                print(f'{k + 1}/{len(jobs)}', flush=True)
    if jobs:
        pickle.dump(full, open(cache, 'wb'))
    os.makedirs(os.path.join(OUT, 'bundle'), exist_ok=True)
    sizes = {}
    for L in LEVELS:
        step, S, q = STEP[L], SPACING[L], HQ[L]
        nx = int(math.ceil(1360 / step)) + 2
        strips = {}
        hmax = bmax = 0.0
        samples = 0
        for n in sorted(full[L]):
            H = np.zeros(nx, np.float32); B = np.zeros(nx, np.float32)
            Wt = np.zeros(nx, np.float32); K = np.zeros(nx, np.float32)
            for i0, (hh, bb, ww, kk) in full[L][n]['h'].items():
                H[i0:i0 + len(hh)] = hh; B[i0:i0 + len(hh)] = bb
                Wt[i0:i0 + len(hh)] = ww; K[i0:i0 + len(hh)] = kk
            on = H > 0.02
            on = on | np.r_[on[1:], False] | np.r_[False, on[:-1]]
            idx = np.flatnonzero(on)
            if len(idx) == 0:
                continue
            br = np.flatnonzero(np.diff(idx) > 1)
            starts, ends = np.r_[idx[0], idx[br + 1]], np.r_[idx[br], idx[-1]]
            runs = [(s, e) for s, e in zip(starts, ends) if e - s >= 3]
            if not runs:
                continue
            hmax, bmax = max(hmax, float(H.max())), max(bmax, float(B.max()))
            rec = bytearray(varint(n) + varint(len(runs)))
            for s, e in runs:
                c = e - s + 1
                samples += int(c)
                rec += varint(int(s)) + varint(int(c))
                for arr, scale in ((H, q), (B, q), (Wt, 255), (K, 255)):
                    rec += np.clip(np.round(arr[s:e + 1] * scale), 0, 255).astype(np.uint8).tobytes()
            acc = sorted(full[L][n]['acc'])
            rec += varint(len(acc)) + b''.join(varint(int(a)) for a in acc)
            j = int((n * S) // STRIP[L])
            strips.setdefault(j, []).append(bytes(rec))
        tiles = {(0, j): varint(len(recs)) + b''.join(recs) for j, recs in strips.items()}
        data = pack(tiles, magic=b'DDB1')
        open(os.path.join(OUT, 'bundle', f'bands_{L}.bin'), 'wb').write(data)
        sizes[f'bands_{L}.bin'] = len(data)
        sizes[f'bands_{L}'] = {'lines': len(full[L]), 'samples': samples, 'hmax': round(hmax, 2), 'bmax': round(bmax, 2),
                               'h_clipped': hmax * q > 255}
        print(L, sizes[f'bands_{L}.bin'], sizes[f'bands_{L}'], flush=True)
    path = os.path.join(OUT, 'sizes.json')
    try:
        allsizes = json.load(open(path))
    except (OSError, ValueError):
        allsizes = {}
    allsizes.update(sizes)
    json.dump(allsizes, open(path, 'w'), indent=1, ensure_ascii=False)


if __name__ == '__main__':
    main()
