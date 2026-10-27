"""한양과 서라벌을 깊이 키웠을 때(8배 위, 끝까지 64배) 그릴 옛 자리와 큰길 (city layer for the deep zoom).

usage: python -I city.py [<south-korea-latest.osm.pbf>]     (기본: data/osm/south-korea-latest.osm.pbf)
writes out/bundle/city.bin  (앱에 넣을 때는 앱의 자료 폴더로 복사한다)

두 도시 둘레(CITIES의 위도·경도 상자)에서 OpenStreetMap의 다음 것을 뽑아 지도 좌표(포인트)로 옮긴다.
  wall   11  도성과 산성: barrier=city_wall 가운데 이름 있는 성벽(서울 한양도성 따위)이나 historic=citywalls와
             이어진 것 (끝 점을 함께 쓰는 조각까지 따라간다). 앱은 옛 지도의 성벽처럼 굵은 선과 톱니(점선)로 긋는다.
  fence  12  궁궐·사당의 담장: 그 밖의 barrier=city_wall, 그리고 historic 표시가 있는 barrier=wall.
  hall   13  옛 건물: building과 historic(또는 building:architecture=joseon)이 함께 있는 면 (전각, 문루, 정자).
  tomb   14  고분: historic=tomb 면 (서라벌의 봉분).
  site   15  유적의 터: historic=archaeological_site 면, 그리고 궁궐(historic=castle|palace) 면과 다각형 관계.
  road   16  큰길(trunk, primary), 17 버금길(secondary, tertiary). 도시 가운데에서 ROAD_KM 안만.
선과 면은 4포인트 조각으로 나누고(geom.pack, 'DDC1'), 점은 1/1024포인트(약 0.5m)로 적는다. 1.5m보다 작은 굴곡은 줄인다.
조각 꼴은 lines_*.bin과 같다 (bundle.py 머리 주석). 사용권: © OpenStreetMap 기여자 (ODbL 1.0).
"""
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from osmpbf import iter_blocks, Block  # noqa: E402
from common import proj, meters_per_point, OUT, DATA  # noqa: E402
from geom import rdp, clip_poly, varint, zz, pack  # noqa: E402
from osm_extract import assemble  # noqa: E402

# 이름: (서경, 남위, 동경, 북위), 큰길을 그릴 가운데(경도, 위도)와 반지름(km)
CITIES = {
    'seoul': ((126.93, 37.52, 127.04, 37.61), (126.984, 37.572), 4.2),
    'gyeongju': ((129.16, 35.74, 129.38, 35.88), (129.222, 35.835), 4.5),
}
CLASS = {'wall': 11, 'fence': 12, 'hall': 13, 'tomb': 14, 'site': 15, 'road': 16, 'road2': 17}
AREAS = {'hall', 'tomb', 'site'}
T = 4.0          # 조각 한 변 (지도 포인트)
Q = 1024.0       # 1포인트를 몇 눈금으로
TOL_M = 1.5
ROADS = {'trunk': 'road', 'primary': 'road', 'secondary': 'road2', 'tertiary': 'road2'}
HALL_HISTORIC = {'yes', 'building', 'castle', 'city_gate', 'shrine', 'temple', 'palace', 'manor', 'heritage_building',
                 'house', 'gate', 'monument', 'technical_monument', 'pavilion', 'church'}


def in_any_box(lo, la):
    m = np.zeros(lo.shape, bool)
    for (x0, y0, x1, y1), _, _ in CITIES.values():
        m |= (lo >= x0) & (lo <= x1) & (la >= y0) & (la <= y1)
    return m


def kind_of(tags):
    """선이나 면의 갈래 (도성/담장 갈래는 나중에 이어짐으로 다시 가른다). 없으면 None."""
    if tags.get('barrier') == 'city_wall' or tags.get('historic') == 'citywalls':
        return 'citywall'
    if tags.get('barrier') == 'wall' and 'historic' in tags:
        return 'fence'
    if 'building' in tags and (tags.get('historic') in HALL_HISTORIC or tags.get('building:architecture') == 'joseon'):
        return 'hall'
    if tags.get('historic') == 'tomb':
        return 'tomb'
    if tags.get('historic') == 'archaeological_site' or (tags.get('historic') in ('castle', 'palace') and 'building' not in tags):
        return 'site'
    if tags.get('highway') in ROADS:
        return ROADS[tags['highway']]
    return None


def split_tiles(pts):
    """선을 조각 경계에서 자른다: {(i, j): [조각, ...]}."""
    out = {}
    for a, b in zip(pts, pts[1:]):
        # 한 마디를 조각 경계마다 나눈다 (조각이 작아 마디가 여러 조각을 지날 수 있다).
        ts = [0.0, 1.0]
        for axis in (0, 1):
            lo, hi = sorted((a[axis], b[axis]))
            k = math.floor(lo / T) + 1
            while k * T < hi:
                ts.append((k * T - a[axis]) / (b[axis] - a[axis]))
                k += 1
        ts.sort()
        for t0, t1 in zip(ts, ts[1:]):
            if t1 - t0 < 1e-9:
                continue
            p = (a[0] + t0 * (b[0] - a[0]), a[1] + t0 * (b[1] - a[1]))
            q = (a[0] + t1 * (b[0] - a[0]), a[1] + t1 * (b[1] - a[1]))
            mid = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
            key = (int(mid[0] // T), int(mid[1] // T))
            pieces = out.setdefault(key, [])
            if pieces and pieces[-1][-1] == p:
                pieces[-1].append(q)
            else:
                pieces.append([p, q])
    return out


def encode(feats, ox, oy):
    out = bytearray()
    varint(len(feats), out)
    for cls, parts in feats:
        varint(CLASS[cls], out)
        varint(len(parts), out)
        for pts in parts:
            varint(len(pts), out)
            px = py = 0
            for x, y in pts:
                ix, iy = int(round((x - ox) * Q)), int(round((y - oy) * Q))
                varint(zz(ix - px), out)
                varint(zz(iy - py), out)
                px, py = ix, iy
    return bytes(out)


def main():
    pbf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, 'osm', 'south-korea-latest.osm.pbf')
    t0 = time.time()
    # 1. 다각형 관계 (궁궐, 유적): 바깥 고리를 이룰 길 번호
    rels = {}
    for raw in iter_blocks(pbf):
        b = Block(raw)
        if 'relation' not in b.kinds():
            continue
        for rid, tags, members in b.relations():
            if tags.get('type') == 'multipolygon' and kind_of(tags) in ('site', 'hall', 'tomb'):
                rels[rid] = (kind_of(tags), [m for t, m, role in members if t == 1 and role != 'inner'], tags.get('name', ''))
    member = {m for _, ms, _ in rels.values() for m in ms}
    print('관계', len(rels), round(time.time() - t0, 1), 's')
    # 2. 두 도시 상자 안의 점
    ids_l, lat_l, lon_l = [], [], []
    for raw in iter_blocks(pbf):
        b = Block(raw)
        if 'dense' not in b.kinds():
            continue
        for ids, la, lo, kv in b.dense():
            m = in_any_box(lo, la)
            if m.any():
                ids_l.append(ids[m]); lat_l.append(la[m]); lon_l.append(lo[m])
    ids = np.concatenate(ids_l)
    order = np.argsort(ids)
    ids, lat, lon = ids[order], np.concatenate(lat_l)[order], np.concatenate(lon_l)[order]
    print('점', ids.size, round(time.time() - t0, 1), 's')
    # 3. 길
    ways = {}
    for raw in iter_blocks(pbf):
        b = Block(raw)
        if 'way' not in b.kinds():
            continue
        for wid, tags, refs in b.ways(lambda w, t: kind_of(t) is not None or w in member):
            if refs is None:
                continue
            idx = np.searchsorted(ids, refs)
            idx[idx >= ids.size] = 0
            if not (ids[idx] == refs).all():
                continue
            ways[wid] = (kind_of(tags), tags, refs)
    print('길', len(ways), round(time.time() - t0, 1), 's')

    def xy(refs):
        idx = np.searchsorted(ids, refs)
        return [proj(float(lon[k]), float(lat[k])) for k in idx]

    def simp(pts):
        if len(pts) < 3:
            return pts
        mpp = meters_per_point(37.0 if pts[0][0] < 700 else 35.8)
        return rdp(pts, TOL_M / mpp)

    # 4. 성벽: 이름 있는 성벽(또는 historic=citywalls)과 끝 점으로 이어진 조각은 도성, 나머지는 담장
    walls = {w: v for w, v in ways.items() if v[0] == 'citywall'}
    by_node = {}
    for w, (_, tags, refs) in walls.items():
        for n in (int(refs[0]), int(refs[-1])):
            by_node.setdefault(n, []).append(w)
    seeds = [w for w, (_, tags, _) in walls.items() if 'name' in tags or tags.get('historic') == 'citywalls']
    city = set(seeds)
    stack = list(seeds)
    while stack:
        w = stack.pop()
        refs = walls[w][2]
        for n in (int(refs[0]), int(refs[-1])):
            for o in by_node.get(n, []):
                if o not in city:
                    city.add(o)
                    stack.append(o)

    feats = []   # (갈래, [고리 또는 선])
    for w, (k, tags, refs) in ways.items():
        if k is None:
            continue
        if k == 'citywall':
            k = 'wall' if w in city else 'fence'
        pts = xy(refs)
        if k in AREAS:
            if refs[0] != refs[-1] or len(refs) < 4:
                continue
            feats.append((k, [simp(pts)]))
        else:
            if k.startswith('road'):
                _, (clon, clat), km = CITIES['seoul' if pts[0][0] < 700 else 'gyeongju']
                cx, cy = proj(clon, clat)
                reach = km * 1000 / meters_per_point(clat)
                if min(math.hypot(x - cx, y - cy) for x, y in pts) > reach:
                    continue
            feats.append((k, [simp(pts)]))
    for rid, (k, ms, name) in rels.items():
        seqs = [ways[m][2] for m in ms if m in ways]
        if not seqs:
            continue
        for ring in assemble(seqs):
            feats.append((k, [simp(xy(ring))]))
    count = {}
    for k, _ in feats:
        count[k] = count.get(k, 0) + 1
    print('갈래별', count)

    tiles = {}
    for k, parts in feats:
        for pts in parts:
            if len(pts) < 2:
                continue
            if k in AREAS:
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                i0, i1, j0, j1 = int(min(xs) // T), int(max(xs) // T), int(min(ys) // T), int(max(ys) // T)
                for j in range(j0, j1 + 1):
                    for i in range(i0, i1 + 1):
                        if i0 == i1 and j0 == j1:
                            c = pts
                        else:
                            c = clip_poly(pts, i * T, j * T, i * T + T, j * T + T)
                        if len(c) >= 3:
                            tiles.setdefault((i, j), []).append((k, [c]))
            else:
                for key, pieces in split_tiles(pts).items():
                    tiles.setdefault(key, []).append((k, pieces))
    payload = {key: encode(v, key[0] * T, key[1] * T) for key, v in tiles.items()}
    data = pack(payload, b'DDC1')
    os.makedirs(os.path.join(OUT, 'bundle'), exist_ok=True)
    path = os.path.join(OUT, 'bundle', 'city.bin')
    open(path, 'wb').write(data)
    print(path, len(data), '바이트,', len(payload), '조각', round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
