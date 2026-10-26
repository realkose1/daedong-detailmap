"""Pull what the detail map needs out of the Geofabrik extracts (South and North Korea).

usage: python -I osm_extract.py <out.pkl> <extract.osm.pbf> [<extract.osm.pbf> ...]

Three passes per file (relations, ways, nodes). Keeps, in app map points:
  coast    natural=coastline ways
  water    waterway=river | stream | canal lines (with name), kind kept
  areas    natural=water / waterway=riverbank polygons (closed ways and multipolygon relations), river flag, name
  pois     named natural=peak|saddle, place=*, buddhist place_of_worship, waterway=waterfall
Features shared by both extracts (the border rivers, the coast near the DMZ) are kept once (by OSM id).
"""
import math
import pickle
import sys
import time
import numpy as np

sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from osmpbf import iter_blocks, Block  # noqa: E402
from common import S, TX, TY  # noqa: E402

POI_KEYS = {'natural', 'place', 'amenity', 'waterway'}
WATER_TAGS = lambda t: t.get('natural') == 'water' or t.get('waterway') == 'riverbank'


def project(lon, lat):
    x = S * np.radians(lon) + TX
    y = TY - S * np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))
    return np.stack([x, y], axis=1).astype(np.float64)


def assemble(seqs):
    """Join node-id sequences into closed rings. seqs: list of int64 arrays."""
    seqs = [list(s) for s in seqs if len(s) >= 2]
    rings = []
    while seqs:
        cur = seqs.pop()
        progress = True
        while cur[0] != cur[-1] and progress:
            progress = False
            for k, s in enumerate(seqs):
                if s[0] == cur[-1]:
                    cur += s[1:]
                elif s[-1] == cur[-1]:
                    cur += s[::-1][1:]
                elif s[-1] == cur[0]:
                    cur = s[:-1] + cur
                elif s[0] == cur[0]:
                    cur = s[::-1][:-1] + cur
                else:
                    continue
                seqs.pop(k)
                progress = True
                break
        if cur[0] == cur[-1] and len(cur) >= 4:
            rings.append(np.array(cur, dtype=np.int64))
    return rings


def main():
    out, files = sys.argv[1], sys.argv[2:]
    rels = {}
    member_ways = set()
    t0 = time.time()
    for f in files:
        for raw in iter_blocks(f):
            b = Block(raw)
            if 'relation' not in b.kinds():
                continue
            for rid, tags, members in b.relations():
                if tags.get('type') == 'multipolygon' and WATER_TAGS(tags):
                    rels[rid] = (tags, [(m, role) for t, m, role in members if t == 1])
                    member_ways.update(m for t, m, role in members if t == 1)
    print('relations', len(rels), 'member ways', len(member_ways), round(time.time() - t0, 1), 's')

    def keep_way(wid, tags):
        return (tags.get('natural') == 'coastline' or tags.get('waterway') in ('river', 'stream', 'canal')
                or WATER_TAGS(tags) or wid in member_ways
                or (tags.get('amenity') == 'place_of_worship' and tags.get('religion') == 'buddhist' and 'name' in tags))
    ways = {}
    for f in files:
        for raw in iter_blocks(f):
            b = Block(raw)
            if 'way' not in b.kinds():
                continue
            for wid, tags, refs in b.ways(keep_way):
                if refs is not None and wid not in ways:
                    ways[wid] = (tags, refs)
    print('ways kept', len(ways), round(time.time() - t0, 1), 's')
    need = np.unique(np.concatenate([r for t, r in ways.values()]))
    lat = np.full(need.size, np.nan)
    lon = np.full(need.size, np.nan)
    pois = {}
    for f in files:
        for raw in iter_blocks(f):
            b = Block(raw)
            if 'dense' not in b.kinds():
                continue
            for ids, la, lo, kv in b.dense():
                idx = np.searchsorted(need, ids)
                idx[idx >= need.size] = 0
                hit = need[idx] == ids
                lat[idx[hit]] = la[hit]
                lon[idx[hit]] = lo[hit]
                for node, tags in b.dense_tags(kv, POI_KEYS):
                    if 'name' not in tags:
                        continue
                    kind = None
                    if tags.get('natural') in ('peak', 'saddle'):
                        kind = tags['natural']
                    elif tags.get('place') in ('city', 'town', 'village', 'hamlet', 'suburb', 'quarter', 'island', 'islet'):
                        kind = 'place'
                    elif tags.get('amenity') == 'place_of_worship' and tags.get('religion') == 'buddhist':
                        kind = 'temple'
                    elif tags.get('waterway') == 'waterfall':
                        kind = 'waterfall'
                    if kind:
                        nid = int(ids[node])
                        ele = tags.get('ele', '')
                        try:
                            ele = float(ele.replace(',', '').replace('m', '').strip())
                        except ValueError:
                            ele = None
                        pois[nid] = {'kind': kind, 'name': tags.get('name:ko') or tags['name'], 'place': tags.get('place'),
                                     'ele': ele, 'wd': bool(tags.get('wikidata') or tags.get('wikipedia')),
                                     'lon': float(lo[node]), 'lat': float(la[node])}
    print('nodes resolved', int(np.isfinite(lat).sum()), 'of', need.size, 'pois', len(pois), round(time.time() - t0, 1), 's')
    xy_all = project(lon, lat)

    def coords(refs):
        k = np.searchsorted(need, refs)
        p = xy_all[k]
        return p[np.isfinite(p).all(axis=1)]
    feats = {'coast': [], 'water': [], 'areas': [], 'pois': []}
    for wid, (tags, refs) in ways.items():
        if tags.get('natural') == 'coastline':
            feats['coast'].append(coords(refs))
        if tags.get('waterway') in ('river', 'stream', 'canal'):
            feats['water'].append({'kind': tags['waterway'], 'name': tags.get('name:ko') or tags.get('name'), 'xy': coords(refs)})
        if WATER_TAGS(tags) and refs[0] == refs[-1] and len(refs) >= 4:
            feats['areas'].append({'river': tags.get('water') in ('river', 'stream', 'canal', 'oxbow') or tags.get('waterway') == 'riverbank',
                                   'name': tags.get('name:ko') or tags.get('name'), 'rings': [coords(refs)]})
        if tags.get('amenity') == 'place_of_worship':
            c = coords(refs).mean(axis=0)
            feats['pois'].append({'kind': 'temple', 'name': tags.get('name:ko') or tags['name'], 'x': float(c[0]), 'y': float(c[1]),
                                  'wd': bool(tags.get('wikidata') or tags.get('wikipedia'))})
    for rid, (tags, members) in rels.items():
        rings = []
        for role in ('outer', 'inner'):
            seqs = [ways[m][1] for m, r in members if (r or 'outer') == role and m in ways]
            for ring in assemble(seqs):
                rings.append((role, coords(ring)))
        outer = [r for role, r in rings if role == 'outer']
        inner = [r for role, r in rings if role == 'inner']
        for o in outer:
            feats['areas'].append({'river': tags.get('water') in ('river', 'stream', 'canal', 'oxbow') or tags.get('waterway') == 'riverbank',
                                   'name': tags.get('name:ko') or tags.get('name'), 'rings': [o] + inner})
    for nid, p in pois.items():
        xy = project(np.array([p['lon']]), np.array([p['lat']]))[0]
        q = dict(p, x=float(xy[0]), y=float(xy[1]))
        feats['pois'].append(q)
    pickle.dump(feats, open(out, 'wb'))
    print({k: len(v) for k, v in feats.items()}, round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
