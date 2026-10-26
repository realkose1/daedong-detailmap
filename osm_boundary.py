"""National territory of South and North Korea from the extracts (admin_level 2 boundary relations).

usage: python -I osm_boundary.py <out.pkl> <extract.osm.pbf> [<extract.osm.pbf> ...]
Writes {'KR': [rings], 'KP': [rings]} in app map points (outer rings; the polygons include territorial waters).
The land border follows the boundary rivers (압록강, 두만강), so it is also the line the detail layer stops at.
"""
import pickle, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from osmpbf import iter_blocks, Block  # noqa: E402
from osm_extract import project, assemble  # noqa: E402


def main():
    out, files = sys.argv[1], sys.argv[2:]
    rels = {}
    for f in files:
        for raw in iter_blocks(f):
            b = Block(raw)
            if 'relation' not in b.kinds():
                continue
            for rid, tags, members in b.relations():
                if tags.get('boundary') == 'administrative' and tags.get('admin_level') == '2':
                    iso = tags.get('ISO3166-1') or tags.get('ISO3166-1:alpha2')
                    if iso in ('KR', 'KP'):
                        rels[iso] = [(m, role) for t, m, role in members if t == 1]
    want = {m for v in rels.values() for m, r in v}
    ways = {}
    for f in files:
        for raw in iter_blocks(f):
            b = Block(raw)
            if 'way' not in b.kinds():
                continue
            for wid, tags, refs in b.ways(lambda wid, t: wid in want):
                if refs is not None:
                    ways[wid] = refs
    need = np.unique(np.concatenate(list(ways.values())))
    lat = np.full(need.size, np.nan)
    lon = np.full(need.size, np.nan)
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
    xy = project(lon, lat)
    res = {}
    for iso, members in rels.items():
        seqs = [ways[m] for m, role in members if role in ('outer', '') and m in ways]
        rings = assemble(seqs)
        res[iso] = [xy[np.searchsorted(need, r)] for r in rings]
        print(iso, 'members', len(members), 'outer rings', len(rings), 'points', sum(len(r) for r in rings))
    pickle.dump(res, open(out, 'wb'))


if __name__ == '__main__':
    main()
