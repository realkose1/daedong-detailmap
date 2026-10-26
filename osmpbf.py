"""Minimal OpenStreetMap PBF reader (no extra packages: zlib + numpy).

Reads the blocks of a .osm.pbf file and yields what the detail-map build needs:
  iter_blocks(path)                     -> PrimitiveBlock bytes
  dense_nodes(block) / ways(block) / relations(block) in decoded form
Packed varint arrays are decoded with numpy, so a country extract takes seconds per pass, not hours.
Format reference: https://wiki.openstreetmap.org/wiki/PBF_Format
"""
import struct
import zlib
import numpy as np


def _varint(buf, i):
    shift = result = 0
    while True:
        b = buf[i]
        i += 1
        result |= (b & 0x7F) << shift
        if b < 0x80:
            return result, i
        shift += 7


def fields(buf, start=0, end=None):
    """Yield (field number, wire type, value) where value is an int (varint/fixed) or (start, end) for bytes."""
    i, end = start, len(buf) if end is None else end
    while i < end:
        key, i = _varint(buf, i)
        fn, wt = key >> 3, key & 7
        if wt == 0:
            v, i = _varint(buf, i)
            yield fn, wt, v
        elif wt == 2:
            ln, i = _varint(buf, i)
            yield fn, wt, (i, i + ln)
            i += ln
        elif wt == 1:
            yield fn, wt, struct.unpack_from('<q', buf, i)[0]
            i += 8
        elif wt == 5:
            yield fn, wt, struct.unpack_from('<i', buf, i)[0]
            i += 4
        else:
            raise ValueError('wire type %d' % wt)


def packed(buf, span, signed=False):
    """Packed varints -> int64 numpy array (zigzag-decoded when signed)."""
    a = np.frombuffer(buf, np.uint8, span[1] - span[0], span[0])
    if a.size == 0:
        return np.zeros(0, np.int64)
    ends = np.flatnonzero(a < 0x80)
    starts = np.empty_like(ends)
    starts[0] = 0
    starts[1:] = ends[:-1] + 1
    lens = ends - starts + 1
    grp = np.repeat(np.arange(ends.size), lens)
    pos = np.arange(a.size) - starts[grp]
    v = (a & 0x7F).astype(np.uint64) << (7 * pos).astype(np.uint64)
    out = np.add.reduceat(v, starts)
    if signed:
        return ((out >> np.uint64(1)).astype(np.int64)) ^ (-(out & np.uint64(1)).astype(np.int64))
    return out.astype(np.int64)


def iter_blocks(path):
    """Yield decompressed OSMData blocks."""
    with open(path, 'rb') as f:
        while True:
            hdr = f.read(4)
            if len(hdr) < 4:
                return
            n = struct.unpack('>I', hdr)[0]
            bh = f.read(n)
            btype, size = None, 0
            for fn, wt, v in fields(bh):
                if fn == 1:
                    btype = bh[v[0]:v[1]].decode()
                elif fn == 3:
                    size = v
            blob = f.read(size)
            if btype != 'OSMData':
                continue
            raw = None
            for fn, wt, v in fields(blob):
                if fn == 1:
                    raw = blob[v[0]:v[1]]
                elif fn == 3:
                    raw = zlib.decompress(blob[v[0]:v[1]])
            if raw is not None:
                yield raw


class Block:
    def __init__(self, raw):
        self.raw = raw
        self.strings = []
        self.groups = []
        self.gran, self.lat_off, self.lon_off = 100, 0, 0
        for fn, wt, v in fields(raw):
            if fn == 1:
                self.strings = [raw[a:b] for f2, w2, (a, b) in fields(raw, *v)]
            elif fn == 2:
                self.groups.append(v)
            elif fn == 17:
                self.gran = v
            elif fn == 19:
                self.lat_off = v
            elif fn == 20:
                self.lon_off = v
        self.s = [b.decode('utf-8', 'replace') for b in self.strings]

    def kinds(self):
        """Set of entity kinds in this block: 'dense', 'node', 'way', 'relation'."""
        k = set()
        for g in self.groups:
            for fn, wt, v in fields(self.raw, *g):
                k.add({1: 'node', 2: 'dense', 3: 'way', 4: 'relation'}.get(fn, 'other'))
                break
        return k

    def dense(self):
        """(ids, lat, lon, tags-per-node list or None) for every DenseNodes group."""
        raw = self.raw
        for g in self.groups:
            for fn, wt, v in fields(raw, *g):
                if fn != 2:
                    continue
                ids = lat = lon = kv = None
                for f2, w2, v2 in fields(raw, *v):
                    if f2 == 1:
                        ids = np.cumsum(packed(raw, v2, True))
                    elif f2 == 8:
                        lat = np.cumsum(packed(raw, v2, True))
                    elif f2 == 9:
                        lon = np.cumsum(packed(raw, v2, True))
                    elif f2 == 10:
                        kv = packed(raw, v2)
                la = 1e-9 * (self.lat_off + self.gran * lat)
                lo = 1e-9 * (self.lon_off + self.gran * lon)
                yield ids, la, lo, kv

    def dense_tags(self, kv, want_keys=None):
        """Split a DenseNodes keys_vals array into per-node (index, {k: v}) for tagged nodes only.
        want_keys: only nodes having one of these keys are returned (cheap pre-filter)."""
        keys = None if want_keys is None else {i for i, t in enumerate(self.s) if t in want_keys}
        out = []
        if kv is None or kv.size == 0:
            return out
        zeros = np.flatnonzero(kv == 0)
        starts = np.empty_like(zeros)
        starts[0] = 0
        starts[1:] = zeros[:-1] + 1
        tagged = np.flatnonzero(zeros > starts)
        kvl = kv.tolist()
        st, zl = starts.tolist(), zeros.tolist()
        s = self.s
        for node in tagged.tolist():
            seg = kvl[st[node]:zl[node]]
            if keys is not None and not any(seg[i] in keys for i in range(0, len(seg), 2)):
                continue
            out.append((node, {s[seg[i]]: s[seg[i + 1]] for i in range(0, len(seg) - 1, 2)}))
        return out

    def ways(self, want_refs):
        """Yield (id, tags, refs or None). want_refs(id, tags) decides whether refs are decoded."""
        raw = self.raw
        for g in self.groups:
            for fn, wt, v in fields(raw, *g):
                if fn != 3:
                    continue
                wid, keys, vals, refs_span = 0, None, None, None
                for f2, w2, v2 in fields(raw, *v):
                    if f2 == 1:
                        wid = v2
                    elif f2 == 2:
                        keys = packed(raw, v2)
                    elif f2 == 3:
                        vals = packed(raw, v2)
                    elif f2 == 8:
                        refs_span = v2
                tags = {self.s[k]: self.s[x] for k, x in zip(keys, vals)} if keys is not None else {}
                refs = None
                if refs_span is not None and want_refs(wid, tags):
                    refs = np.cumsum(packed(raw, refs_span, True))
                yield wid, tags, refs

    def relations(self):
        """Yield (id, tags, [(type, ref, role)])."""
        raw = self.raw
        for g in self.groups:
            for fn, wt, v in fields(raw, *g):
                if fn != 4:
                    continue
                rid, keys, vals, roles, mem, types = 0, None, None, None, None, None
                for f2, w2, v2 in fields(raw, *v):
                    if f2 == 1:
                        rid = v2
                    elif f2 == 2:
                        keys = packed(raw, v2)
                    elif f2 == 3:
                        vals = packed(raw, v2)
                    elif f2 == 8:
                        roles = packed(raw, v2)
                    elif f2 == 9:
                        mem = np.cumsum(packed(raw, v2, True))
                    elif f2 == 10:
                        types = packed(raw, v2)
                tags = {self.s[k]: self.s[x] for k, x in zip(keys, vals)} if keys is not None else {}
                members = []
                if mem is not None:
                    members = [(int(t), int(m), self.s[r]) for t, m, r in zip(types, mem, roles)]
                yield rid, tags, members
