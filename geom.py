"""Geometry helpers shared by the build: simplification, clipping to tiles, compact tile encoding."""
import math
import struct
import zlib

from common import meters_per_point, unproj

CLASS_ID = {'coast': 1, 'river': 2, 'stream': 3, 'water': 4, 'riverarea': 5, 'ridge': 6, 'land': 7}


def rdp(pts, eps):
    if len(pts) < 3:
        return pts
    keep = [False] * len(pts)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        ax, ay = pts[a]
        bx, by = pts[b]
        dx, dy = bx - ax, by - ay
        L = math.hypot(dx, dy)
        md, mi = 0.0, -1
        for i in range(a + 1, b):
            if L < 1e-9:  # closed ring: measure from the shared end point
                d = math.hypot(pts[i][0] - ax, pts[i][1] - ay)
            else:
                d = abs(dy * (pts[i][0] - ax) - dx * (pts[i][1] - ay)) / L
            if d > md:
                md, mi = d, i
        if md > eps:
            keep[mi] = True
            stack += [(a, mi), (mi, b)]
    return [p for p, k in zip(pts, keep) if k]


def clip_poly(poly, x0, y0, x1, y1):
    def inter(a, b, axis, v):
        t = (v - a[axis]) / (b[axis] - a[axis])
        return (v, a[1] + t * (b[1] - a[1])) if axis == 0 else (a[0] + t * (b[0] - a[0]), v)
    out = poly
    for axis, v, keep_ge in ((0, x0, True), (0, x1, False), (1, y0, True), (1, y1, False)):
        inp, out = out, []
        if not inp:
            break
        for i in range(len(inp)):
            c, p = inp[i], inp[i - 1]
            cin = c[axis] >= v if keep_ge else c[axis] <= v
            pin = p[axis] >= v if keep_ge else p[axis] <= v
            if cin:
                if not pin:
                    out.append(inter(p, c, axis, v))
                out.append(c)
            elif pin:
                out.append(inter(p, c, axis, v))
    return out


def clip_line(pts, x0, y0, x1, y1):
    """Liang-Barsky per segment; returns list of pieces."""
    pieces, cur = [], []
    for i in range(len(pts) - 1):
        (ax, ay), (bx, by) = pts[i], pts[i + 1]
        dx, dy = bx - ax, by - ay
        t0, t1, ok = 0.0, 1.0, True
        for p, q in ((-dx, ax - x0), (dx, x1 - ax), (-dy, ay - y0), (dy, y1 - ay)):
            if p == 0:
                if q < 0:
                    ok = False
                    break
            else:
                t = q / p
                if p < 0:
                    t0 = max(t0, t)
                else:
                    t1 = min(t1, t)
        if not ok or t0 > t1:
            if len(cur) > 1:
                pieces.append(cur)
            cur = []
            continue
        a = (ax + t0 * dx, ay + t0 * dy)
        b = (ax + t1 * dx, ay + t1 * dy)
        if not cur or cur[-1] != a:
            if len(cur) > 1:
                pieces.append(cur)
            cur = [a]
        cur.append(b)
        if t1 < 1.0:
            pieces.append(cur)
            cur = []
    if len(cur) > 1:
        pieces.append(cur)
    return pieces


def length_m(pts):
    s = 0.0
    for (ax, ay), (bx, by) in zip(pts, pts[1:]):
        s += math.hypot(bx - ax, by - ay) * meters_per_point(unproj(ax, (ay + by) / 2)[1])
    return s


def area_km2(ring):
    s = 0.0
    for i in range(len(ring)):
        (ax, ay), (bx, by) = ring[i - 1], ring[i]
        s += ax * by - bx * ay
    m = meters_per_point(unproj(ring[0][0], ring[0][1])[1])
    return abs(s) / 2 * m * m / 1e6


def centroid(ring):
    xs = [p[0] for p in ring]
    ys = [p[1] for p in ring]
    return sum(xs) / len(xs), sum(ys) / len(ys)


def varint(n, out):
    while True:
        b = n & 0x7F
        n >>= 7
        if n:
            out.append(b | 0x80)
        else:
            out.append(b)
            return


def zz(v):
    return (v << 1) ^ (v >> 63)


def encode_tile(feats, ox, oy, q):
    out = bytearray()
    varint(len(feats), out)
    for cls, parts in feats:
        varint(CLASS_ID[cls], out)
        varint(len(parts), out)
        for pts in parts:
            varint(len(pts), out)
            px, py = 0, 0
            for x, y in pts:
                ix, iy = int(round((x - ox) * q)), int(round((y - oy) * q))
                varint(zz(ix - px), out)
                varint(zz(iy - py), out)
                px, py = ix, iy
    return bytes(out)




def deflate(data):
    """Raw DEFLATE (no zlib header), which the app inflates with NSData.decompressed(using: .zlib)."""
    if not data:
        return b''
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    return c.compress(data) + c.flush()


def pack(tiles, magic=b'DDM1'):
    """Container: magic, tile count, then (i, j, offset, length) per tile (int16, int16, uint32, uint32), then the
    raw-deflated tile payloads. Tiles are addressed by integer tile coordinates (tile size depends on the file).
    An empty payload (length 0) is allowed and means "full" where a file says so (area tiles)."""
    keys = sorted(tiles)
    blobs = [deflate(tiles[k]) for k in keys]
    head = bytearray(magic) + struct.pack('<I', len(keys))
    off = len(head) + 12 * len(keys)
    for k, b in zip(keys, blobs):
        head += struct.pack('<hhII', k[0], k[1], off, len(b))
        off += len(b)
    return bytes(head) + b''.join(blobs)
