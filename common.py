"""Shared constants: the app's map projection, preview frames, DEM tile paths.

The app projects with d3.geoMercator().fitExtent([[30,30],[1330,2370]], MultiPoint(124.0 32.85, 131.95 43.15))
(the app's general-map build). d3's spherical Mercator is EPSG:3857 up to scale and offset:
    x = S*lambda + TX,  y = TY - S*ln(tan(pi/4 + phi/2))      (lambda, phi in radians)
    x = (S/R)*X3857 + TX, y = TY - (S/R)*Y3857
"""
import math

R = 6378137.0
W, H, M = 1360.0, 2400.0, 30.0


def _merc(lat):
    return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))


_l0, _l1 = math.radians(124.0), math.radians(131.95)
_m0, _m1 = _merc(32.85), _merc(43.15)   # m0 < m1
S = min((W - 2 * M) / (_l1 - _l0), (H - 2 * M) / (_m1 - _m0))
# d3 fitExtent: x = x0 + (w - k*(bx1+bx0))/2 with bx = S*lambda; by = -S*merc
TX = M + ((W - 2 * M) - S * (_l0 + _l1)) / 2
TY = M + ((H - 2 * M) + S * (_m0 + _m1)) / 2


def proj(lon, lat):
    return S * math.radians(lon) + TX, TY - S * _merc(lat)


def unproj(x, y):
    lon = math.degrees((x - TX) / S)
    lat = math.degrees(2 * math.atan(math.exp((TY - y) / S)) - math.pi / 2)
    return lon, lat


def app_to_3857(x, y):
    return (x - TX) * R / S, (TY - y) * R / S


def meters_per_point(lat):
    return R / S * math.cos(math.radians(lat))


SCREEN_W, SCREEN_H, DEV = 466.0, 678.0, 3.0

# name: (center lon, center lat, zoom)
FRAMES = {
    'jiri_8': (127.70, 35.33, 8.0),
    'halla_8': (126.55, 33.37, 8.0),
    'geumgang_8': (128.10, 38.64, 8.0),
    'myohyang_8': (126.25, 40.00, 8.0),
    'seorak_2.5': (128.58, 38.25, 2.5),
    'seorak_4': (128.50, 38.15, 4.0),
    'seorak_8': (128.53, 38.15, 8.0),
    'seorak_16': (128.47, 38.15, 16.0),
    'odae_8': (128.57, 37.76, 8.0),
    'bukhan_16': (126.91, 37.60, 16.0),
}


def frame_rect(name):
    lon, lat, z = FRAMES[name]
    cx, cy = proj(lon, lat)
    w, h = SCREEN_W / z, SCREEN_H / z
    return cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2, z


if __name__ == '__main__':
    print('S', S, 'TX', TX, 'TY', TY)
    for n, ll, app in [('낙산사', (128.6276, 38.1246), (787, 1212)), ('월정사', (128.5925, 37.7317), (781, 1293)),
                       ('해인사', (128.0981, 35.8013), (700, 1687)), ('통도사', (129.0640, 35.4877), (858, 1750)),
                       ('법주사', (127.8331, 36.5420), (657, 1537)), ('만월대', (126.5550, 37.9890), (447, 1242))]:
        x, y = proj(*ll)
        print(f'{n}: computed ({x:.1f}, {y:.1f}) app {app}  diff ({x - app[0]:+.1f}, {y - app[1]:+.1f}) pt')
    for n in FRAMES:
        x0, y0, x1, y1, z = frame_rect(n)
        a, b = unproj(x0, y1), unproj(x1, y0)
        print(n, f'pt [{x0:.1f},{y0:.1f},{x1:.1f},{y1:.1f}]', f'lon {a[0]:.3f}..{b[0]:.3f} lat {a[1]:.3f}..{b[1]:.3f}')


import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
DATA = _os.path.join(HERE, 'data')
WORK = _os.path.join(HERE, 'work')
OUT = _os.path.join(HERE, 'out')


def dem_tiles():
    """{(lat, lon): path} of the downloaded Copernicus GLO-30 tiles."""
    d = _os.path.join(DATA, 'dem')
    out = {}
    for f in _os.listdir(d) if _os.path.isdir(d) else []:
        if f.startswith('Copernicus_DSM_COG_10_N') and f.endswith('.tif'):
            lat = int(f.split('_N')[1][:2])
            lon = int(f.split('_E')[1][:3])
            out[(lat, lon)] = _os.path.join(d, f)
    return out


def dem_paths_for(lon0, lat0, lon1, lat1):
    tiles = dem_tiles()
    return [p for (la, lo), p in tiles.items() if la < lat1 and la + 1 > lat0 and lo < lon1 and lo + 1 > lon0]


def mapbase_path():
    """Path of the general-map scene (MapBase.json) that korea_mask.py, land.py and bundle.py read.

    In order: $DETAILMAP_MAPBASE; ../../Daedong/Resources/Map/MapBase.json (the app's own copy, when this folder sits
    inside the app repo); data/MapBase.json. Exits with an explanation when no file is found."""
    env = _os.environ.get('DETAILMAP_MAPBASE')
    if env:
        if _os.path.isfile(env):
            return env
        raise SystemExit(f'DETAILMAP_MAPBASE points to a file that does not exist: {env}')
    for p in (_os.path.join(HERE, '..', '..', 'Daedong', 'Resources', 'Map', 'MapBase.json'),
              _os.path.join(DATA, 'MapBase.json')):
        if _os.path.isfile(p):
            return p
    raise SystemExit('general-map scene (MapBase.json) not found. Set DETAILMAP_MAPBASE=/path/to/MapBase.json '
                     '(or put the file at data/MapBase.json). What the scene must contain: see README.md.')
