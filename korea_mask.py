"""Raster mask of Korean land (the app's Natural Earth 1:10m Korea polygon, ~1.5 km buffer), 4 cells per map point.
usage: python -I korea_mask.py [<MapBase.json>] <out.npz>
Without the scene path, common.mapbase_path() picks it ($DETAILMAP_MAPBASE, ...).
Used to keep summits and ridges on the Korean side of the Yalu and Tumen (the DEM tiles also hold China)."""
import json, os, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFilter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import mapbase_path  # noqa: E402
P = 4
args = sys.argv[1:]
if len(args) not in (1, 2):
    raise SystemExit('usage: python -I korea_mask.py [<MapBase.json>] <out.npz>')
scene = args[0] if len(args) == 2 else mapbase_path()
if not os.path.isfile(scene):
    raise SystemExit(f'general-map scene not found: {scene}')
geo = json.load(open(scene))['geoms']['korea']
im = Image.new('L', (1360 * P, 2400 * P), 0)
dr = ImageDraw.Draw(im)
for ring in geo:
    pts = [(ring[i] * P, ring[i + 1] * P) for i in range(0, len(ring) - 1, 2)]
    if len(pts) > 2:
        dr.polygon(pts, fill=255)
im = im.filter(ImageFilter.MaxFilter(13))
np.savez_compressed(args[-1], mask=np.array(im) > 0, p=P)
print('korea mask', args[-1])
