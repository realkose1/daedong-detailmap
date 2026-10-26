"""Copernicus GLO-30 into the app's map-point grid (EPSG:3857 up to scale and offset, see common.py)."""
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import Affine
from common import S, TX, TY, R


def load_dem(x0, y0, w, h, paths, d=16, resampling=Resampling.bilinear):
    """Grid of w x h cells, d cells per map point, top-left at map point (x0, y0). nan where no data."""
    dst = np.full((h, w), np.nan, dtype=np.float32)
    px = R / (S * d)
    transform = Affine(px, 0, (x0 - TX) * R / S, 0, -px, (TY - y0) * R / S)
    for p in paths:
        with rasterio.open(p) as ds:
            tmp = np.full((h, w), np.nan, dtype=np.float32)
            reproject(source=rasterio.band(ds, 1), destination=tmp, src_transform=ds.transform, src_crs=ds.crs,
                      dst_transform=transform, dst_crs='EPSG:3857', resampling=resampling,
                      src_nodata=None, dst_nodata=np.nan)
            m = ~np.isnan(tmp)
            dst[m] = tmp[m]
    return dst


def smoothstep(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0, 1)
    return t * t * (3 - 2 * t)
