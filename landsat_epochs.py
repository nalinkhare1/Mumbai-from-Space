"""
Mumbai from Space, 1988-2026: download dry-season Landsat composites for each epoch.

For each epoch, takes up to N_DAYS least-cloudy Jan-Apr Landsat Collection 2 Level-2 scenes
(Landsat 5, 7, 8, 9) that cover all of Greater Mumbai, masks cloud/shadow/fill, and saves per-scene
MNDWI, median NDVI and a median true-colour image to epochs/<label>.npz.

pip install -r requirements.txt
python landsat_epochs.py
"""
import json
import os
from collections import defaultdict

import numpy as np
import pystac_client
import planetary_computer
from odc.stac import load
from shapely.geometry import box, shape

BOUNDARY = shape(json.load(open("greater_mumbai.geojson"))["geometry"])
lon0, lat0, lon1, lat1 = BOUNDARY.bounds
BBOX = (lon0 - 0.01, lat0 - 0.01, lon1 + 0.01, lat1 + 0.01)
CRS, RES = "EPSG:32643", 30

# label -> (first year, last year); early epochs pool a few years because scenes are sparse
EPOCHS = {
    "1988-91": (1988, 1991),
    "1993-95": (1993, 1995),
    "1997-99": (1997, 1999),
    "2001-03": (2001, 2003),
    "2005-07": (2005, 2007),
    "2009-11": (2009, 2011),
    "2014": (2014, 2014),
    "2018": (2018, 2018),
    "2022": (2022, 2022),
    "2026": (2026, 2026),
}
N_DAYS = 8
BANDS = ["green", "swir16", "red", "nir08", "blue", "qa_pixel"]
# qa_pixel bits: 0 fill, 1 dilated cloud, 3 cloud, 4 cloud shadow
QA_BAD = (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4)

cat = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",
                                modifier=planetary_computer.sign_inplace)


def pick_scenes(y0, y1):
    """Least-cloudy Jan-Apr scenes that cover the whole area, one per date."""
    items = cat.search(collections=["landsat-c2-l2"], bbox=BBOX,
                       datetime=f"{y0}-01-01/{y1}-04-30",
                       query={"eo:cloud_cover": {"lt": 25}}).item_collection()
    by_day = defaultdict(list)
    for it in items:
        if it.datetime.month <= 4 and shape(it.geometry).contains(box(*BBOX)):
            by_day[it.datetime.strftime("%Y-%m-%d")].append(it)
    best = [min(its, key=lambda i: i.properties["eo:cloud_cover"]) for its in by_day.values()]
    return sorted(best, key=lambda i: i.properties["eo:cloud_cover"])[:N_DAYS]


def reflectance(dn):
    # Collection 2 Level-2 surface reflectance scaling; 0 = no data
    return np.where(dn > 0, dn * 0.0000275 - 0.2, np.nan).astype("float32")


os.makedirs("epochs", exist_ok=True)
for label, (y0, y1) in EPOCHS.items():
    out = f"epochs/{label}.npz"
    if os.path.exists(out):
        print(f"{label}: cached")
        continue
    scenes = pick_scenes(y0, y1)
    mndwi, ndvi, rgb, used = [], [], [], []
    for it in scenes:
        ds = load([it], bands=BANDS, bbox=BBOX, crs=CRS, resolution=RES,
                  resampling={b: "nearest" if b == "qa_pixel" else "bilinear" for b in BANDS}
                  ).isel(time=0)
        f = {b: reflectance(ds[b].values) for b in BANDS if b != "qa_pixel"}
        ok = ((ds["qa_pixel"].values.astype("uint32") & QA_BAD) == 0)
        ok &= np.isfinite(f["green"]) & np.isfinite(f["swir16"])
        g, s = np.clip(f["green"], 1e-4, None), np.clip(f["swir16"], 1e-4, None)
        mndwi.append(np.where(ok, (g - s) / (g + s), np.nan).astype("float16"))
        nir, red = np.clip(f["nir08"], 1e-4, None), np.clip(f["red"], 1e-4, None)
        ndvi.append(np.where(ok, (nir - red) / (nir + red), np.nan).astype("float16"))
        rgb.append((it.properties["platform"] == "landsat-7",
                    np.where(ok[..., None], np.stack([f["red"], f["green"], f["blue"]], -1), np.nan)
                    .astype("float16")))
        used.append(f"{it.datetime:%Y-%m-%d} {it.properties['platform']} {ok.mean():.0%} clear")
        xs, ys = ds.x.values, ds.y.values
    # true colour from Landsat 5/8/9 when possible: Landsat 7 after 2003 has scan-line gaps (stripes)
    # (pixels that no Landsat 5/8/9 scene saw clearly fall back to all scenes)
    rgb_all = np.nanmedian(np.stack([r for _, r in rgb]).astype("float32"), 0)
    rgb_pick = [r for is7, r in rgb if not is7]
    rgb_out = (np.nanmedian(np.stack(rgb_pick).astype("float32"), 0) if len(rgb_pick) >= 2
               else rgb_all)
    rgb_out = np.where(np.isnan(rgb_out), rgb_all, rgb_out)
    np.savez_compressed(out, mndwi=np.stack(mndwi),
                        ndvi=np.nanmedian(np.stack(ndvi).astype("float32"), 0),
                        rgb=rgb_out,
                        xs=xs, ys=ys, scenes=np.array(used))
    print(f"{label}: {len(used)} scenes\n  " + "\n  ".join(used))
