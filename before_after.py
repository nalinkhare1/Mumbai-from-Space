"""
Mumbai from Space: quick visual check.

Pulls dry-season (Jan-Apr) Sentinel-2 imagery of the Coastal Road (South) stretch
for 2018 (before reclamation began) and 2026, and saves them side by side.

pip install -r requirements.txt
python before_after.py
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import pystac_client
import planetary_computer
from odc.stac import load
from shapely.geometry import box, shape

# Marine Drive -> Worli Sea Link (lon_min, lat_min, lon_max, lat_max)
BBOX = (72.795, 18.940, 72.835, 19.030)
PERIODS = {"2018": "2018-01-01/2018-04-30", "2026": "2026-01-01/2026-04-30"}
N_SCENES = 5  # least-cloudy scenes per period, median-composited

# Tried in order; falls through if a catalogue has no scenes or fails to load.
CATALOGUES = [
    ("Planetary Computer", "https://planetarycomputer.microsoft.com/api/stac/v1",
     "sentinel-2-l2a", planetary_computer.sign_inplace, ["B04", "B03", "B02"]),
    ("Earth Search", "https://earth-search.aws.element84.com/v1",
     "sentinel-2-l2a", None, ["red", "green", "blue"]),
]


def composite(period):
    for name, url, coll, modifier, bands in CATALOGUES:
        try:
            cat = pystac_client.Client.open(url, modifier=modifier)
            items = list(cat.search(collections=[coll], bbox=BBOX, datetime=period,
                                    query={"eo:cloud_cover": {"lt": 20}}).items())
            # keep only scenes whose footprint covers the whole area (Mumbai sits on a tile edge)
            items = [i for i in items if shape(i.geometry).contains(box(*BBOX))]
            if not items:
                print(f"  {name}: no scenes")
                continue
            items = sorted(items, key=lambda i: i.properties["eo:cloud_cover"])[:N_SCENES]
            ds = load(items, bands=bands, bbox=BBOX, crs="EPSG:32643",
                      resolution=10, groupby="solar_day")
            rgb = np.stack([ds[b].where(ds[b] > 0).median("time").values for b in bands], -1)
            dates = sorted({str(t)[:10] for t in ds.time.values})
            print(f"  {name}: {len(dates)} scenes {dates}")
            return rgb, name, dates
        except Exception as e:
            print(f"  {name}: failed ({e.__class__.__name__}: {e})")
    raise RuntimeError(f"No usable imagery for {period}")


def stretch(rgb):
    # per-image 2-98% stretch; also absorbs the 2022 Sentinel-2 radiometric offset
    lo, hi = np.nanpercentile(rgb, [2, 98])
    return np.clip(np.nan_to_num((rgb - lo) / (hi - lo)), 0, 1)


fig, axes = plt.subplots(1, 2, figsize=(10, 10))
for ax, (year, period) in zip(axes, PERIODS.items()):
    print(f"{year}:")
    rgb, src, dates = composite(period)
    ax.imshow(stretch(rgb))
    ax.set_title(f"{year} (Jan-Apr, {len(dates)}-scene median)", fontsize=13)
    ax.axis("off")
fig.suptitle("Mumbai Coastal Road (South) - Sentinel-2, 10 m", fontsize=15)
fig.text(0.5, 0.02, "Contains modified Copernicus Sentinel data", ha="center", fontsize=9)
plt.tight_layout(rect=(0, 0.03, 1, 0.97))
os.makedirs("images", exist_ok=True)
plt.savefig("images/coastal_road_2018_vs_2026.png", dpi=150)
print("Saved images/coastal_road_2018_vs_2026.png")
