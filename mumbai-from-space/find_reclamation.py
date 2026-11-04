"""
Mumbai from Space: find every patch of land reclaimed from the sea across Greater Mumbai,
dry season 2018 vs dry season 2026, from free Sentinel-2 imagery.

A 10 m pixel counts as reclaimed when:
  sea in 2018   = wet (MNDWI > 0) in at least one clear scene, in a water body connected to
                  open sea/creeks (lakes like Powai and Vihar don't reach the edge, so are excluded)
  land in 2026  = dry in every clear scene (at least 2 clear looks)
  bare in 2026  = NDVI < 0.2 (on Mumbai's creeks, "new land" that is green is mangroves
                  advancing over mudflats, not reclamation)
Patches under 5 ha are dropped (tide and beach-edge noise); 2-5 ha patches are listed separately.

Boundary: Mumbai City + Mumbai Suburban districts (geoBoundaries gbOpen IND ADM2, CC BY 4.0).

pip install -r requirements.txt
python find_reclamation.py
"""
import json
import os
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from shapely.geometry import box, shape
from shapely.ops import unary_union
from shapely import contains_xy
from pyproj import Transformer
import pystac_client
import planetary_computer
from odc.stac import load

BOUNDARY = shape(json.load(open("greater_mumbai.geojson"))["geometry"])
lon0, lat0, lon1, lat1 = BOUNDARY.bounds
BBOX = (lon0 - 0.01, lat0 - 0.01, lon1 + 0.01, lat1 + 0.01)
PERIODS = {"2018": "2018-01-01/2018-04-30", "2026": "2026-01-01/2026-04-30"}
N_DAYS = 6
CRS = "EPSG:32643"
BANDS = dict(green="B03", swir="B11", nir="B08", red="B04", blue="B02", scl="SCL")
BAD_SCL = [0, 1, 3, 8, 9, 10]  # nodata, saturated, cloud shadow, clouds, cirrus
MIN_PATCH_PX = 500  # 5 ha at 10 m
SMALL_PATCH_PX = 200  # 2 ha: listed separately for checking by eye
BARE_NDVI = 0.2  # Coastal Road fill reads ~0.04; mangroves advancing on the creeks read ~0.65
CACHE = "city_cache.npz"  # delete to re-download

cat = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",
                                modifier=planetary_computer.sign_inplace)


def best_days(period):
    """Days whose tiles together cover the whole city, least cloudy first."""
    items = cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=period,
                       query={"eo:cloud_cover": {"lt": 20}}).item_collection()
    by_day = defaultdict(list)
    for it in items:
        by_day[it.datetime.strftime("%Y-%m-%d")].append(it)
    full = [(np.mean([i.properties["eo:cloud_cover"] for i in its]), d, its)
            for d, its in by_day.items()
            if unary_union([shape(i.geometry) for i in its]).contains(box(*BBOX))]
    return [(d, its) for _, d, its in sorted(full)[:N_DAYS]]


def period_stack(period):
    mndwi, ndvi, rgbs = [], [], []
    for day, its in best_days(period):
        ds = load(its, bands=list(BANDS.values()), bbox=BBOX, crs=CRS, resolution=10,
                  groupby="solar_day", resampling={"B11": "bilinear"}).isel(time=0)
        f = {k: ds[v].values.astype("float32") for k, v in BANDS.items() if k != "scl"}
        ok = ~np.isin(ds["SCL"].values, BAD_SCL) & (f["green"] > 0) & (f["swir"] > 0)
        # ESA processing baseline 04.00+ (Jan 2022 on) stores reflectance with a +1000 offset
        if its[0].properties["s2:processing_baseline"] >= "04.00":
            f = {k: np.maximum(v - 1000, 1) for k, v in f.items()}
        mw = (f["green"] - f["swir"]) / (f["green"] + f["swir"])
        mndwi.append(np.where(ok, mw, np.nan).astype("float16"))
        nd = (f["nir"] - f["red"]) / np.maximum(f["nir"] + f["red"], 1)
        ndvi.append(np.where(ok, nd, np.nan).astype("float16"))
        if len(rgbs) < 3:
            rgbs.append(np.stack([f["red"], f["green"], f["blue"]], -1).astype("uint16"))
        print(f"  {day}: {ok.mean():.0%} clear")
        xs, ys = ds.x.values, ds.y.values
    rgb = np.median(np.stack(rgbs), axis=0)
    return np.stack(mndwi), np.nanmedian(np.stack(ndvi).astype("float32"), 0), rgb, xs, ys


try:
    c = np.load(CACHE)
    m18, m26, ndvi26, rgb18, rgb26, xs, ys = (c[k] for k in
                                              ["m18", "m26", "ndvi26", "rgb18", "rgb26", "xs", "ys"])
    print("Using cached scenes")
except (FileNotFoundError, KeyError):
    print("2018:"); m18, _, rgb18, xs, ys = period_stack(PERIODS["2018"])
    print("2026:"); m26, ndvi26, rgb26, _, _ = period_stack(PERIODS["2026"])
    assert m18.shape[1:] == m26.shape[1:], "2018 and 2026 grids differ"
    np.savez_compressed(CACHE, m18=m18, m26=m26, ndvi26=ndvi26, rgb18=rgb18, rgb26=rgb26,
                        xs=xs, ys=ys)

m18, m26 = m18.astype("float32"), m26.astype("float32")

# --- sea in 2018: wet at least once, in a water body that reaches the image edge ---
wet18 = (m18 > 0).any(axis=0)
lab, _ = ndimage.label(wet18)
edge_ids = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
sea18 = np.isin(lab, edge_ids[edge_ids > 0])

# --- bare land in 2026: dry in every clear look, and not vegetated ---
# (green "new land" on the creeks is mangroves advancing over mudflats, not reclamation)
looks26 = (~np.isnan(m26)).sum(0)
land26 = ~(m26 > 0).any(0) & (looks26 >= 2) & (ndvi26 < BARE_NDVI)

change = ndimage.binary_opening(sea18 & land26)
lab, n = ndimage.label(change)
sizes = ndimage.sum(change, lab, range(1, n + 1))

# --- inside Greater Mumbai? (pixel centres) ---
to_ll = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
X, Y = np.meshgrid(xs, ys)
LON, LAT = to_ll.transform(X, Y)
inside = contains_xy(BOUNDARY.buffer(0.003), LON, LAT)  # ~300 m: new land sits past the old coastline

sites = []
for i in 1 + np.flatnonzero(sizes >= SMALL_PATCH_PX):
    m = lab == i
    if inside[m].mean() < 0.5:
        continue  # Navi Mumbai / Thane side of the creek
    sites.append(dict(id=int(i), ha=round(m.sum() * 0.01, 1),
                      lat=round(float(LAT[m].mean()), 4), lon=round(float(LON[m].mean()), 4)))
sites.sort(key=lambda s: -s["ha"])
big = [s for s in sites if s["ha"] * 100 >= MIN_PATCH_PX]
small = [s for s in sites if s["ha"] * 100 < MIN_PATCH_PX]

print(f"\nLand reclaimed from the sea, Greater Mumbai, 2018 -> 2026 (patches >= 5 ha): "
      f"{sum(s['ha'] for s in big):.1f} ha in {len(big)} patches")
for s in big:
    print(f"  {s['ha']:6.1f} ha  ({s['lat']}, {s['lon']})")
print(f"Smaller patches (2-5 ha), check by eye before believing: {len(small)}")
for s in small:
    print(f"  {s['ha']:6.1f} ha  ({s['lat']}, {s['lon']})")
json.dump(dict(reclaimed=big, smaller_unverified=small), open("sites.json", "w"), indent=1)


# --- map: whole city + zoom on the biggest change ---
def stretch(rgb):
    lo, hi = np.nanpercentile(rgb, [2, 98])
    return np.clip(np.nan_to_num((rgb - lo) / (hi - lo)), 0, 1)


mask = np.isin(lab, [s["id"] for s in big])
img = stretch(rgb26)
img[mask] = 0.25 * img[mask] + 0.75 * np.array([1.0, 0.15, 0.1])
rows, cols = np.nonzero(mask)
r0, r1 = rows.min() - 150, rows.max() + 150
c0, c1 = max(cols.min() - 250, 0), cols.max() + 250

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 13), gridspec_kw=dict(width_ratios=[1, 0.9]))
ax1.imshow(img)
ax1.add_patch(plt.Rectangle((c0, r0), c1 - c0, r1 - r0, fill=False, ec="white", lw=1.5))
ax1.set_title("Greater Mumbai, 2026", fontsize=14)
ax2.imshow(img[r0:r1, c0:c1])
ax2.set_title(f"Reclaimed from the sea since 2018: {mask.sum() * 0.01:.0f} ha", fontsize=14)
for ax in (ax1, ax2):
    ax.axis("off")
fig.suptitle("The only large land reclamation in Greater Mumbai since 2018: the Coastal Road",
             fontsize=16)
fig.text(0.5, 0.015, "Sentinel-2, 10 m, dry season 2018 vs 2026, patches >= 5 ha. "
         "Contains modified Copernicus Sentinel data", ha="center", fontsize=9)
plt.tight_layout(rect=(0, 0.03, 1, 0.96))
os.makedirs("images", exist_ok=True)
plt.savefig("images/greater_mumbai_new_land.png", dpi=130)
print("Saved images/greater_mumbai_new_land.png and sites.json")


# --- before/after GIF of the zoom, for the README ---
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager

font = ImageFont.truetype(font_manager.findfont("DejaVu Sans:bold"), 34)


def frame(arr, label):
    im = Image.fromarray((arr * 255).astype("uint8"))
    h = 760
    im = im.resize((round(im.width * h / im.height), h), Image.LANCZOS)
    d = ImageDraw.Draw(im)
    d.text((18, 14), label, fill="white", font=font, stroke_width=3, stroke_fill="black")
    return im


crop = (slice(r0, r1), slice(c0, c1))
a18, a26 = stretch(rgb18[crop]), stretch(rgb26[crop])
over = a26.copy()
over[mask[crop]] = 0.25 * over[mask[crop]] + 0.75 * np.array([1.0, 0.15, 0.1])
frames = [frame(a18, "2018")] * 6
frames += [frame((1 - t) * a18 + t * a26, "") for t in np.linspace(0.15, 0.85, 6)]
frames += [frame(a26, "2026")] * 6
frames += [frame(over, f"+{mask.sum() * 0.01:.0f} ha")] * 8
frames[0].save("images/coastal_road_before_after.gif", save_all=True, append_images=frames[1:],
               duration=250, loop=0, optimize=True)
print("Saved images/coastal_road_before_after.gif")
