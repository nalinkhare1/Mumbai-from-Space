"""
Mumbai from Space, 1988-2026: how much land has Greater Mumbai taken from the sea, and when?

Needs epochs/*.npz from landsat_epochs.py. A 30 m pixel counts as new land when:
  sea in 1988-91  = wet (MNDWI > 0) in at least half the clear scenes, in water connected to open
                    sea/creeks (stricter than the Sentinel-2 rule: at 30 m, rocky shore and flooded
                    mangrove read wet at high tide)
  land in 2026    = dry in every clear 2026 scene (at least 2 clear looks)
New land is split by its 2026 NDVI: bare/built (< 0.2, reclamation) vs vegetated (mostly mangroves
spreading over mudflats). Each pixel is dated to the first epoch from which it stayed land.
Patches under 5 ha are dropped.

python coastline_1988_2026.py
"""
import glob
import json
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager
from pyproj import Transformer
from scipy import ndimage
from shapely import contains_xy
from shapely.geometry import shape

PIX_HA = 0.09  # 30 m x 30 m
MIN_PATCH_PX = 56  # ~5 ha
BARE_NDVI = 0.2
SEA_FREQ = 0.5  # sea in 1988-91 = water in at least half of the clear scenes

files = sorted(glob.glob("epochs/*.npz"))
labels = [os.path.basename(f)[:-4] for f in files]
E = [np.load(f) for f in files]
xs, ys = E[0]["xs"], E[0]["ys"]
assert all(e["mndwi"].shape[1:] == E[0]["mndwi"].shape[1:] for e in E), "grids differ"

wet = [(e["mndwi"].astype("float32") > 0).any(0) for e in E]
looks = [np.isfinite(e["mndwi"].astype("float32")).sum(0) for e in E]
land = [~w & (n >= 2) for w, n in zip(wet, looks)]
unknown = [n < 2 for n in looks]

# --- sea in the first epoch: wet in most clear looks, connected to the image edge ---
# (at 30 m, "wet at least once" catches rocky shore and flooded mangrove at high tide)
m0 = E[0]["mndwi"].astype("float32")
water_freq0 = (m0 > 0).sum(0) / np.maximum(np.isfinite(m0).sum(0), 1)
lab, _ = ndimage.label(water_freq0 >= SEA_FREQ)
edge = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
sea0 = np.isin(lab, edge[edge > 0])

# --- new land: sea then, land now ---
change = ndimage.binary_opening(sea0 & land[-1])
lab, n = ndimage.label(change)
sizes = ndimage.sum(change, lab, range(1, n + 1))

# --- when did it appear? first epoch from which the pixel stayed land (unknown looks pass) ---
stays = np.ones_like(land[-1])
first = np.full(land[-1].shape, -1, dtype=int)
for i in range(len(E) - 1, 0, -1):
    stays = stays & (land[i] | unknown[i])
    first[stays & land[i]] = i

ndvi_now = E[-1]["ndvi"].astype("float32")
bare = ndvi_now < BARE_NDVI

boundary = shape(json.load(open("greater_mumbai.geojson"))["geometry"]).buffer(0.003)
to_ll = Transformer.from_crs("EPSG:32643", "EPSG:4326", always_xy=True)
LON, LAT = to_ll.transform(*np.meshgrid(xs, ys))
inside = contains_xy(boundary, LON, LAT)

sites = []
keep = np.zeros_like(change)
for i in 1 + np.flatnonzero(sizes >= MIN_PATCH_PX):
    m = lab == i
    if inside[m].mean() < 0.5:
        continue
    keep |= m
    dated = first[m][first[m] > 0]
    sites.append(dict(
        ha=round(m.sum() * PIX_HA, 1),
        bare_ha=round((m & bare).sum() * PIX_HA, 1),
        lat=round(float(LAT[m].mean()), 4), lon=round(float(LON[m].mean()), 4),
        appeared=labels[int(np.median(dated))] if dated.size else "?"))
sites.sort(key=lambda s: -s["ha"])

new_bare, new_green = keep & bare, keep & ~bare
print(f"Epochs: {', '.join(labels)}")
print(f"New land in Greater Mumbai since {labels[0]} (patches >= 5 ha): {keep.sum() * PIX_HA:.0f} ha")
print(f"  bare / built (reclamation): {new_bare.sum() * PIX_HA:.0f} ha")
print(f"  vegetated (mostly mangroves on former mudflat): {new_green.sum() * PIX_HA:.0f} ha")
print("\nBare/built new land by the epoch it appeared:")
for i in range(1, len(E)):
    print(f"  {labels[i]:>8}: {(new_bare & (first == i)).sum() * PIX_HA:6.0f} ha")
print("\nSites (largest first):")
for s in sites:
    print(f"  {s['ha']:6.1f} ha ({s['bare_ha']:6.1f} bare)  appeared ~{s['appeared']:>8}  "
          f"({s['lat']}, {s['lon']})")
json.dump(dict(epochs=labels, sites=sites), open("sites_1988_2026.json", "w"), indent=1)

# --- visuals: two zoomed panels per epoch (Coastal Road, Thane creek), then the new land ---
FRAME_EPOCHS = ["1988-91", "1993-95", "1997-99", "2001-03", "2014", "2018", "2022", "2026"]
# (2005-07 and 2009-11 are skipped in the visuals only: Landsat 7's scan-line gaps leave stripes)
PANELS = {"Coastal Road": (18.950, 19.030, 72.790, 72.835),   # lat0, lat1, lon0, lon1
          "Thane creek": (19.000, 19.150, 72.915, 72.995)}
PANEL_H = 900
font = ImageFont.truetype(font_manager.findfont("DejaVu Sans:bold"), 44)
small = ImageFont.truetype(font_manager.findfont("DejaVu Sans:bold"), 24)
to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32643", always_xy=True)


def window(lat0, lat1, lon0, lon1):
    (x0, x1), (y0, y1) = to_utm.transform([lon0, lon1], [lat0, lat1])
    cols = np.flatnonzero((xs >= x0) & (xs <= x1))
    rows = np.flatnonzero((ys >= y0) & (ys <= y1))
    return slice(rows.min(), rows.max() + 1), slice(cols.min(), cols.max() + 1)


wins = {k: window(*v) for k, v in PANELS.items()}


def stretch(rgb):
    gap = np.isnan(rgb).any(-1)
    if gap.any():  # pixels never seen clearly: borrow the nearest clear pixel
        idx = ndimage.distance_transform_edt(gap, return_distances=False, return_indices=True)
        rgb = rgb[idx[0], idx[1]]
    lo, hi = np.nanpercentile(rgb, [2, 98])
    return np.clip((rgb - lo) / (hi - lo), 0, 1)


def overlay(img, sl):
    img = img.copy()
    g, b = new_green[sl], new_bare[sl]
    img[g] = 0.35 * img[g] + 0.65 * np.array([0.25, 1.0, 0.4])
    img[b] = 0.2 * img[b] + 0.8 * np.array([1.0, 0.15, 0.1])
    return img


def panel(img, name):
    im = Image.fromarray((img * 255).astype("uint8"))
    im = im.resize((round(im.width * PANEL_H / im.height), PANEL_H), Image.LANCZOS)
    ImageDraw.Draw(im).text((14, PANEL_H - 40), name, fill="white", font=small,
                            stroke_width=2, stroke_fill="black")
    return im


def compose(imgs, title, note=""):
    w = sum(i.width for i in imgs) + 12 * (len(imgs) - 1)
    out = Image.new("RGB", (w, PANEL_H + 100), "black")
    x = 0
    for im in imgs:
        out.paste(im, (x, 100))
        x += im.width + 12
    d = ImageDraw.Draw(out)
    d.text((14, 6), title, fill="white", font=font)
    d.text((14, 62), note, fill="white", font=small)
    return out


frames = []
for label in FRAME_EPOCHS:
    rgb = E[labels.index(label)]["rgb"].astype("float32")
    frames.append(compose([panel(stretch(rgb[sl]), k) for k, sl in wins.items()], label,
                          "Landsat, dry season (Jan-Apr)"))
rgb_now = E[-1]["rgb"].astype("float32")
frames.append(compose([panel(overlay(stretch(rgb_now[sl]), sl), k) for k, sl in wins.items()],
                      "1988-91 \u2192 2026", "red = built on water   green = mangroves"))

os.makedirs("images", exist_ok=True)
frames[0].save("images/coastline_1988_2026.gif", save_all=True, append_images=frames[1:],
               duration=[1000] * (len(frames) - 1) + [4000], loop=0, optimize=True)
frames[-1].save("images/new_land_zooms.png")

# whole city, outside Greater Mumbai dimmed
city = overlay(stretch(rgb_now), (slice(None), slice(None)))
city[~inside] *= 0.35
Image.fromarray((city * 255).astype("uint8")).save("images/new_land_1988_2026.png")
print("\nSaved images/coastline_1988_2026.gif, images/new_land_zooms.png, images/new_land_1988_2026.png")
