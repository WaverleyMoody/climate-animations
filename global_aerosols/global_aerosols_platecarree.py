"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: global_aerosols_platecarree.py
Description: Generates the global aerosols animation from NASA MERRA-2 hourly
    aerosol extinction AOT (dust, carbonaceous, sulfate, sea salt) and total
    precipitable water vapor (2018-08-01 to 2020-01-31), rendered in the
    PlateCarree projection.
Note: For the data download, see download_aerosols_merra2.py in the
    global_aerosols scripts folder.

Run:
    conda activate climate
    # 1) Set PREVIEW_ONLY = True, check the PNG against the reference frame
    # 2) Set PREVIEW_ONLY = False for the full render
    caffeinate -i python global_aerosols_platecarree.py
"""

import re
import shutil
import sys
from pathlib import Path

import cartopy.io.shapereader as shpreader
import dask
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
import xarray as xr
from matplotlib.animation import FFMpegWriter
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, LogNorm, Normalize
from PIL import Image
from scipy.ndimage import map_coordinates
from shapely.ops import unary_union

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
DATA_DIR = Path("/Volumes/CLIMATEDATA/global_aerosols/merra2_hourly")
BLUE_MARBLE_DIR = Path("/Volumes/CLIMATEDATA/blue_marble/2004")

LOCAL_DIR = Path.home() / "CLIMATE ANIMATIONS" / "global_aerosols" / "_render"   # local APFS
FINAL_DIR = Path("/Volumes/CLIMATEDATA/global_aerosols/animations")              # exFAT
OUT_NAME = "global_aerosols_2018_2020_platecarree.mp4"

PREVIEW_ONLY = False                       # True: save one PNG and stop
PREVIEW_TIME = "2019-08-02T14:30"         # matches the reference frame

HOUR_STEP = 3                             # every 3rd hour -> 4,392 frames
FPS = 30                                  # 4,392 frames -> ~2:26 (original 2:13)
BITRATE = 14000

# Composited map raster (2:1 equirectangular) and full video frame, in pixels
WIDTH, HEIGHT = 1920, 960
FIG_W, FIG_H = 1920, 1200
DPI = 100

TITLE = "Global Aerosols, August 2018 – January 2020"
SUBTITLE = "NASA MERRA-2 hourly aerosol optical thickness (550 nm) by species, with total precipitable water vapor"
CREDIT = "Data: NASA GMAO MERRA-2 (M2T1NXAER, M2T1NXSLV)  ·  Base map: NASA Blue Marble Next Generation (2004)"

LON_TICKS = np.arange(-180, 181, 90)
LAT_TICKS = np.arange(-90, 91, 45)
GRIDLINES = True

OCEAN_COLOR = "#0b2550"                   # deep-ocean base under the sea salt layer
LAND_KEY_GRAY = 0.55                      # backdrop used to draw the land-layer color keys
SEA_ICE_LUM = 0.35                        # Blue Marble ocean pixels brighter than this = sea ice

# Layers, composited bottom -> top.
#   vars:   summed MERRA-2 variables
#   scale:  "log" or "linear" between vmin and vmax
#   colors: color ramp from low to high
#   alpha:  maximum opacity; gamma shapes how fast opacity rises
#   ocean_only: draw only over ocean pixels
#   label, ticks, key_bg: color-key text, tick values, and backdrop ("ocean" or "land")
LAYERS = [
    dict(name="sea salt", vars=["SSEXTTAU"], scale="log", vmin=0.02, vmax=0.40,
         colors=["#0d3a6e", "#2a7fb8", "#8fd4d9", "#f0eebf"],
         alpha=0.95, gamma=0.7, ocean_only=True,
         label="Sea salt AOT", ticks=[0.02, 0.05, 0.1, 0.2, 0.4], key_bg="ocean"),
    dict(name="water vapor", vars=["TQV"], scale="linear", vmin=20.0, vmax=65.0,
         colors=["#c8ccd0", "#e6e8ea"],
         alpha=0.25, gamma=1.5, ocean_only=False,
         label="Water vapor (kg m$^{-2}$)", ticks=[20, 35, 50, 65], key_bg="ocean"),
    dict(name="sulfate", vars=["SUEXTTAU"], scale="log", vmin=0.05, vmax=0.60,
         colors=["#d9d9d9", "#ffffff"],
         alpha=0.75, gamma=1.0, ocean_only=False,
         label="Sulfate AOT", ticks=[0.05, 0.1, 0.2, 0.6], key_bg="land"),
    dict(name="carbon / smoke", vars=["BCEXTTAU", "OCEXTTAU"], scale="log", vmin=0.04, vmax=3.0,
         colors=["#3f9e4a", "#7cc25a", "#c8331b", "#ff5a1a", "#ffd27a", "#fff6e0"],
         alpha=0.95, gamma=0.6, ocean_only=False,
         label="Smoke / carbon AOT", ticks=[0.04, 0.1, 0.3, 1, 3], key_bg="land"),
    dict(name="dust", vars=["DUEXTTAU"], scale="log", vmin=0.05, vmax=2.0,
         colors=["#c56a1c", "#ef9a2e", "#f9c35a", "#fde9b0"],
         alpha=0.92, gamma=0.7, ocean_only=False,
         label="Dust AOT", ticks=[0.05, 0.1, 0.2, 0.5, 1, 2], key_bg="land"),
]
KEY_ORDER = ["dust", "carbon / smoke", "sulfate", "sea salt", "water vapor"]   # left -> right
ALL_VARS = sorted({v for layer in LAYERS for v in layer["vars"]})

dask.config.set(scheduler="synchronous")
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()


# ---------------------------------------------------------------------------
# Static fields (built once)
# ---------------------------------------------------------------------------
def _target_grid():
    lat_c = 90.0 - (np.arange(HEIGHT) + 0.5) * 180.0 / HEIGHT
    lon_c = -180.0 + (np.arange(WIDTH) + 0.5) * 360.0 / WIDTH
    return lat_c, lon_c


def _interp_coords(src_lat, src_lon):
    """Fractional MERRA-2 indices for every output pixel (lon wraps via a padded column)."""
    lat_c, lon_c = _target_grid()
    dlat = src_lat[1] - src_lat[0]
    dlon = src_lon[1] - src_lon[0]
    rows = (lat_c - src_lat[0]) / dlat
    cols = (lon_c - src_lon[0]) / dlon
    rr, cc = np.meshgrid(rows, cols, indexing="ij")
    return np.stack([rr, cc]).astype(np.float32)


def _land_mask():
    cache = LOCAL_DIR / f"land_mask_{WIDTH}x{HEIGHT}.npy"
    if cache.exists():
        return np.load(cache)
    print("Building land mask (one-time)...")
    path = shpreader.natural_earth(resolution="50m", category="physical", name="land")
    land = unary_union([rec.geometry for rec in shpreader.Reader(path).records()])
    shapely.prepare(land)
    lat_c, lon_c = _target_grid()
    lon2, lat2 = np.meshgrid(lon_c, lat_c)
    mask = shapely.contains_xy(land, lon2, lat2)
    np.save(cache, mask)
    return mask


def _base_maps(land):
    """Per-month base image: desaturated Blue Marble on land/sea ice, navy ocean."""
    ocean_rgb = np.array(matplotlib.colors.to_rgb(OCEAN_COLOR), dtype=np.float32)
    bases = {}
    for m in range(1, 13):
        f = BLUE_MARBLE_DIR / f"world.2004{m:02d}.3x5400x2700.jpg"
        img = Image.open(f).convert("RGB").resize((WIDTH, HEIGHT), Image.LANCZOS)
        rgb = np.asarray(img, dtype=np.float32) / 255.0
        lum = 0.299 * rgb[..., 0] + 0.587 * rgb[..., 1] + 0.114 * rgb[..., 2]
        gray = (0.35 + 0.65 * lum)[..., None] * np.ones(3, dtype=np.float32)
        is_ocean_water = (~land) & (lum < SEA_ICE_LUM)
        base = np.where(is_ocean_water[..., None], ocean_rgb, gray)
        bases[m] = base.astype(np.float32)
    return bases


def _build_ramps():
    return [LinearSegmentedColormap.from_list(L["name"], L["colors"]) for L in LAYERS]


# ---------------------------------------------------------------------------
# Per-frame compositing
# ---------------------------------------------------------------------------
def _upsample(field, coords):
    padded = np.concatenate([field, field[:, :1]], axis=1)   # wrap 180E -> 180W
    return map_coordinates(padded, coords, order=1, mode="nearest")


def _intensity(x, layer):
    x = np.maximum(x, 1e-6)
    if layer["scale"] == "log":
        lo, hi = np.log10(layer["vmin"]), np.log10(layer["vmax"])
        t = (np.log10(x) - lo) / (hi - lo)
    else:
        t = (x - layer["vmin"]) / (layer["vmax"] - layer["vmin"])
    return np.clip(t, 0.0, 1.0)


def _composite(fields, base, ocean, coords, ramps):
    out = base.copy()
    for layer, ramp in zip(LAYERS, ramps):
        raw = sum(fields[v] for v in layer["vars"])
        t = _intensity(_upsample(raw, coords), layer)
        color = ramp(t)[..., :3].astype(np.float32)
        a = (layer["alpha"] * t ** layer["gamma"]).astype(np.float32)
        if layer["ocean_only"]:
            a = a * ocean
        a = a[..., None]
        out = out * (1.0 - a) + color * a
    return (np.clip(out, 0.0, 1.0) * 255).astype(np.uint8)


# ---------------------------------------------------------------------------
# Figure: map, axes labels, title, color keys
# ---------------------------------------------------------------------------
def _key_cmap(layer, ramp):
    """Colormap showing the layer as it actually appears: ramp blended over its backdrop."""
    t = np.linspace(0.0, 1.0, 256)
    color = ramp(t)[:, :3]
    a = (layer["alpha"] * t ** layer["gamma"])[:, None]
    if layer["key_bg"] == "ocean":
        bg = np.array(matplotlib.colors.to_rgb(OCEAN_COLOR))
    else:
        bg = np.full(3, LAND_KEY_GRAY)
    return ListedColormap(bg * (1.0 - a) + color * a, name=f"{layer['name']}_key")


def _lon_label(v):
    v = int(v)
    if v in (-180, 180):
        return "180°"
    if v == 0:
        return "0°"
    return f"{abs(v)}°{'E' if v > 0 else 'W'}"


def _lat_label(v):
    v = int(v)
    if v == 0:
        return "0°"
    return f"{abs(v)}°{'N' if v > 0 else 'S'}"


def _build_figure(ramps):
    fig = plt.figure(figsize=(FIG_W / DPI, FIG_H / DPI), dpi=DPI, facecolor="white")

    # Map axes sized to an exact 2:1 pixel box so PlateCarree is undistorted
    map_w = 0.86
    map_h = map_w * (FIG_W / FIG_H) / 2.0
    ax = fig.add_axes([0.07, 0.215, map_w, map_h])
    im = ax.imshow(np.zeros((HEIGHT, WIDTH, 3), dtype=np.uint8), origin="upper",
                   extent=[-180, 180, -90, 90], interpolation="antialiased", aspect="auto")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(LON_TICKS)
    ax.set_xticklabels([_lon_label(v) for v in LON_TICKS], fontsize=16)
    ax.set_yticks(LAT_TICKS)
    ax.set_yticklabels([_lat_label(v) for v in LAT_TICKS], fontsize=16)
    ax.tick_params(direction="out", length=4, width=0.8, top=True, right=True,
                   labeltop=False, labelright=False)
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    if GRIDLINES:
        ax.set_axisbelow(False)
        ax.grid(color="white", alpha=0.22, linewidth=0.6)

    # Title, subtitle, timestamp
    top = 0.215 + map_h
    fig.text(0.07, top + 0.050, TITLE, ha="left", va="bottom", fontsize=32, fontweight="bold")
    fig.text(0.07, top + 0.018, SUBTITLE, ha="left", va="bottom", fontsize=18, color="#333333")
    stamp = fig.text(0.93, top + 0.050, "", ha="right", va="bottom", fontsize=28,
                     family="monospace")

    # Color keys: one horizontal bar per layer
    by_name = {L["name"]: (L, r) for L, r in zip(LAYERS, ramps)}
    n = len(KEY_ORDER)
    gap = 0.05
    bar_w = (map_w - gap * (n - 1)) / n
    for i, name in enumerate(KEY_ORDER):
        layer, ramp = by_name[name]
        norm = (LogNorm(layer["vmin"], layer["vmax"]) if layer["scale"] == "log"
                else Normalize(layer["vmin"], layer["vmax"]))
        cax = fig.add_axes([0.07 + i * (bar_w + gap), 0.105, bar_w, 0.022])
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=_key_cmap(layer, ramp)),
                          cax=cax, orientation="horizontal")
        cb.set_ticks(layer["ticks"])
        cb.set_ticklabels([f"{v:g}" for v in layer["ticks"]])
        cb.ax.minorticks_off()
        cb.ax.tick_params(labelsize=15)
        cb.outline.set_linewidth(0.6)
        cax.set_title(layer["label"], fontsize=16, pad=7)

    fig.text(0.93, 0.015, CREDIT, ha="right", va="bottom", fontsize=13, color="#555555")
    return fig, im, stamp


# ---------------------------------------------------------------------------
# Files, verification
# ---------------------------------------------------------------------------
def _day_files():
    files = []
    for f in sorted(DATA_DIR.iterdir()):
        if f.name.startswith("._"):
            continue
        if re.fullmatch(r"merra2_aerosols_\d{8}\.nc", f.name):
            files.append(f)
    return files


def _stamp(t):
    return pd.Timestamp(t).strftime("%Y-%m-%d %H:%M UTC")


def _verify_mp4(path, expected_frames):
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"{path.name} is missing or empty")
    n_frames, secs = imageio_ffmpeg.count_frames_and_secs(str(path))
    if n_frames != expected_frames:
        raise ValueError(f"{path.name}: {n_frames} frames, expected {expected_frames}")
    print(f"Verified {path.name}: {n_frames} frames, {secs:.1f} s")


def _load_day(path):
    with xr.open_dataset(path) as ds:
        ds = ds[ALL_VARS].load()
    return ds


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if not DATA_DIR.exists():
        sys.exit(f"{DATA_DIR} not found. Is /Volumes/CLIMATEDATA mounted?")
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)

    files = _day_files()
    print(f"Found {len(files)} daily files")

    first = _load_day(files[0])
    src_lat, src_lon = first.lat.values, first.lon.values
    assert np.isclose(src_lon[0], -180.0) and np.isclose(src_lat[0], -90.0), "Unexpected MERRA-2 grid"
    coords = _interp_coords(src_lat, src_lon)

    land = _land_mask()
    ocean = (~land).astype(np.float32)
    bases = _base_maps(land)
    ramps = _build_ramps()
    fig, im, stamp = _build_figure(ramps)

    if PREVIEW_ONLY:
        t = pd.Timestamp(PREVIEW_TIME)
        day = DATA_DIR / f"merra2_aerosols_{t:%Y%m%d}.nc"
        ds = _load_day(day).sel(time=t, method="nearest")
        fields = {v: ds[v].values for v in ALL_VARS}
        im.set_data(_composite(fields, bases[t.month], ocean, coords, ramps))
        stamp.set_text(_stamp(ds.time.values))
        out = LOCAL_DIR / "preview_frame.png"
        fig.savefig(out, dpi=DPI)
        print(f"Preview saved: {out}")
        return

    local_mp4 = LOCAL_DIR / OUT_NAME
    writer = FFMpegWriter(fps=FPS, codec="libx264", bitrate=BITRATE,
                          extra_args=["-pix_fmt", "yuv420p"])
    n_written = 0
    with writer.saving(fig, str(local_mp4), dpi=DPI):
        for d, path in enumerate(files, 1):
            ds = _load_day(path)
            for i in range(0, ds.sizes["time"], HOUR_STEP):
                frame = ds.isel(time=i)
                fields = {v: frame[v].values for v in ALL_VARS}
                month = pd.Timestamp(frame.time.values).month
                im.set_data(_composite(fields, bases[month], ocean, coords, ramps))
                stamp.set_text(_stamp(frame.time.values))
                writer.grab_frame()
                n_written += 1
            ds.close()
            if d % 10 == 0 or d == len(files):
                print(f"  day {d}/{len(files)}  {path.name}  ({n_written} frames)")
    plt.close(fig)

    _verify_mp4(local_mp4, n_written)

    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    final_mp4 = FINAL_DIR / OUT_NAME
    if final_mp4.exists():
        final_mp4.unlink()
    shutil.move(str(local_mp4), str(final_mp4))
    print(f"Saved {final_mp4}")


if __name__ == "__main__":
    main()