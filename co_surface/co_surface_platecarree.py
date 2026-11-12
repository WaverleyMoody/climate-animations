"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: co_surface_platecarree.py
Description: Generates the surface carbon monoxide animation from NASA MERRA-2
    hourly CO surface concentration COSC (2006-12-01 to 2006-12-31), rendered
    in the Plate Carree projection over the December Blue Marble.
Note: For the data download, see download_merra2_cosc_2006_12.py in the
    co_surface scripts folder.

Look:
    Low CO is transparent, so the Blue Marble shows through (navy ocean, white
    Antarctica). CO fades in as teal/cyan from ~70 ppb and is fully opaque from
    ~160 ppb, then runs yellow -> red -> black -> grey -> white up to 1600 ppb,
    matching the original's colorbar (40-1600 ppb, non-linear ticks).

Run:
    conda activate climate
    # 1) Set PREVIEW_ONLY = True, check the PNG
    # 2) Set PREVIEW_ONLY = False for the full render
    caffeinate -i python co_surface_platecarree.py
"""

import shutil
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import dask
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FFMpegWriter
from matplotlib.cm import ScalarMappable
from matplotlib.colors import FuncNorm, LinearSegmentedColormap, to_rgb
from PIL import Image

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
VARIABLE = "CO_surface"                   # ppbv, from download_merra2_cosc_2006_12.py
DRIVE = Path("/Volumes/CLIMATEDATA")
DATA_PATH = DRIVE / "co" / "merra2_2006_12" / "MERRA2_COSC_hourly_2006-12.nc"
BM_DIR = DRIVE / "blue_marble" / "2004"
BM_MONTH = 12                             # December snow cover, like the original

LOCAL_DIR = Path.home() / "CLIMATE ANIMATIONS" / "co_surface" / "_render"   # APFS
FINAL_DIR = DRIVE / "co" / "animations"                                     # exFAT
OUT_NAME = "co_surface_dec2006_platecarree.mp4"

PREVIEW_ONLY = False                       # True: save one PNG and stop
PREVIEW_TIME = "2006-12-07T03:30"         # matches the reference frame (0300)

CENTRAL_LON = 165.0                       # Pacific-centered, as in the original
SHOW_COASTLINES = True
BM_BRIGHTNESS = 0.9

# Colorbar ticks from the original; equally spaced on the bar (non-linear scale)
TICKS = [40, 100, 160, 200, 400, 600, 800, 1200, 1600]

# (value ppb, color, alpha). Alpha ramps in so clean air shows the Blue Marble.
ANCHORS = [
    (40,   "#0b1622", 0.00),
    (70,   "#12405a", 0.45),
    (100,  "#1b6f8c", 0.75),
    (140,  "#3cb4cf", 0.90),
    (160,  "#7fd6d0", 1.00),
    (175,  "#f2e03b", 1.00),
    (190,  "#f5a623", 1.00),
    (200,  "#ec5a1c", 1.00),
    (300,  "#c0230f", 1.00),
    (400,  "#1c0705", 1.00),
    (600,  "#3d3d3d", 1.00),
    (800,  "#7d7d7d", 1.00),
    (1200, "#c9c9c9", 1.00),
    (1600, "#ffffff", 1.00),
]

CB_TITLE = "Carbon Monoxide"
CB_SUBTITLE = "(surface concentration in ppb)"

# Frame layout: 1920 x 1200 px = black top bar, 2:1 map, black bottom bar
W, H = 1920, 1200
TOP_PX, MAP_PX = 50, 960
DPI = 100
FPS = 30                                  # 744 frames -> 24.8 s (original 0:25)
BITRATE = 14000

dask.config.set(scheduler="synchronous")
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_TICK_POS = np.linspace(0, 1, len(TICKS))


def _fwd(v):
    return np.interp(v, TICKS, _TICK_POS, left=-0.01, right=1.01)


def _inv(p):
    return np.interp(p, _TICK_POS, TICKS)


def _build_cmaps():
    """Map colormap (with alpha ramp) and an opaque twin for the colorbar."""
    pos = _fwd(np.array([a[0] for a in ANCHORS]))
    rgba = [(*to_rgb(c), a) for _, c, a in ANCHORS]
    opaque = [(*to_rgb(c), 1.0) for _, c, _ in ANCHORS]
    cmap = LinearSegmentedColormap.from_list("co", list(zip(pos, rgba)), N=1024)
    cmap.set_under((0, 0, 0, 0))
    cmap.set_over((1, 1, 1, 1))
    cmap.set_bad((0, 0, 0, 0))
    cmap_bar = LinearSegmentedColormap.from_list("co_bar", list(zip(pos, opaque)), N=1024)
    norm = FuncNorm((_fwd, _inv), vmin=TICKS[0], vmax=TICKS[-1])
    return cmap, cmap_bar, norm


def _load_blue_marble():
    paths = sorted(p for p in BM_DIR.glob(f"world.2004{BM_MONTH:02d}*.jpg")
                   if not p.name.startswith("._"))
    if not paths:
        sys.exit(f"No Blue Marble image for month {BM_MONTH:02d} in {BM_DIR}")
    img = Image.open(paths[0]).convert("RGB")
    # Recenter on CENTRAL_LON, then resize once to the on-screen map size
    arr = np.asarray(img)
    shift = int(round(((CENTRAL_LON - 180.0) + 180.0) / 360.0 * arr.shape[1]))
    arr = np.roll(arr, -shift, axis=1)
    img = Image.fromarray(arr).resize((W, MAP_PX), Image.LANCZOS)
    arr = np.asarray(img).astype(np.float32) / 255.0 * BM_BRIGHTNESS
    return np.clip(arr, 0, 1)


def _stamp(t):
    return pd.Timestamp(t).strftime("%Y-%m-%d %H%M")


def _build_figure(lon_shift_extent):
    cmap, cmap_bar, norm = _build_cmaps()
    proj = ccrs.PlateCarree(central_longitude=CENTRAL_LON)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor="black")
    map_bottom = 1 - (TOP_PX + MAP_PX) / H
    ax = fig.add_axes([0, map_bottom, 1, MAP_PX / H], projection=proj)
    ax.set_extent([-180, 180, -90, 90], crs=proj)
    ax.spines["geo"].set_visible(False)
    ax.set_facecolor("black")

    ax.imshow(_load_blue_marble(), origin="upper", extent=[-180, 180, -90, 90],
              transform=proj, interpolation="antialiased", zorder=0)

    im = ax.imshow(np.full((2, 2), np.nan), origin="lower", extent=lon_shift_extent,
                   transform=proj, cmap=cmap, norm=norm, interpolation="bilinear",
                   zorder=1)

    if SHOW_COASTLINES:
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.45,
                       edgecolor="black", zorder=2)

    # Bottom bar: timestamp (left), colorbar block (right)
    stamp = fig.text(0.115, 0.075, "", color="white", fontsize=24,
                     ha="left", va="center")

    fig.text(0.835, 0.135, CB_TITLE, color="white", fontsize=18,
             fontweight="bold", ha="center", va="center")
    fig.text(0.835, 0.108, CB_SUBTITLE, color="#9a9a9a", fontsize=12,
             ha="center", va="center")
    cax = fig.add_axes([0.71, 0.062, 0.25, 0.026])
    sm = ScalarMappable(norm=norm, cmap=cmap_bar)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=TICKS)
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, pad=6)
    cb.ax.set_xticklabels([str(t) for t in TICKS])
    for i, lab in enumerate(cb.ax.get_xticklabels()):
        end = i in (0, len(TICKS) - 1)
        lab.set_color("white" if end else "#8a8a8a")
        lab.set_fontsize(15 if end else 12)

    return fig, im, stamp


def _verify_mp4(path, expected_frames):
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"{path.name} is missing or empty")
    n_frames, secs = imageio_ffmpeg.count_frames_and_secs(str(path))
    if n_frames != expected_frames:
        raise ValueError(f"{path.name}: {n_frames} frames, expected {expected_frames}")
    print(f"Verified {path.name}: {n_frames} frames, {secs:.1f} s")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if not DATA_PATH.exists():
        sys.exit(f"{DATA_PATH} not found. Is /Volumes/CLIMATEDATA mounted?")
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)

    ds = xr.open_dataset(DATA_PATH, chunks={"time": 1})
    da = ds[VARIABLE].transpose("time", "lat", "lon")
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.isel(lat=slice(None, None, -1))
    lon, lat = da.lon.values, da.lat.values
    if not np.all(np.diff(lon) > 0):
        sys.exit("Longitude axis is not ascending; check the input file")

    # Roll so column 0 is the left edge of the Pacific-centered map
    dlon, dlat = lon[1] - lon[0], lat[1] - lat[0]
    left_edge = CENTRAL_LON - 180.0                       # geographic lon
    start = int(np.argmin(np.abs(((lon - left_edge + 180) % 360) - 180)))
    x0 = ((lon[start] - CENTRAL_LON + 180) % 360) - 180   # projected x of col 0
    extent = [x0 - dlon / 2, x0 - dlon / 2 + lon.size * dlon,
              lat[0] - dlat / 2, lat[-1] + dlat / 2]

    def field(i):
        return np.roll(da.isel(time=i).values, -start, axis=1)

    times = da.time.values
    n = times.size
    print(f"Loaded {VARIABLE}: {n} hourly frames, {lat.size} x {lon.size}")

    fig, im, stamp = _build_figure(extent)

    if PREVIEW_ONLY:
        i = int(np.argmin(np.abs(times - np.datetime64(PREVIEW_TIME))))
        f = field(i)
        p = np.percentile(f, [1, 50, 90, 99, 99.9])
        print(f"{_stamp(times[i])}: p1={p[0]:.0f}  median={p[1]:.0f}  p90={p[2]:.0f}  "
              f"p99={p[3]:.0f}  p99.9={p[4]:.0f}  max={f.max():.0f} ppb")
        im.set_data(f)
        stamp.set_text(_stamp(times[i]))
        out = LOCAL_DIR / "preview_frame.png"
        fig.savefig(out, dpi=DPI, facecolor="black")
        print(f"Preview saved: {out}")
        return

    local_mp4 = LOCAL_DIR / OUT_NAME
    writer = FFMpegWriter(fps=FPS, codec="libx264", bitrate=BITRATE,
                          extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, str(local_mp4), dpi=DPI):
        for i in range(n):
            im.set_data(field(i))
            stamp.set_text(_stamp(times[i]))
            writer.grab_frame(facecolor="black")
            if i % 24 == 0:
                print(f"  frame {i + 1}/{n}  {_stamp(times[i])}")
    plt.close(fig)
    ds.close()

    _verify_mp4(local_mp4, n)

    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    final_mp4 = FINAL_DIR / OUT_NAME
    if final_mp4.exists():
        final_mp4.unlink()
    shutil.move(str(local_mp4), str(final_mp4))
    print(f"Saved {final_mp4}")


if __name__ == "__main__":
    main()
