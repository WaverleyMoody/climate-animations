"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: animate_fires_platecarree.py

Description: Generates the fire carbon emissions time-series animation from
             GFED5.1 weekly totals (March 2003 - January 2022), rendered in the
             Plate Carree projection in the style of the NASA SVS "Fires - a
             global perspective" animation.

Note: For the Robinson, Foucaut, and Nicolosi projections, see the other scripts
      in the fires scripts folder.

Run inside the climate env:  conda activate climate
  1. RENDER_MODE = "test"  -> writes a few PNG frames for review (fast)
  2. RENDER_MODE = "full"  -> caffeinate -i python animate_fires_platecarree.py
If you change any basemap colors, delete BG_CACHE so the background is rebuilt.
"""

import io
import shutil
import sys
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio_ffmpeg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import xarray as xr
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image
from tqdm import tqdm

# ---------------------------------------------------------------- settings
RENDER_MODE = "full"            # "test" or "full"
TEST_DATES = ["2004-03-06", "2010-08-14", "2015-09-26", "2019-12-28"]

WORK = Path.home() / "CLIMATE ANIMATIONS" / "fires"                    # local APFS
WEEKLY = Path("/Volumes/CLIMATEDATA/fires/GFED5.1_weekly")
BASEMAP_DIR = Path("/Volumes/CLIMATEDATA/fires/basemap")
STAGING = WORK / "_staging"
TEST_DIR = WORK / "test_frames"
BG_CACHE = WORK / "fires_background_platecarree.png"
FINAL_VIDEO = Path("/Volumes/CLIMATEDATA/fires_2003_2022_platecarree.mp4")

WIDTH, HEIGHT, DPI = 1920, 1080, 100
FPS = 12                        # 987 weeks -> ~82 s (NASA original runs 1:24)

# Basemap colors (NASA SVS look)
OCEAN_COLOR = "#4d5868"
LAND_DARK, LAND_LIGHT = 0.10, 0.36   # gray range for shaded relief (0 = black, 1 = white)
LAND_FLAT = "#2e2e2e"               # used only if the relief raster is unavailable
BORDER_COLOR = "#d9d9d9"
RELIEF_URLS = [
    "https://naciscdn.org/naturalearth/50m/raster/GRAY_50M_SR.zip",
    "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/GRAY_50M_SR.zip",
]

# Data colors: green -> orange -> white, 0-20 g C m-2 week-1, saturating above 20
VMIN, VMAX = 0.0, 20.0
CMAP = LinearSegmentedColormap.from_list("nasa_fire", [
    (0.00, "#5e8c3a"),
    (0.25, "#a3b257"),
    (0.50, "#e5b260"),
    (0.75, "#f0d5a4"),
    (1.00, "#f3f1ec"),
])
CMAP.set_over("#f3f1ec")
NORM = Normalize(VMIN, VMAX)
ALPHA_MIN, ALPHA_FULL = 0.05, 1.0   # transparent below ALPHA_MIN, opaque above ALPHA_FULL
TEXT_COLOR = "#e6e6e6"
CREDIT = "Data: GFED5.1 (van der Werf et al.)  |  SDSU Climate Informatics Lab"


# ---------------------------------------------------------------- basemap
def load_relief() -> np.ndarray | None:
    """Natural Earth 1:50m gray shaded relief as uint8 (lat, lon), north up; None if unavailable."""
    BASEMAP_DIR.mkdir(parents=True, exist_ok=True)
    tif = BASEMAP_DIR / "GRAY_50M_SR.tif"
    if not tif.exists():
        for url in RELIEF_URLS:
            try:
                print(f"Downloading shaded relief: {url}")
                r = requests.get(url, timeout=120)
                r.raise_for_status()
                with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
                    member = next(m for m in zf.namelist()
                                  if m.lower().endswith(".tif") and not Path(m).name.startswith("._"))
                    tif.write_bytes(zf.read(member))
                break
            except Exception as e:  # noqa: BLE001
                print(f"  failed ({type(e).__name__}: {e})")
    if not tif.exists():
        print("Shaded relief unavailable; using flat land color.")
        return None
    Image.MAX_IMAGE_PIXELS = None
    arr = np.asarray(Image.open(tif))
    if arr.ndim == 3:
        arr = arr[..., :3].mean(axis=2)
    return arr[::2, ::2].astype(np.uint8)          # 10800x5400 -> 5400x2700 is plenty for 1920 px


def build_background() -> np.ndarray:
    """Render the static basemap once at output size; returns (HEIGHT, WIDTH, 3) float array."""
    if BG_CACHE.exists():
        return plt.imread(BG_CACHE)[..., :3]

    # Cartopy keeps Plate Carree at 2:1, so render at 2:1 (2x output width for
    # sharpness), then stretch to the 16:9 frame as in the NASA original.
    print("Building background ...")
    bg_w, bg_h = 2 * WIDTH, WIDTH
    fig = plt.figure(figsize=(bg_w / DPI, bg_h / DPI), dpi=DPI, facecolor=OCEAN_COLOR)
    ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.spines["geo"].set_visible(False)
    ax.set_facecolor(OCEAN_COLOR)

    relief = load_relief()
    if relief is not None:
        shade = LAND_DARK + (LAND_LIGHT - LAND_DARK) * (relief.astype(np.float32) / 255.0)
        ax.imshow(shade, cmap="gray", vmin=0, vmax=1, origin="upper",
                  extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), zorder=1)
    else:
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "50m"),
                       facecolor=LAND_FLAT, edgecolor="none", zorder=1)

    ax.add_feature(cfeature.NaturalEarthFeature("physical", "ocean", "50m"),
                   facecolor=OCEAN_COLOR, edgecolor="none", zorder=2)
    ax.add_feature(cfeature.NaturalEarthFeature("physical", "lakes", "50m"),
                   facecolor=OCEAN_COLOR, edgecolor="none", zorder=2)
    ax.add_feature(cfeature.NaturalEarthFeature("cultural", "admin_0_boundary_lines_land", "50m"),
                   facecolor="none", edgecolor=BORDER_COLOR, linewidth=0.9, alpha=0.6, zorder=3)

    fig.canvas.draw()
    big = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    if big.shape[:2] != (bg_h, bg_w):
        sys.exit(f"Background rendered at {big.shape[:2]}, expected {(bg_h, bg_w)}")
    bg = np.asarray(Image.fromarray(big).resize((WIDTH, HEIGHT), Image.LANCZOS))
    WORK.mkdir(parents=True, exist_ok=True)
    plt.imsave(BG_CACHE, bg)
    return bg.astype(np.float32) / 255.0


# ---------------------------------------------------------------- frames
def fire_rgba(field: np.ndarray) -> np.ndarray:
    """Map a weekly field to RGBA with an alpha ramp so near-zero cells stay transparent."""
    rgba = CMAP(NORM(field))
    rgba[..., 3] = np.clip((field - ALPHA_MIN) / (ALPHA_FULL - ALPHA_MIN), 0.0, 1.0)
    return rgba


def setup_figure(bg: np.ndarray, lat: np.ndarray, lon: np.ndarray):
    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI, facecolor=OCEAN_COLOR)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.imshow(bg, extent=[-180, 180, -90, 90], origin="upper", aspect="auto",
              interpolation="nearest", zorder=0)
    origin = "lower" if lat[0] < lat[-1] else "upper"
    fire_im = ax.imshow(np.zeros((len(lat), len(lon), 4)), extent=[-180, 180, -90, 90],
                        origin=origin, aspect="auto", interpolation="nearest", zorder=1)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

    # Legend block, bottom-left (NASA layout)
    cax = fig.add_axes([0.03, 0.185, 0.22, 0.022])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                      orientation="horizontal", extend="max", ticks=[0, 5, 10, 15, 20])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=TEXT_COLOR, labelsize=14, length=4)
    cb.set_label("grams of carbon per square meter per week", color=TEXT_COLOR,
                 fontsize=14, labelpad=6)
    xc = 0.03 + 0.22 / 2
    date_txt = fig.text(xc, 0.255, "", ha="center", va="bottom", color=TEXT_COLOR,
                        fontsize=22, fontweight="bold")
    fig.text(xc, 0.222, "Carbon Emissions", ha="center", va="bottom",
             color=TEXT_COLOR, fontsize=18)
    fig.text(0.99, 0.012, CREDIT, ha="right", va="bottom", color="#bdbdbd", fontsize=11)
    return fig, fire_im, date_txt


def week_index() -> list[tuple[Path, int, pd.Timestamp]]:
    files = sorted(p for p in WEEKLY.glob("GFED5.1_C_weekly_*.nc") if not p.name.startswith("._"))
    if not files:
        sys.exit(f"No weekly files in {WEEKLY}")
    index = []
    for p in files:
        with xr.open_dataset(p) as ds:
            index += [(p, k, t) for k, t in enumerate(pd.DatetimeIndex(ds["time"].values))]
    return index


# ---------------------------------------------------------------- main
def main() -> None:
    if not Path("/Volumes/CLIMATEDATA").is_dir():
        sys.exit("External drive /Volumes/CLIMATEDATA is not mounted.")
    WORK.mkdir(parents=True, exist_ok=True)
    STAGING.mkdir(parents=True, exist_ok=True)

    index = week_index()
    print(f"{len(index)} weekly frames: {index[0][2].date()} to {index[-1][2].date()}")
    with xr.open_dataset(index[0][0]) as ds:
        lat = ds["lat"].values
        lon = ds["lon"].values

    bg = build_background()
    fig, fire_im, date_txt = setup_figure(bg, lat, lon)

    if RENDER_MODE == "test":
        TEST_DIR.mkdir(parents=True, exist_ok=True)
        times = pd.DatetimeIndex([t for _, _, t in index])
        for d in TEST_DATES:
            i = int(np.abs(times - pd.Timestamp(d)).argmin())
            path, k, t = index[i]
            with xr.open_dataset(path) as ds:
                field = ds["C_weekly"].isel(time=k).values
            fire_im.set_data(fire_rgba(field))
            date_txt.set_text(t.strftime("%b  %Y"))
            out = TEST_DIR / f"fires_platecarree_{t.date()}.png"
            fig.savefig(out, dpi=DPI, facecolor=fig.get_facecolor())
            print(f"  wrote {out}")
        plt.close(fig)
        print("Review the PNGs, then set RENDER_MODE = 'full'.")
        return

    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    writer = FFMpegWriter(fps=FPS, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    tmp = STAGING / FINAL_VIDEO.name
    current_path, ds = None, None
    with writer.saving(fig, str(tmp), dpi=DPI):
        for path, k, t in tqdm(index, desc="Rendering", unit="frame"):
            if path != current_path:
                if ds is not None:
                    ds.close()
                ds, current_path = xr.open_dataset(path), path
            fire_im.set_data(fire_rgba(ds["C_weekly"].isel(time=k).values))
            date_txt.set_text(t.strftime("%b  %Y"))
            writer.grab_frame()
    ds.close()
    plt.close(fig)

    nframes, secs = imageio_ffmpeg.count_frames_and_secs(str(tmp))
    if nframes != len(index):
        sys.exit(f"Video has {nframes} frames, expected {len(index)}; left at {tmp}")
    shutil.move(str(tmp), str(FINAL_VIDEO))
    print(f"Done: {FINAL_VIDEO}  ({nframes} frames, {secs:.1f} s)")


if __name__ == "__main__":
    main()