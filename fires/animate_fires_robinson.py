"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: animate_fires_robinson.py

Description: Generates the fire carbon emissions time-series animation from
             GFED5.1 weekly totals (March 2003 - January 2022), rendered in the
             Robinson projection in the style of the NASA SVS "Fires - a global
             perspective" animation.

Note: For the Plate Carree, Foucaut, and Nicolosi projections, see the other
      scripts in the fires scripts folder.

Run inside the climate env:  conda activate climate
  1. RENDER_MODE = "test"  -> writes a few PNG frames for review (fast)
  2. RENDER_MODE = "full"  -> caffeinate -i python animate_fires_robinson.py
If you change any basemap colors or the layout, delete BG_CACHE and GEOM_CACHE.

"""

import io
import json
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
from matplotlib import patheffects
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
BG_CACHE = WORK / "fires_background_robinson.png"
GEOM_CACHE = WORK / "fires_background_robinson.json"
FINAL_VIDEO = Path("/Volumes/CLIMATEDATA/fires_2003_2022_robinson.mp4")

WIDTH, HEIGHT, DPI = 1920, 1080, 100
FPS = 12                        # 987 weeks -> ~82 s (NASA original runs 1:24)
PROJ = ccrs.Robinson(central_longitude=0)
MAP_MARGIN_X = 24               # px on each side; height follows the Robinson aspect

# Basemap colors (NASA SVS look)
FRAME_COLOR = "#1c1e22"             # outside the map outline
OCEAN_COLOR = "#4d5868"
LAND_DARK, LAND_LIGHT = 0.10, 0.36   # gray range for shaded relief (0 = black, 1 = white)
LAND_FLAT = "#2e2e2e"               # used only if the relief raster is unavailable
BORDER_COLOR = "#d9d9d9"
OUTLINE_COLOR = "#8b919b"
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
TEXT_HALO = [patheffects.withStroke(linewidth=3, foreground="#1c1e22")]  # legibility over land
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
    return arr[::2, ::2].astype(np.uint8)          # 10800x5400 -> 5400x2700


def map_rect() -> list[float]:
    """Axes rectangle (figure fraction) whose aspect matches the Robinson outline."""
    xr_ = PROJ.x_limits[1] - PROJ.x_limits[0]
    yr_ = PROJ.y_limits[1] - PROJ.y_limits[0]
    w_px = WIDTH - 2 * MAP_MARGIN_X
    h_px = w_px * yr_ / xr_
    return [MAP_MARGIN_X / WIDTH, (HEIGHT - h_px) / 2 / HEIGHT, w_px / WIDTH, h_px / HEIGHT]


def build_background() -> tuple[np.ndarray, dict]:
    """
    Render the static basemap once; returns ((HEIGHT, WIDTH, 3) float array, geometry),
    where geometry records the map's pixel box and projected limits for the lookup.
    """
    if BG_CACHE.exists() and GEOM_CACHE.exists():
        return plt.imread(BG_CACHE)[..., :3], json.loads(GEOM_CACHE.read_text())

    print("Building background ...")
    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI, facecolor=FRAME_COLOR)
    ax = fig.add_axes(map_rect(), projection=PROJ)
    ax.set_global()
    ax.set_facecolor(OCEAN_COLOR)
    ax.spines["geo"].set_edgecolor(OUTLINE_COLOR)
    ax.spines["geo"].set_linewidth(0.8)

    relief = load_relief()
    if relief is not None:
        shade = LAND_DARK + (LAND_LIGHT - LAND_DARK) * (relief.astype(np.float32) / 255.0)
        ax.imshow(shade, cmap="gray", vmin=0, vmax=1, origin="upper",
                  extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(),
                  regrid_shape=3000, zorder=1)
    else:
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "50m"),
                       facecolor=LAND_FLAT, edgecolor="none", zorder=1)

    ax.add_feature(cfeature.NaturalEarthFeature("physical", "ocean", "50m"),
                   facecolor=OCEAN_COLOR, edgecolor="none", zorder=2)
    ax.add_feature(cfeature.NaturalEarthFeature("physical", "lakes", "50m"),
                   facecolor=OCEAN_COLOR, edgecolor="none", zorder=2)
    ax.add_feature(cfeature.NaturalEarthFeature("cultural", "admin_0_boundary_lines_land", "50m"),
                   facecolor="none", edgecolor=BORDER_COLOR, linewidth=0.45, alpha=0.6, zorder=3)

    fig.canvas.draw()
    bg = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    pos = ax.get_position()                     # after Cartopy's aspect adjustment
    geom = {
        "x0": pos.x0 * WIDTH, "x1": pos.x1 * WIDTH,     # display px, origin bottom-left
        "y0": pos.y0 * HEIGHT, "y1": pos.y1 * HEIGHT,
        "xlim": list(ax.get_xlim()), "ylim": list(ax.get_ylim()),
    }
    plt.close(fig)
    if bg.shape[:2] != (HEIGHT, WIDTH):
        sys.exit(f"Background rendered at {bg.shape[:2]}, expected {(HEIGHT, WIDTH)}")
    WORK.mkdir(parents=True, exist_ok=True)
    plt.imsave(BG_CACHE, bg)
    GEOM_CACHE.write_text(json.dumps(geom, indent=2))
    return bg.astype(np.float32) / 255.0, geom


def build_lookup(geom: dict, lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    For every output pixel inside the Robinson outline, find the source grid cell.
    Returns (flat pixel indices, flat source-cell indices).
    """
    cols = np.arange(WIDTH) + 0.5
    rows_up = HEIGHT - (np.arange(HEIGHT) + 0.5)          # image row 0 is the top
    X = geom["xlim"][0] + (cols - geom["x0"]) / (geom["x1"] - geom["x0"]) * np.diff(geom["xlim"])[0]
    Y = geom["ylim"][0] + (rows_up - geom["y0"]) / (geom["y1"] - geom["y0"]) * np.diff(geom["ylim"])[0]
    XX, YY = np.meshgrid(X, Y)
    xx, yy = XX.ravel(), YY.ravel()

    ll = ccrs.PlateCarree().transform_points(PROJ, xx, yy)
    lonp, latp = ll[:, 0], ll[:, 1]
    back = PROJ.transform_points(ccrs.PlateCarree(), lonp, latp)    # round-trip check
    valid = (np.isfinite(lonp) & np.isfinite(latp)
             & (np.abs(back[:, 0] - xx) < 1000) & (np.abs(back[:, 1] - yy) < 1000))

    dlat, dlon = float(lat[1] - lat[0]), float(lon[1] - lon[0])
    ilat = np.clip(np.rint((latp[valid] - lat[0]) / dlat), 0, len(lat) - 1).astype(np.int64)
    ilon = np.clip(np.rint((lonp[valid] - lon[0]) / dlon), 0, len(lon) - 1).astype(np.int64)
    pix = np.flatnonzero(valid)
    print(f"Lookup: {pix.size:,} map pixels")
    return pix, ilat * len(lon) + ilon


# ---------------------------------------------------------------- frames
def fire_rgba_u8(field: np.ndarray) -> np.ndarray:
    """Weekly field -> flat (cells, 4) uint8 RGBA with an alpha ramp for near-zero cells."""
    rgba = CMAP(NORM(field), bytes=True)
    alpha = np.clip((field - ALPHA_MIN) / (ALPHA_FULL - ALPHA_MIN), 0.0, 1.0)
    rgba[..., 3] = (alpha * 255).astype(np.uint8)
    return rgba.reshape(-1, 4)


class Projector:
    """Turns a weekly field into a full-frame RGBA overlay using the pixel lookup."""

    def __init__(self, pix: np.ndarray, src: np.ndarray):
        self.pix, self.src = pix, src
        self.buf = np.zeros((HEIGHT * WIDTH, 4), dtype=np.uint8)

    def __call__(self, field: np.ndarray) -> np.ndarray:
        self.buf[self.pix] = fire_rgba_u8(field)[self.src]
        return self.buf.reshape(HEIGHT, WIDTH, 4)


def setup_figure(bg: np.ndarray):
    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI, facecolor=FRAME_COLOR)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.imshow(bg, origin="upper", aspect="auto", interpolation="nearest", zorder=0)
    fire_im = ax.imshow(np.zeros((HEIGHT, WIDTH, 4), dtype=np.uint8), origin="upper",
                        aspect="auto", interpolation="nearest", zorder=1)
    ax.set_xlim(-0.5, WIDTH - 0.5)
    ax.set_ylim(HEIGHT - 0.5, -0.5)

    # Legend block, bottom-left (NASA layout; sits over the open South Pacific corner)
    cax = fig.add_axes([0.03, 0.185, 0.22, 0.022])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                      orientation="horizontal", extend="max", ticks=[0, 5, 10, 15, 20])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=TEXT_COLOR, labelsize=14, length=4)
    cb.set_label("grams of carbon per square meter per week", color=TEXT_COLOR,
                 fontsize=14, labelpad=6, path_effects=TEXT_HALO)
    for lbl in cb.ax.get_xticklabels():
        lbl.set_path_effects(TEXT_HALO)
    xc = 0.03 + 0.22 / 2
    date_txt = fig.text(xc, 0.255, "", ha="center", va="bottom", color=TEXT_COLOR,
                        fontsize=22, fontweight="bold", path_effects=TEXT_HALO)
    fig.text(xc, 0.222, "Carbon Emissions", ha="center", va="bottom",
             color=TEXT_COLOR, fontsize=18, path_effects=TEXT_HALO)
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
        lat, lon = ds["lat"].values, ds["lon"].values

    bg, geom = build_background()
    project = Projector(*build_lookup(geom, lat, lon))
    fig, fire_im, date_txt = setup_figure(bg)

    if RENDER_MODE == "test":
        TEST_DIR.mkdir(parents=True, exist_ok=True)
        times = pd.DatetimeIndex([t for _, _, t in index])
        for d in TEST_DATES:
            i = int(np.abs(times - pd.Timestamp(d)).argmin())
            path, k, t = index[i]
            with xr.open_dataset(path) as ds:
                field = ds["C_weekly"].isel(time=k).values
            fire_im.set_data(project(field))
            date_txt.set_text(t.strftime("%b  %Y"))
            out = TEST_DIR / f"fires_robinson_{t.date()}.png"
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
            fire_im.set_data(project(ds["C_weekly"].isel(time=k).values))
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