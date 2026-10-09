"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: saharan_dust_platecarree.py
Description: Generates the Saharan dust optical depth animation from NASA
    MERRA-2 hourly DUEXTTAU (2020-06-01 to 2020-06-30), rendered in the
    PlateCarree projection.
Note: For the data download, see download_merra2_dust_2020_06.py in the
    saharan_dust scripts folder.

Run:
    conda activate climate
    # 1) Set PREVIEW_ONLY = True, check the PNG
    # 2) Set PREVIEW_ONLY = False for the full render
    caffeinate -i python saharan_dust_platecarree.py
"""

import shutil
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import dask
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
VARIABLE = "DUEXTTAU"
DATA_PATH = Path("/Volumes/CLIMATEDATA/dust/merra2_2020_06/MERRA2_DUEXTTAU_hourly_2020-06.nc")

LOCAL_DIR = Path.home() / "CLIMATE ANIMATIONS" / "saharan_dust" / "_render"   # local APFS
FINAL_DIR = Path("/Volumes/CLIMATEDATA/dust/animations")                     # exFAT
OUT_NAME = "saharan_dust_june2020_platecarree.mp4"

PREVIEW_ONLY = False                       # True: save one PNG and stop
PREVIEW_TIME = "2020-06-06T10:30"         # matches the reference frame

# Map extent (lon_min, lon_max, lat_min, lat_max) — limited by the data subset
MAP_EXTENT = [-128.0, 42.0, 0.0, 55.0]

# Colorbar levels from the original (log-spaced, extended both ends)
LEVELS = [0.10, 0.15, 0.21, 0.31, 0.45, 0.66, 0.97, 1.40, 2.10, 3.00]
CBAR_LABEL = "Dust optical depth (550 nm)"

BG_COLOR = "#fdf5e8"                      # cream ocean/land background
CMAP_ANCHORS = ["#fbe9cf", "#f6cf9c", "#eeac72", "#e2864f",
                "#cd5f34", "#ad3c21", "#842312", "#5a0e06"]
OVER_COLOR = "#3d0703"

SHOW_US_CITIES = True
CITY_MIN_POP = 50_000

TITLE = "Saharan Dust, June 2020"
SUBTITLE = "NASA MERRA-2 hourly dust optical depth (550 nm)"

LON_TICKS = list(range(-120, 41, 20))     # 120°W ... 40°E
LAT_TICKS = list(range(0, 51, 10))        # 0° ... 50°N
SHOW_GRIDLINES = True                     # False: labels only, no lines

FPS = 24                                  # 720 frames -> 30 s, same as original
FIGSIZE = (14, 5.2)                       # 1400 x 520 px at DPI 100
DPI = 100
BITRATE = 6000

dask.config.set(scheduler="synchronous")
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_cmap():
    smooth = LinearSegmentedColormap.from_list("dust", CMAP_ANCHORS)
    n_bins = len(LEVELS) - 1
    cmap = ListedColormap(smooth(np.linspace(0, 1, n_bins)), name="dust_binned")
    cmap.set_under(BG_COLOR)
    cmap.set_over(OVER_COLOR)
    cmap.set_bad(BG_COLOR)
    return cmap, BoundaryNorm(LEVELS, n_bins)


def _us_city_points():
    try:
        path = shpreader.natural_earth(resolution="10m", category="cultural",
                                       name="populated_places")
    except Exception as exc:
        print(f"City layer skipped ({exc})")
        return np.array([]), np.array([])
    lons, lats = [], []
    for rec in shpreader.Reader(path).records():
        a = {k.upper(): v for k, v in rec.attributes.items()}
        if a.get("ADM0_A3") == "USA" and (a.get("POP_MAX") or 0) >= CITY_MIN_POP:
            lons.append(rec.geometry.x)
            lats.append(rec.geometry.y)
    return np.array(lons), np.array(lats)


def _stamp(t):
    return pd.Timestamp(t).strftime("%Y-%m-%d %H:%Mz")


def _build_figure(da):
    lon, lat = da.lon.values, da.lat.values
    dlon, dlat = lon[1] - lon[0], lat[1] - lat[0]
    img_extent = [lon[0] - dlon / 2, lon[-1] + dlon / 2,
                  lat[0] - dlat / 2, lat[-1] + dlat / 2]

    cmap, norm = _build_cmap()
    proj = ccrs.PlateCarree()

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI, facecolor="white")
    ax = fig.add_axes([0.045, 0.08, 0.83, 0.72], projection=proj)
    ax.set_extent(MAP_EXTENT, crs=proj)
    ax.set_facecolor(BG_COLOR)
    ax.spines["geo"].set_visible(False)

    im = ax.imshow(np.full((lat.size, lon.size), np.nan), origin="lower",
                   extent=img_extent, transform=proj, cmap=cmap, norm=norm,
                   interpolation="bilinear", zorder=1)

    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.4,
                   edgecolor="#555555", zorder=2)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.3,
                   edgecolor="#777777", zorder=2)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.3,
                   edgecolor="#777777", zorder=2)

    if SHOW_US_CITIES:
        clon, clat = _us_city_points()
        if clon.size:
            ax.scatter(clon, clat, s=1.5, c="black", linewidths=0,
                       transform=proj, zorder=3)

    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=0.3,
                      color="#999999", alpha=0.6, linestyle="--", zorder=2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = SHOW_GRIDLINES
    gl.ylines = SHOW_GRIDLINES
    gl.xlocator = mticker.FixedLocator(LON_TICKS)
    gl.ylocator = mticker.FixedLocator(LAT_TICKS)
    gl.xlabel_style = {"size": 8, "color": "#333333"}
    gl.ylabel_style = {"size": 8, "color": "#333333"}

    fig.text(0.46, 0.925, TITLE, ha="center", va="bottom",
             fontsize=15, fontweight="bold")
    fig.text(0.46, 0.875, SUBTITLE, ha="center", va="bottom",
             fontsize=10, color="#444444")

    cax = fig.add_axes([0.915, 0.12, 0.012, 0.64])
    cb = fig.colorbar(im, cax=cax, extend="both", ticks=LEVELS)
    cb.ax.set_yticklabels([f"{v:.2f}" for v in LEVELS], fontsize=8)
    cb.set_label(CBAR_LABEL, fontsize=9)
    cb.ax.yaxis.set_label_position("left")

    stamp = fig.text(0.985, 0.025, "", ha="right", va="bottom", fontsize=9)
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
    da = ds[VARIABLE].sel(lon=slice(MAP_EXTENT[0] - 1, MAP_EXTENT[1] + 1),
                          lat=slice(MAP_EXTENT[2], MAP_EXTENT[3] + 1))
    times = da.time.values
    n = times.size
    print(f"Loaded {VARIABLE}: {n} hourly frames, {da.sizes['lat']} x {da.sizes['lon']}")

    fig, im, stamp = _build_figure(da)

    if PREVIEW_ONLY:
        t = da.sel(time=PREVIEW_TIME, method="nearest")
        im.set_data(t.values)
        stamp.set_text(_stamp(t.time.values))
        out = LOCAL_DIR / "preview_frame.png"
        fig.savefig(out, dpi=DPI)
        print(f"Preview saved: {out}")
        return

    local_mp4 = LOCAL_DIR / OUT_NAME
    writer = FFMpegWriter(fps=FPS, codec="libx264", bitrate=BITRATE,
                          extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, str(local_mp4), dpi=DPI):
        for i in range(n):
            im.set_data(da.isel(time=i).values)
            stamp.set_text(_stamp(times[i]))
            writer.grab_frame()
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