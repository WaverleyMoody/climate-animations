"""
SDSU Climate Informatics Lab 
San Diego State University 
by Waverley Moody 
Supervised by Distinguished Professor Samuel Shen 
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: animate_biosphere_australia.py

Description: Generates the Global Biosphere animation from the merged land NDVI + ocean CHL dataset
    (2000-2010, 0.05 deg grid). Two acts: (1) a full 360-degree world tour that plays through the
    entire 2000-2010 monthly record while rotating, so the establishing shot shows real data
    evolving; (2) a zoom into a fixed, close-up view of the focus region
    (currently Australia) that replays the full 2000-2010 record a second time. Rendering is
    resumable -- each frame is written to its own PNG and already-rendered frames are skipped on
    rerun, so the render can be safely interrupted (e.g. closing the laptop lid) and picked up
    again later without losing progress. 

Note: Land NDVI and ocean CHL are two independently-sourced variables in the same merged NetCDF,
    each NaN over the other's domain, so they composite directly via pcolormesh layering with no
    additional masking logic needed. Residual CHL gaps (sea ice, quality-flag exclusions) are real, 
    expected features of ocean color data, not defects, and are rendered as a neutral ocean base color 
    rather than being filled in. Zoom is implemented via cartopy's NearsidePerspective projection 
    (a simulated satellite view), since Orthographic has no zoom control. 
"""

import glob
import os
import subprocess
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import dask
import imageio_ffmpeg
from pathlib import Path
from tqdm import tqdm

# --- Config -------------------------------------------------------------------

MERGED_DIR = "/Volumes/CLIMATEDATA/global_biosphere/merged"
FRAMES_DIR = "/Users/waverleymoody/CLIMATE ANIMATIONS/global_biosphere/frames"
OUTPUT_PATH = "/Users/waverleymoody/CLIMATE ANIMATIONS/global_biosphere/animate_biosphere_australia.mp4"

# How many consecutive frames each month is held for before advancing to the
# next month. Camera motion (rotation, zoom) still updates every single frame
# for smoothness -- only the DATA changes at this slower rate. Raise this to
# slow the animation down further; each unit here multiplies total runtime.
HOLD_FRAMES_PER_MONTH = 3

FPS = 8   # lowered from the original 12 -- combined with HOLD_FRAMES_PER_MONTH,
          # each month is now on screen for HOLD_FRAMES_PER_MONTH / FPS seconds

EASE_FRAMES = 40   # transition frames between the world tour and the focus-region zoom

ROTATION_CENTRAL_LATITUDE = 10.0   # slight tilt during the world tour, not a flat equatorial view

# Australia -- fixed, close-up focus region for the second act.
FOCUS_LON = 135.0
FOCUS_LAT = -25.0

# Simulated satellite altitude (meters). ~35.8 million m is geostationary
# altitude -- a whole-Earth-disk view similar to GOES full-disk imagery, used
# as the "world tour" framing. The focus-region zoom eases down to a much
# lower altitude for a close-up regional view. UNTESTED -- tune after first render.
WORLD_TOUR_SATELLITE_HEIGHT = 35_785_831
FOCUS_SATELLITE_HEIGHT = 4_000_000

# NDVI: brown-to-green vegetation ramp. Bare soil/sparse vegetation near -0.2,
# dense vegetation near 1.0. NaN (ocean) is fully transparent so the ocean
# layer underneath shows through.
NDVI_CMAP = plt.get_cmap("YlGn").copy()
NDVI_CMAP.set_bad(alpha=0)
NDVI_VMIN, NDVI_VMAX = -0.2, 1.0

# CHL: log-scaled, since chlorophyll concentration spans orders of magnitude
# (open-ocean "biological deserts" near 0.01 mg/m3 up to coastal blooms near
# 20+ mg/m3). NaN (land, or genuine data gaps -- see script docstring) is
# fully transparent, so the neutral ocean base feature color shows through
# rather than a harsh blank/black void.
CHL_CMAP = plt.get_cmap("nipy_spectral").copy()
CHL_CMAP.set_bad(alpha=0)
CHL_VMIN, CHL_VMAX = 0.01, 20.0

OCEAN_BASE_COLOR = "#0a1a2f"   # neutral dark navy for ocean pixels with no CHL data
LAND_BASE_COLOR = "#3a3a3a"    # neutral gray for land pixels with no NDVI data (rare)

dask.config.set(scheduler="synchronous")


# --- Data loading ---------------------------------------------------------------

def load_dataset() -> xr.Dataset:
    files = sorted(glob.glob(f"{MERGED_DIR}/merged_biosphere_*.nc"))
    return xr.open_mfdataset(files, combine="by_coords", chunks={"time": 1})


# --- Helpers --------------------------------------------------------------------

def _smoothstep(t: float) -> float:
    """Standard smoothstep easing (3t^2 - 2t^3), matching the SST animation's easing convention."""
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3 - 2 * t)


def _format_date(np_datetime64) -> str:
    """Platform-safe month/year formatting (avoids strftime's non-portable no-leading-zero codes)."""
    dt = np_datetime64.astype("datetime64[M]").astype(object)
    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    return f"{month_names[dt.month - 1]} {dt.year}"


def _frame_plan(n_months: int):
    """
    Builds the full frame-by-frame plan as a list of dicts, each describing one
    frame's camera position, satellite height, and data month index. Three acts:
      1. World tour: one full 360-degree rotation while playing through every
         month in the record (frame count = n_months * HOLD_FRAMES_PER_MONTH).
      2. Ease: camera pans and zooms from wherever the tour ended to the fixed,
         close-up focus-region view (EASE_FRAMES long, data held on month 0 so
         the second act starts the record over from the beginning).
      3. Focus zoom: fixed close-up view, replaying the full record again.
    """
    plan = []

    tour_frames = n_months * HOLD_FRAMES_PER_MONTH
    for f in range(tour_frames):
        lon = 180 - 360 * (f / tour_frames)
        month_index = f // HOLD_FRAMES_PER_MONTH
        plan.append({
            "lon": lon, "lat": ROTATION_CENTRAL_LATITUDE,
            "height": WORLD_TOUR_SATELLITE_HEIGHT,
            "month_index": month_index, "show_date": True,
        })

    tour_end_lon = plan[-1]["lon"]
    for f in range(EASE_FRAMES):
        t = _smoothstep(f / EASE_FRAMES)
        lon = tour_end_lon + t * (FOCUS_LON - tour_end_lon)
        lat = ROTATION_CENTRAL_LATITUDE + t * (FOCUS_LAT - ROTATION_CENTRAL_LATITUDE)
        height = WORLD_TOUR_SATELLITE_HEIGHT + t * (FOCUS_SATELLITE_HEIGHT - WORLD_TOUR_SATELLITE_HEIGHT)
        plan.append({
            "lon": lon, "lat": lat, "height": height,
            "month_index": 0, "show_date": False,
        })

    focus_frames = n_months * HOLD_FRAMES_PER_MONTH
    for f in range(focus_frames):
        month_index = f // HOLD_FRAMES_PER_MONTH
        plan.append({
            "lon": FOCUS_LON, "lat": FOCUS_LAT,
            "height": FOCUS_SATELLITE_HEIGHT,
            "month_index": month_index, "show_date": True,
        })

    return plan


# --- Rendering ---------------------------------------------------------------

def render_frame(fig, ds, spec: dict) -> None:
    fig.clf()

    ax = fig.add_subplot(
        1, 1, 1,
        projection=ccrs.NearsidePerspective(
            central_longitude=spec["lon"],
            central_latitude=spec["lat"],
            satellite_height=spec["height"],
        ),
    )
    ax.set_global()
    ax.add_feature(cfeature.OCEAN, facecolor=OCEAN_BASE_COLOR, zorder=0)
    ax.add_feature(cfeature.LAND, facecolor=LAND_BASE_COLOR, zorder=0)

    # Load just this one month's slice -- chunks={"time": 1} keeps this to a
    # single frame's worth of data in memory, not the full 131-month record.
    month_index = spec["month_index"]
    ndvi = ds["NDVI"].isel(time=month_index).load()
    chl = ds["CHL"].isel(time=month_index).load()

    ax.pcolormesh(
        ds["lon"], ds["lat"], chl,
        transform=ccrs.PlateCarree(),
        cmap=CHL_CMAP, norm=LogNorm(vmin=CHL_VMIN, vmax=CHL_VMAX),
        shading="auto", zorder=1,
    )
    ax.pcolormesh(
        ds["lon"], ds["lat"], ndvi,
        transform=ccrs.PlateCarree(),
        cmap=NDVI_CMAP, vmin=NDVI_VMIN, vmax=NDVI_VMAX,
        shading="auto", zorder=2,
    )

    ax.coastlines(linewidth=0.4, color="#888888", zorder=3)

    fig.text(0.02, 0.96, "Global Biosphere", fontsize=16, fontweight="bold",
              ha="left", va="top", color="white")

    if spec["show_date"]:
        date_label = _format_date(ds["time"].values[month_index])
        fig.text(0.98, 0.96, date_label, fontsize=14,
                  ha="right", va="top", color="white")

    fig.patch.set_facecolor("black")


def assemble_video(frames_dir: str, output_path: str, fps: int) -> None:
    """
    Stitches the rendered PNG frames into the final MP4 via ffmpeg (through
    imageio_ffmpeg's bundled binary). Folded in here from the former separate
    assemble_biosphere_video.py -- same ffmpeg invocation, just called
    automatically once main()'s render loop confirms every frame is present.
    """
    frame_files = sorted(glob.glob(f"{frames_dir}/frame_*.png"))
    if not frame_files:
        print(f"No frames found in {frames_dir} -- nothing to assemble.")
        return

    print(f"Assembling {len(frame_files)} frames into {output_path} at {fps} fps...")

    ffmpeg_path = imageio_ffmpeg.get_ffmpeg_exe()
    frame_pattern = str(Path(frames_dir) / "frame_%05d.png")

    subprocess.run([
        ffmpeg_path,
        "-y",
        "-framerate", str(fps),
        "-i", frame_pattern,
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        output_path,
    ], check=True)

    print(f"Saved video to {output_path}")


def main() -> None:
    Path(FRAMES_DIR).mkdir(parents=True, exist_ok=True)

    ds = load_dataset()
    n_months = ds.sizes["time"]
    plan = _frame_plan(n_months)
    total_frames = len(plan)

    fig = plt.figure(figsize=(10, 10), dpi=150)

    already_done = len([
        f for f in os.listdir(FRAMES_DIR)
        if f.startswith("frame_") and f.endswith(".png")
    ])
    if already_done > 0:
        print(f"Resuming: {already_done} frame(s) already rendered, skipping those.")

    for frame_num, spec in enumerate(tqdm(plan, desc="Rendering Global Biosphere frames")):
        frame_path = os.path.join(FRAMES_DIR, f"frame_{frame_num:05d}.png")
        if os.path.exists(frame_path):
            continue  # already rendered in a previous run -- safe to interrupt and resume

        render_frame(fig, ds, spec)
        fig.savefig(frame_path, facecolor="black")

    plt.close(fig)

    rendered_count = len([
        f for f in os.listdir(FRAMES_DIR)
        if f.startswith("frame_") and f.endswith(".png")
    ])
    print(f"All {total_frames} frames rendered to {FRAMES_DIR}")
    print(f"Estimated runtime at {FPS} fps: {total_frames / FPS:.0f} seconds")

    if rendered_count == total_frames:
        assemble_video(FRAMES_DIR, OUTPUT_PATH, FPS)
    else:
        print(f"WARNING: expected {total_frames} frames but found {rendered_count} -- "
              f"skipping assembly. Rerun this script to fill in the missing frame(s) first.")


if __name__ == "__main__":
    main()