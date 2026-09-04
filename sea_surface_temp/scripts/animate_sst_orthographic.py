"""
SDSU Climate Informatics Lab / San Diego State University / by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations Library by Professor John Michael Wallace.

Script: animate_sst_orthographic.py
Description: Generates the sea surface temperature (SST) daily time-series animation from JPL MUR25 (2016-2020),
             rendered in the orthographic projection as a spinning globe that gradually tilts to reveal the
             South Pole over the course of the animation.
Note: This is currently the only projection implemented for SST; other projections may be added later
      following the pattern used for climatology variables.
"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colorbar
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
IN_PATH = "/Volumes/CLIMATEDATA/sst_mur25/timeseries/sst_daily_timeseries_2016_2020_025deg_celsius.nc"
OUT_DIR = "/Users/waverleymoody/CLIMATE ANIMATIONS/sea_surface_temp/_local_output"
OUT_PATH = os.path.join(OUT_DIR, "sst_daily_2016_2020_orthographic.mp4")

VAR_NAME = "analysed_sst"

VMIN, VMAX = -2, 32           # degrees C
CMAP = "RdYlBu_r"

START_LAT = 20                # initial viewing latitude (tilt of the globe)
END_LAT = -55                 # final viewing latitude, tipped toward South Pole
TILT_START_FRAC = 0.6         # tilt spans the whole animation for a gradual reveal

ROTATIONS = 2                 # number of full 360-degree spins over the whole animation

FIGSIZE = (8, 8)
DPI = 150

TARGET_RUNTIME_SEC = 104      # 1:44

TITLE_TEXT = "Sea Surface Temperature"
# ---------------------------------------------------------------------------

os.makedirs(OUT_DIR, exist_ok=True)


def _verify_readable(path):
    try:
        with xr.open_dataset(path) as check:
            check[VAR_NAME].isel(time=0).load()
        return True
    except Exception as e:
        print(f"Verification failed: {e}")
        return False


def _format_date(dt64):
    """Format a numpy datetime64 as 'Jan 4, 2018' without relying on
    platform-specific strftime flags (%-d breaks on some systems)."""
    dt = dt64.astype("datetime64[s]").astype(object)  # -> python datetime
    return f"{dt.strftime('%b')} {dt.day}, {dt.year}"


def main():
    print("Opening dataset (chunked, lazy)...")
    ds = xr.open_dataset(IN_PATH, chunks={"time": 1})
    sst = ds[VAR_NAME]

    n_frames = sst.sizes["time"]
    times = ds.time.values
    lons = ds.lon.values
    lats = ds.lat.values

    fps = n_frames / TARGET_RUNTIME_SEC
    print(f"Total frames: {n_frames}")
    print(f"FPS set to {fps:.2f} to hit target runtime of {TARGET_RUNTIME_SEC}s")

    # ------------------------------------------------------------------
    # Precompute per-frame longitude (continuous spin) and latitude (tilt)
    # ------------------------------------------------------------------
    lon_per_frame = (ROTATIONS * 360.0) / n_frames
    central_lons = [(-i * lon_per_frame) % 360 for i in range(n_frames)]

    tilt_start_frame = int(n_frames * TILT_START_FRAC)
    central_lats = []
    for i in range(n_frames):
        if i < tilt_start_frame:
            central_lats.append(START_LAT)
        else:
            progress = (i - tilt_start_frame) / max(1, (n_frames - 1 - tilt_start_frame))
            eased = 0.5 - 0.5 * np.cos(np.pi * progress)  # smoothstep easing
            lat = START_LAT + eased * (END_LAT - START_LAT) 
            central_lats.append(lat)

    # ------------------------------------------------------------------
    # Figure setup: globe axes recreated per frame, colorbar fixed separately
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=FIGSIZE, facecolor="black")
    ax = None  # created fresh each frame

    norm = mcolors.Normalize(vmin=VMIN, vmax=VMAX)
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.03])  # fixed colorbar strip at bottom
    cbar = matplotlib.colorbar.ColorbarBase(
        cbar_ax, cmap=plt.get_cmap(CMAP), norm=norm,
        orientation="horizontal", extend="both",
    )
    cbar.set_label("Sea Surface Temperature (°C)", color="white")
    cbar.ax.xaxis.set_tick_params(color="white")
    plt.setp(plt.getp(cbar.ax.axes, "xticklabels"), color="white")

    # Bold title pinned top-left, date pinned top-right
    title_text = fig.text(
        0.06, 0.92, TITLE_TEXT, color="white", fontsize=13,
        ha="left", va="center", fontweight="bold",
    )
    date_text = fig.text(
        0.94, 0.92, "", color="white", fontsize=13,
        ha="right", va="center", fontweight="bold",
    )

    def update(frame_idx):
        nonlocal ax
        if ax is not None:
            ax.remove()

        proj = ccrs.Orthographic(
            central_longitude=central_lons[frame_idx],
            central_latitude=central_lats[frame_idx],
        )
        ax = fig.add_axes([0.05, 0.12, 0.9, 0.8], projection=proj)
        ax.set_global()
        ax.set_facecolor("black")
        ax.add_feature(cfeature.LAND, facecolor="#D3D3D3", zorder=1)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=2)

        data = sst.isel(time=frame_idx).load().values
        mesh = ax.pcolormesh(
            lons, lats, data,
            transform=ccrs.PlateCarree(),
            cmap=CMAP, norm=norm,
            shading="auto", zorder=0,
        )

        date_text.set_text(_format_date(times[frame_idx]))

        return mesh, ax, title_text, date_text

    print(f"Rendering animation to {OUT_PATH}...")
    anim = animation.FuncAnimation(
        fig, update, frames=n_frames,
        blit=False, cache_frame_data=False,
    )

    # ------------------------------------------------------------------
    # Progress bar tied to matplotlib's own save progress callback, so it
    # reflects frames actually encoded to disk, not just rendered in memory.
    # ------------------------------------------------------------------
    pbar = tqdm(total=n_frames, desc="Encoding frames", unit="frame")

    def _progress_callback(current_frame, total_frames):
        pbar.update(1)

    writer = animation.FFMpegWriter(fps=fps, bitrate=4000)
    anim.save(OUT_PATH, writer=writer, dpi=DPI, progress_callback=_progress_callback)
    pbar.close()
    plt.close(fig)

    print("Verifying output file...")
    if os.path.exists(OUT_PATH) and os.path.getsize(OUT_PATH) > 0:
        print(f"Done. Output saved: {OUT_PATH}")
    else:
        print("WARNING: output file missing or empty — check for errors above.")


if __name__ == "__main__":
    main()