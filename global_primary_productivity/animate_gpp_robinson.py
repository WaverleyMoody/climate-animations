"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: animate_gpp_robinson.py
Description: Generates the Gross Primary Productivity climatology
    animation from MOD17A2HGF (2000-2009), rendered in the Robinson
    projection.
Note: For the Plate Carrée projection, see animate_gpp_platecarree.py
    in the gpp scripts folder. Data excludes 2000-01-01 through
    2000-02-10 (early Terra commissioning gap) - the animation
    effectively starts 2000-02-18, 454 frames total.
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FuncAnimation
from matplotlib.colors import PowerNorm
from tqdm import tqdm

# ---- Config -------------------------------------------------------------
DATA_PATH = "/Volumes/CLIMATEDATA/gpp_2000_2009.nc"
OUTPUT_PATH = "/Volumes/CLIMATEDATA/gpp_2000_2009_robinson.mp4"
PROJECTION = ccrs.Robinson()
FPS = 4.7  # matches the original UW clip: ~454 frames / 98 seconds
CMAP = "YlGn"
VMIN, VMAX = 0, 0.12  # matches the observed value range from spot-checking
GAMMA = 0.6  # <1 darkens mid-range values, pulling more of the colormap's
             # saturated greens into the visible range instead of clustering
             # everything in the pale end
OCEAN_COLOR = "#0a1a3c"  # dark navy, matches the original clip's ocean
LAND_BASE_COLOR = "#f2f2ea"  # off-white/gray base for non-vegetated land

# GPP colorbar placement/style. Figure-fraction [left, bottom, width, height].
# Centered horizontally (left + width/2 == 0.5) and sized to leave room below
# for its tick labels and axis label, which are reserved via MAP_RECT below
# rather than left to matplotlib's default margins.
COLORBAR_RECT = [0.28, 0.085, 0.44, 0.035]
COLORBAR_LABEL = "GPP (kg C/m$^2$/8-day)"
COLORBAR_TICKS = [0, 0.03, 0.06, 0.09, 0.12]
COLORBAR_TEXT_COLOR = "black"

# Map axes placement, figure-fraction [left, bottom, width, height]. Bottom
# sits just above COLORBAR_RECT's top (0.085 + 0.035 = 0.12) with a thin gap,
# rather than the wide gap that was eating into the top margin (map height is
# fixed at 0.80, so a lower bottom = more headroom above for the title).
MAP_RECT = [0.01, 0.13, 0.98, 0.80]


def load_data():
    ds = xr.open_dataset(DATA_PATH, chunks={"time": 1})
    return ds["gpp"]


def render_frame(ax, da_frame):
    ax.clear()
    ax.set_global()
    ax.add_feature(cfeature.OCEAN, facecolor=OCEAN_COLOR, zorder=0)
    ax.add_feature(cfeature.LAND, facecolor=LAND_BASE_COLOR, zorder=1)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor="gray", zorder=3)

    # NaN-transparent overlay: ocean/no-data pixels show the base layers
    # underneath rather than a solid color, matching the original clip's
    # look of GPP-colored land against a plain navy ocean.
    im = ax.pcolormesh(
        da_frame.x, da_frame.y, da_frame.values,
        transform=ccrs.PlateCarree(),
        cmap=CMAP, norm=PowerNorm(gamma=GAMMA, vmin=VMIN, vmax=VMAX),
        zorder=2,
    )
    date_str = pd.Timestamp(da_frame.time.values).strftime("%b %-d, %Y")
    ax.set_title("Gross Primary Productivity", loc="left", fontsize=12)
    ax.set_title(date_str, loc="right", fontsize=12)
    return im


def _add_gpp_colorbar(fig, mappable, existing_cax=None):
    """
    Draws the GPP colorbar into a fixed figure-fraction axes. Removes any
    prior colorbar axes first, since only `ax` is cleared each frame (not
    the whole figure), and rebuilds fresh from this frame's own pcolormesh
    mappable so it can't drift out of sync with the GPP coloring.
    """
    if existing_cax is not None:
        fig.delaxes(existing_cax)

    cax = fig.add_axes(COLORBAR_RECT)
    cbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")

    cbar.set_ticks(COLORBAR_TICKS)
    cbar.set_label(COLORBAR_LABEL, color=COLORBAR_TEXT_COLOR, fontsize=11, labelpad=6)
    cbar.ax.tick_params(labelsize=9, colors=COLORBAR_TEXT_COLOR)
    cbar.outline.set_edgecolor(COLORBAR_TEXT_COLOR)
    cbar.outline.set_linewidth(0.6)

    return cax


def main():
    da = load_data()
    n_frames = da.sizes["time"]
    print(f"Rendering {n_frames} frames at {FPS} fps...")

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_axes(MAP_RECT, projection=PROJECTION)

    # cache_frame_data=False - OOM discipline on the 8GB Air; don't let
    # FuncAnimation hold references to every frame's data simultaneously
    progress = tqdm(total=n_frames)
    cbar_ax = {"cax": None}  # mutable container so update() can update it

    def update(i):
        frame = da.isel(time=i).load()
        im = render_frame(ax, frame)
        cbar_ax["cax"] = _add_gpp_colorbar(fig, im, existing_cax=cbar_ax["cax"])
        progress.update(1)
        return [im]

    anim = FuncAnimation(
        fig, update, frames=n_frames, cache_frame_data=False,
    )

    anim.save(OUTPUT_PATH, fps=FPS, writer="ffmpeg", dpi=150)
    progress.close()
    print(f"\nSaved animation to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()