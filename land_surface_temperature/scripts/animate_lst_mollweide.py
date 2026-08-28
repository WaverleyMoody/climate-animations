"""
SDSU Climate Informatics Lab 
San Diego State University 
by Waverley Moody 
Supervised by Distinguished Professor Samuel Shen 
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: animate_lst_mollweide.py

Description: Generates the Land Surface Temperature (Day) animation from NASA Terra/MODIS
MOD11C3 v6.1 (Feb 2000-Dec 2013, native 0.05 deg CMG grid), rendered in the
Mollweide projection. Unlike the other ERA5-based variables in this library,
this is a 167-frame real monthly time series, not a 48-frame climatological composite.

Note: Mollweide has a working inverse transform via PROJ/Cartopy, so this follows the
standard Cartopy GeoAxes rendering path (like Foucaut), not the manual forward-projection
pipeline required for Nicolosi.
"""

import os
import sys
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

IN_PATH = "/Volumes/CLIMATEDATA/modis_lst_monthly/lst_day_timeseries_2000_2013_native05deg_celsius.nc"
OUT_DIR = "/Volumes/CLIMATEDATA/modis_lst_monthly/animation_output/"
OUT_MP4 = os.path.join(OUT_DIR, "lst_day_mollweide_2000_2013.mp4")

VAR_NAME = "LST_Day"
VMIN, VMAX = -50, 55    # degrees C, physically sane global land LST range (converted from 220-330 K)

# Custom colormap: light blue -> dark blue -> purple -> bright purple -> red -> orange -> yellow
CUSTOM_COLORS = [
    "#ADD8E6",  # light blue
    "#00008B",  # dark blue
    "#800080",  # purple
    "#DA00FF",  # bright purple
    "#FF0000",  # red
    "#FFA500",  # orange
    "#FFFF00",  # yellow
]
CMAP = mcolors.LinearSegmentedColormap.from_list("lst_custom", CUSTOM_COLORS, N=256)

OCEAN_GRAY = "#808080"  # ocean / no-data background color

FPS = 6                 # slower than a 48-frame climatology loop since this has 167 real frames
DPI = 150


def _verify_readable(path):
    """Confirms the input time-series file is a real, uncorrupted read before animating -- project convention."""
    try:
        with xr.open_dataset(path) as ds:
            _ = ds[VAR_NAME].load()
        return True, None
    except Exception as e:
        return False, str(e)


def main():
    if not os.path.isfile(IN_PATH):
        print(f"ERROR: input time-series file not found: {IN_PATH}")
        print("Run concat_lst_timeseries.py first.")
        sys.exit(1)

    ok, reason = _verify_readable(IN_PATH)
    if not ok:
        print(f"ERROR: input file failed verification read: {reason}")
        sys.exit(1)

    os.makedirs(OUT_DIR, exist_ok=True)

    ds = xr.open_dataset(IN_PATH, chunks={"time": 1})  # chunked/lazy load -- avoids materializing the full array in memory at once
    da = ds[VAR_NAME]  # already in Celsius (see convert_lst_to_celsius.py)
    n_frames = da.sizes["time"]
    times = da["time"].values
    print(f"Loaded {n_frames} frames, {str(times.min())[:10]} to {str(times.max())[:10]}")

    lon = da["lon"].values
    lat = da["lat"].values

    fig = plt.figure(figsize=(12, 7))
    proj = ccrs.Mollweide(central_longitude=0)
    # Slightly taller/wider axes than the default single-axes fill, since removing the large
    # centered ax.set_title() frees up vertical space at the top for the map itself.
    ax = fig.add_axes([0.02, 0.03, 0.96, 0.90], projection=proj)
    ax.set_global()

    # Gray background fills the ocean and any masked/no-data land pixels (NaN cells are not
    # drawn by pcolormesh, so the axes facecolor shows through them)
    ax.set_facecolor(OCEAN_GRAY)

    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="black")
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="black", alpha=0.5)

    first_frame = da.isel(time=0).values  # already Celsius
    mesh = ax.pcolormesh(
        lon, lat, first_frame,
        transform=ccrs.PlateCarree(),
        cmap=CMAP, vmin=VMIN, vmax=VMAX,
        shading="auto",
    )
    cbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.05, shrink=0.7)
    cbar.set_label("Land Surface Temperature (\u00b0C)")

    # Static label (top-left) and changing date (top-right), in figure coordinates rather than
    # a single centered ax title -- smaller font frees up vertical space for a larger map.
    label_text = fig.text(
        0.03, 0.95, "MODIS Land Surface Temperature\n(Daytime)",
        fontsize=11, fontweight="bold", ha="left", va="top",
    )
    date_text = fig.text(
        0.97, 0.95, "", fontsize=11, fontweight="bold", ha="right", va="top",
    )

    def _format_date(dt64):
        # datetime64[M] truncation + .item() yields a plain datetime.date (day=1); strftime
        # then gives "April 2000" style formatting without a leading zero on the month.
        d = dt64.astype("datetime64[M]").item()
        return d.strftime("%B %Y")

    def update(frame_idx):
        data = da.isel(time=frame_idx).values  # already Celsius
        mesh.set_array(data.ravel())
        date_text.set_text(_format_date(times[frame_idx]))
        return mesh, label_text, date_text

    print(f"Rendering {n_frames} frames...")
    anim = animation.FuncAnimation(
        fig, update, frames=n_frames, interval=1000 // FPS, blit=False, cache_frame_data=False
    )

    writer = animation.FFMpegWriter(fps=FPS, bitrate=4000)
    anim.save(OUT_MP4, writer=writer, dpi=DPI)
    plt.close(fig)

    print(f"\nSaved animation to {OUT_MP4}")

    if os.path.isfile(OUT_MP4) and os.path.getsize(OUT_MP4) > 0:
        print("PASS -- output file exists and is non-empty.")
    else:
        print("FAIL -- output file missing or empty.")
        sys.exit(1)


if __name__ == "__main__":
    main()