"""
SDSU Climate Informatics Lab 
San Diego State University
By Waverley Moody 
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations Library
by Professor John Michael Wallace.

Script: render_sea_ice_motion.py

Description: Renders the "Sea ice motion" animation - daily Arctic sea ice concentration
             shading (NOAA/NSIDC G02202 V6) with IABP buoy tracks (red circles, 30-day
             green fading tails) drawn on top, in the sea ice grid's native NSIDC Polar
             Stereographic North projection.

Data: All inputs are provided in the GitHub Release `cryosphere_v1.0_data`:
        - sea_ice_concentration_1979_1989.zip, _1990_1999.zip, _2000_2009.zip, _2010_2022.zip
        - sea_ice_concentration_filled.zip
        - iabp_daily_positions_1979_2022_projected.parquet
      Unzip everything into one folder and set DATA_ROOT below to that folder. Expected layout:
        DATA_ROOT/
            sea_ice_concentration/<year>/sic_psn25_YYYYMMDD_*.nc
            sea_ice_concentration_filled/<year>/sic_psn25_YYYYMMDD_*.nc
            iabp_buoys/iabp_daily_positions_1979_2022_projected.parquet
      The finished video is written to DATA_ROOT/sea_ice_motion/.

Note: TEST_MODE renders a single year (2020) to validate colors, tail rendering, and
      frame timing before committing to the full 1979-2022 range (~16,000 frames).
      Set TEST_MODE = False and adjust YEAR_RANGE to run the full animation.

Note: The whole scene is rotated 90 degrees counterclockwise from the sea ice grid's
      native orientation, by shifting the projection's central meridian (lon_0) by -90
      degrees rather than rotating pixels/arrays - for a polar stereographic projection
      this is an exact rotation (no distortion), so ice shading, buoy tracks, coastlines,
      and the geographic labels all stay perfectly aligned with each other.

Note: Days where the source record has no valid ice concentration (SMMR/SSM/I outages,
      e.g. Jul-Aug 1984) are replaced by the linearly interpolated files in
      sea_ice_concentration_filled/. Those frames are labeled on screen.
"""

import re
import shutil
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio_ffmpeg
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pyproj
import xarray as xr

TEST_MODE = False
YEAR_RANGE = range(2020, 2021) if TEST_MODE else range(1979, 2023)
TAIL_DAYS = 30

# --- Paths: set DATA_ROOT to the folder where you unzipped the release assets ---
DATA_ROOT = Path("/Volumes/CLIMATEDATA")

ICE_DIR = DATA_ROOT / "sea_ice_concentration"
FILLED_DIR = DATA_ROOT / "sea_ice_concentration_filled"
BUOY_PARQUET = DATA_ROOT / "iabp_buoys" / "iabp_daily_positions_1979_2022_projected.parquet"

LOCAL_STAGING = Path.home() / "CLIMATE ANIMATIONS" / "staging" / "sea_ice_motion"  # local disk scratch
FINAL_DEST = DATA_ROOT / "sea_ice_motion"
OUT_NAME = "sea_ice_motion_test_2020.mp4" if TEST_MODE else "sea_ice_motion_1979_2022.mp4"

FPS = 52  # ~365 daily frames / 7 seconds per year
FIGSIZE = (10, 7)  # landscape - the grid's bounding box is wider than tall once rotated
DPI = 120

TAIL_COLOR = "#39FF14"  # bright/neon green, chosen for visibility against the blue-white ice
INTERP_LABEL = "Ice interpolated (satellite data gap)"

# Representative points for on-map geographic labels (lon, lat) - placed with a plain
# PlateCarree transform, so cartopy positions them correctly regardless of the rotation.
LABELS = {
    "ALASKA": (-155.0, 67.0),
    "GREENLAND": (-42.0, 72.0),
    "NORWAY": (11.0, 64.0),
    "SIBERIA": (105.0, 72.0),
}


def rotate_crs_90ccw(crs: pyproj.CRS) -> pyproj.CRS:
    """
    Build a new CRS identical to `crs` except with its central meridian (lon_0) shifted
    by -90 degrees. For a polar stereographic projection this is mathematically an exact
    90-degree counterclockwise rotation of the whole projected plane about the pole - not
    an approximation - so reusing already-projected x/y numbers under this new CRS (rather
    than re-projecting lat/lon through it) gives identical results, with no distortion.
    """
    proj4 = crs.to_proj4()
    match = re.search(r"\+lon_0=(-?\d+\.?\d*)", proj4)
    if not match:
        raise ValueError(f"Could not find +lon_0= in this CRS's proj4 string: {proj4}")
    old_lon0 = float(match.group(1))
    new_lon0 = old_lon0 - 90
    new_proj4 = re.sub(r"\+lon_0=-?\d+\.?\d*", f"+lon_0={new_lon0}", proj4)
    return pyproj.CRS.from_proj4(new_proj4)


def resolve_ice_path(path: Path):
    """Prefer the gap-filled copy of a day's ice file if one exists in FILLED_DIR."""
    filled = FILLED_DIR / path.parent.name / path.name
    if filled.exists() and not filled.name.startswith("._"):
        return filled, True
    return path, False


def _get_ice_var_name(ds: xr.Dataset) -> str:
    for candidate in ("cdr_seaice_conc", "seaice_conc_cdr", "nsidc_nt_seaice_conc"):
        if candidate in ds.variables:
            return candidate
    # fall back to the first 2D+ data variable that isn't a flag/qc field
    for name, var in ds.data_vars.items():
        if "qc" not in name.lower() and "flag" not in name.lower() and var.ndim >= 2:
            return name
    raise ValueError(f"Could not identify the ice concentration variable among: {list(ds.data_vars)}")


def load_ice_frame(path: Path, var_name_cache: dict):
    ds = xr.open_dataset(path)
    if "var_name" not in var_name_cache:
        var_name_cache["var_name"] = _get_ice_var_name(ds)
    var_name = var_name_cache["var_name"]

    conc = ds[var_name]
    if conc.ndim == 3:  # (time, y, x) with a length-1 time dim, common for these granules
        conc = conc.isel({conc.dims[0]: 0})

    values = conc.values.astype(float)

    # Defensive rescale: NSIDC CDR concentration should be a 0-1 fraction. If the file
    # instead stores 0-100 percent (varies by product version), rescale so the colormap
    # and vmin/vmax below stay correct either way.
    finite = values[np.isfinite(values)]
    if finite.size and np.nanmax(finite) > 1.5:
        values = values / 100.0

    x = ds["x"].values
    y = ds["y"].values
    ds.close()
    return values, x, y


def build_ice_colormap():
    cmap = matplotlib.colormaps["Blues_r"].copy()
    cmap.set_bad(color="#dcdcdc")  # land / missing shows as neutral gray, not white (=100% ice)
    return cmap


def main():
    if not ICE_DIR.exists():
        raise SystemExit(
            f"{ICE_DIR} not found - set DATA_ROOT at the top of this script to the folder "
            "where you unzipped the cryosphere_v1.0_data release assets."
        )
    if not BUOY_PARQUET.exists():
        raise SystemExit(
            f"{BUOY_PARQUET} not found - place iabp_daily_positions_1979_2022_projected.parquet "
            f"in {BUOY_PARQUET.parent}."
        )

    LOCAL_STAGING.mkdir(parents=True, exist_ok=True)
    FINAL_DEST.mkdir(parents=True, exist_ok=True)

    print("Collecting daily sea ice files for the requested year range...")
    ice_files = []
    for year in YEAR_RANGE:
        year_dir = ICE_DIR / str(year)
        if not year_dir.exists():
            print(f"  [WARN] {year_dir} does not exist, skipping {year}")
            continue
        ice_files.extend(sorted(
            p for p in year_dir.glob("sic_psn25_*_v*.nc") if not p.name.startswith("._")
        ))
    if not ice_files:
        raise FileNotFoundError(f"No sea ice files found for {list(YEAR_RANGE)} under {ICE_DIR}")
    n_filled = sum(resolve_ice_path(p)[1] for p in ice_files)
    print(f"  {len(ice_files)} daily ice files found ({n_filled} replaced by gap-filled copies)")
    if not FILLED_DIR.exists():
        print(
            f"  [WARN] {FILLED_DIR} not found - unzip sea_ice_concentration_filled.zip "
            "into DATA_ROOT, or gap days render blank"
        )

    print(f"Loading buoy positions from {BUOY_PARQUET} ...")
    buoys = pd.read_parquet(BUOY_PARQUET)
    buoys["date"] = pd.to_datetime(buoys["date"])
    buoys = buoys.sort_values(["BuoyID", "date"], kind="mergesort")
    # Rotate buoy positions 90 degrees CCW to match the rotated ice grid below:
    # (x, y) -> (-y, x) is the exact rotation matrix for +90 degrees counterclockwise.
    buoys["x_disp"] = -buoys["y_m"]
    buoys["y_disp"] = buoys["x_m"]
    print(f"  {len(buoys)} buoy-day rows loaded")

    # Establish the grid CRS/extent once from the first file.
    first_values, x_coords, y_coords = load_ice_frame(ice_files[0], {})
    extent = (float(x_coords.min()), float(x_coords.max()), float(y_coords.min()), float(y_coords.max()))
    ds0 = xr.open_dataset(ice_files[0])
    grid_mapping_var = None
    for candidate_name in ("crs", "polar_stereographic", "projection", "Polar_Stereographic_Grid"):
        if candidate_name in ds0.variables:
            grid_mapping_var = ds0.variables[candidate_name]
            break
    crs_wkt = grid_mapping_var.attrs.get("crs_wkt") or grid_mapping_var.attrs.get("spatial_ref")
    grid_crs = pyproj.CRS.from_wkt(crs_wkt)
    ds0.close()

    # Rotate the display 90 degrees CCW (see rotate_crs_90ccw docstring for why this is
    # exact for a polar stereographic grid). Everything plotted from here on - ice
    # shading, buoy tracks, coastlines, labels - uses `proj`, the rotated CRS, as its
    # transform/axes projection, so it all stays in registration.
    proj = ccrs.Projection(rotate_crs_90ccw(grid_crs))
    xmin, xmax, ymin, ymax = extent
    rotated_extent = (-ymax, -ymin, xmin, xmax)
    # cartopy's generic Projection wrapper derives .bounds from the CRS's own
    # "area of use" metadata (or leaves it None), not from the actual data grid - for
    # this CRS that gives a nonsensical small box (or nothing at all), which breaks
    # .boundary / .x_limits / .y_limits. Overwrite it with the (rotated) grid's real
    # x/y extent, which is the actual boundary that matters here. This is a documented
    # public attribute on cartopy's Projection class, not a private hack.
    proj.bounds = rotated_extent
    print(f"  using projection: {grid_crs.name} (rotated 90 deg CCW)")

    cmap = build_ice_colormap()
    var_name_cache = {}

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    # Reserve a strip of blank space above the map for the title/date, so that text
    # renders outside the map area instead of overlapping the ice/coastlines.
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.03)
    ax = plt.axes(projection=proj)

    writer = None

    try:
        for i, ice_path in enumerate(ice_files):
            date_str = ice_path.name.split("_")[2]  # sic_psn25_{YYYYMMDD}_...
            frame_date = pd.Timestamp(date_str)

            load_path, is_filled = resolve_ice_path(ice_path)
            values, x, y = load_ice_frame(load_path, var_name_cache)

            # Rotate the ice grid's own coordinate mesh the same way as the buoys
            # (x, y) -> (-y, x) - the data array itself (`values`) is untouched, only
            # the x/y location each cell is plotted at changes.
            x_mesh, y_mesh = np.meshgrid(x, y)
            x_mesh_disp = -y_mesh
            y_mesh_disp = x_mesh

            ax.clear()
            ax.set_extent(rotated_extent, crs=proj)
            ax.pcolormesh(
                x_mesh_disp, y_mesh_disp, values, transform=proj,
                cmap=cmap, vmin=0, vmax=1, shading="auto", zorder=0,
            )
            ax.add_feature(cfeature.COASTLINE, edgecolor="black", linewidth=0.5, zorder=2)
            ax.add_feature(cfeature.LAND, facecolor="#dcdcdc", zorder=1)

            gl = ax.gridlines(
                crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.6, color="gray",
                alpha=0.6, linestyle="--", zorder=2.5, x_inline=False, y_inline=True,
            )
            gl.xlocator = mticker.FixedLocator(range(-180, 181, 30))
            gl.ylocator = mticker.FixedLocator([60, 70, 80])
            gl.xlabel_style = {"size": 7, "color": "gray"}
            gl.ylabel_style = {"size": 7, "color": "gray"}

            for name, (lon, lat) in LABELS.items():
                ax.text(
                    lon, lat, name, transform=ccrs.PlateCarree(),
                    fontsize=9, fontweight="bold", color="black", ha="center", va="center",
                    zorder=5, path_effects=[matplotlib.patheffects.withStroke(
                        linewidth=2.5, foreground="white"
                    )],
                )

            window_start = frame_date - pd.Timedelta(days=TAIL_DAYS)
            window = buoys[(buoys["date"] > window_start) & (buoys["date"] <= frame_date)]

            for buoy_id, track in window.groupby("BuoyID"):
                track = track.sort_values("date")
                if len(track) < 2:
                    continue
                n = len(track)
                # Fading tail: oldest segment most transparent, most recent most opaque.
                for seg_i in range(n - 1):
                    alpha = 0.15 + 0.7 * (seg_i / max(n - 2, 1))
                    ax.plot(
                        track["x_disp"].iloc[seg_i:seg_i + 2],
                        track["y_disp"].iloc[seg_i:seg_i + 2],
                        transform=proj, color=TAIL_COLOR, linewidth=1.2, alpha=alpha, zorder=3,
                    )
                # Current position as a red circle.
                ax.plot(
                    track["x_disp"].iloc[-1], track["y_disp"].iloc[-1],
                    transform=proj, marker="o", markersize=3, color="red",
                    markeredgecolor="darkred", linestyle="none", zorder=4,
                )

            # Placed just above the axes (y > 1, va="bottom"), in the margin reserved by
            # fig.subplots_adjust above, so the title sits outside the map, not on it.
            ax.text(
                0.0, 1.05, "International Arctic Buoy Program", transform=ax.transAxes,
                fontsize=13, fontweight="bold", ha="left", va="bottom", zorder=6,
            )
            ax.text(
                1.0, 1.05, f"{frame_date.month}/{frame_date.day}/{frame_date.year}",
                transform=ax.transAxes, fontsize=13, ha="right", va="bottom", zorder=6,
            )
            if is_filled:
                # Small gray note under the date, still above the map, so viewers know
                # the ice field on this frame is interpolated rather than observed.
                ax.text(
                    1.0, 1.015, INTERP_LABEL, transform=ax.transAxes,
                    fontsize=8, color="0.4", ha="right", va="bottom", zorder=6,
                )

            fig.canvas.draw()
            frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
            frame = np.ascontiguousarray(frame)

            if writer is None:
                # Determine the real rendered pixel size from the canvas itself, rather
                # than trusting FIGSIZE*DPI arithmetic, so it can't drift from what
                # ffmpeg is actually told to expect.
                height_px, width_px = frame.shape[0], frame.shape[1]
                writer = imageio_ffmpeg.write_frames(
                    str(LOCAL_STAGING / OUT_NAME), size=(width_px, height_px), fps=FPS
                )
                writer.send(None)  # prime the generator

            writer.send(frame.tobytes())

            if (i + 1) % 30 == 0 or i == len(ice_files) - 1:
                print(f"  rendered {i + 1}/{len(ice_files)} frames")
    finally:
        if writer is not None:
            writer.close()
        plt.close(fig)

    local_out = LOCAL_STAGING / OUT_NAME
    final_out = FINAL_DEST / OUT_NAME
    shutil.move(str(local_out), str(final_out))
    print(f"Done: {final_out}")


if __name__ == "__main__":
    main()