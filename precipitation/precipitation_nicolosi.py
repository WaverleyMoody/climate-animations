"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: precipitation_nicolosi.py
Description: Generates the precipitation climatology animation from ERA5 reanalysis
data (1979-2000), rendered as a double-hemisphere Nicolosi Globular projection. 
Weekly totals (cm water equivalent) computed server-side by the derived-era5-single-levels-daily-statistics 
CDS dataset, averaged into a 22-year climatological mean per 48-frame weekly cycle.
Note: For the Plate Carrée, Robinson, and Foucaut projections, see the other scripts in the precipitation scripts folder.
"""

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pyproj import Transformer
from shapely.geometry import box
from shapely.ops import transform as shp_transform
import cartopy.io.shapereader as shpreader
import imageio.v2 as imageio
import imageio_ffmpeg
import io
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ---- config ------------------------------------------------------------
NC_PATH = "/Volumes/CLIMATEDATA/precip/climatology/climatology_precip_48frame.nc"
OUTPUT_DIR = Path("/Users/waverleymoody/Downloads/precipitation_animation")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
OUT_MP4_PATH = str(OUTPUT_DIR / "precipitation_climatology_nicolosi.mp4")
WEST_LON = -90.0
EAST_LON = 90.0
FPS = 4
DPI = 130

# Same color scale/gradient as the Plate Carrée, Robinson, and Foucaut
# precipitation scripts, so this animation matches the others on the site.
# White at 0 cm blends into blue over just the first 0.3 cm, then continues
# through blue -> green -> yellow -> orange -> red -> purple up to 21 cm.
VMIN, VMAX = 0, 21
TICK_STEP = 3
color_stops = [
    (0 / VMAX, '#ffffff'),    # white           (0 cm)
    (0.3 / VMAX, '#d0e8f0'),  # pale blue       (0.3 cm)
    (2 / VMAX, '#7ab8e0'),    # light blue      (2 cm)
    (3 / VMAX, '#3060c0'),    # blue            (3 cm)
    (6 / VMAX, '#00994d'),    # green           (6 cm)
    (8 / VMAX, '#66cc00'),    # yellow-green    (8 cm)
    (9 / VMAX, '#ffff00'),    # yellow          (9 cm)
    (10 / VMAX, '#ff9900'),   # orange          (10 cm)
    (11 / VMAX, '#ff3300'),   # red-orange      (11 cm)
    (12 / VMAX, '#990000'),   # dark red        (12 cm)
    (21 / VMAX, '#660066'),   # purple          (21 cm)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_precip', color_stops, N=256)
BOUNDS = np.arange(VMIN, VMAX + 1, TICK_STEP)  # 0, 3, 6, ..., up to VMAX
if BOUNDS[-1] != VMAX:
    BOUNDS = np.append(BOUNDS, VMAX)  # always label the true top of the range


def find_precip_var(ds):
    candidates = [v for v in ds.data_vars if v != "frame_label"]
    if "tp_weekly_total_cm" in candidates:
        return "tp_weekly_total_cm"
    if len(candidates) == 1:
        return candidates[0]
    raise ValueError(f"Couldn't auto-pick the variable. Found: {candidates}.")


def build_hemisphere(lons, lats, central_lon):
    """Precompute everything that doesn't depend on frame data:
    the hemisphere mask, projected grid coords, boundary circle,
    and clipped/projected coastlines."""
    proj_str = f"+proj=nicol +lon_0={central_lon} +R=6371000"
    fwd = Transformer.from_crs("EPSG:4326", proj_str, always_xy=True)

    lon_diff = ((lons - central_lon + 180) % 360) - 180  # angular offset from center
    mask_1d = np.abs(lon_diff) <= 90
    if mask_1d.sum() < 2:
        raise RuntimeError(f"Hemisphere mask kept almost no columns for lon_0={central_lon}")

    # Sort the selected columns by their angular offset from the hemisphere
    # center, rather than trusting the original array's storage order. This
    # is what actually fixes the wraparound bug: lon_diff has no discontinuity
    # within a single +/-90deg window (its only jump is at +/-180 from center,
    # which is always excluded by the mask above), so this ordering is safe
    # for ANY central longitude and ANY input longitude convention -- unlike
    # a single global sort, which just relocates the seam onto whichever
    # hemisphere boundary happens to land on it.
    candidate_idx = np.where(mask_1d)[0]
    order = np.argsort(lon_diff[candidate_idx])
    col_index = candidate_idx[order]

    lons_hemi = lons[col_index]
    lon2d, lat2d = np.meshgrid(lons_hemi, lats)
    X, Y = fwd.transform(lon2d, lat2d)

    edge_lats = np.linspace(-90, 90, 400)
    bx1, by1 = fwd.transform(np.full_like(edge_lats, central_lon + 90), edge_lats)
    bx2, by2 = fwd.transform(np.full_like(edge_lats, central_lon - 90), edge_lats)
    boundary_x = np.concatenate([bx1, bx2[::-1]])
    boundary_y = np.concatenate([by1, by2[::-1]])

    shp_path = shpreader.natural_earth(resolution="110m", category="physical", name="coastline")
    reader = shpreader.Reader(shp_path)
    hemi_bbox = box(central_lon - 90, -90, central_lon + 90, 90)

    coastlines = []
    for record in reader.geometries():
        clipped = record.intersection(hemi_bbox)
        if clipped.is_empty:
            continue
        proj_geom = shp_transform(lambda x, y: fwd.transform(x, y), clipped)
        coastlines.append(proj_geom)

    # Country borders -- separate Natural Earth layer from coastlines
    # (admin_0_boundary_lines_land gives line geometries for country
    # borders directly, same structure as the coastline layer above, so
    # the same clip-then-reproject pipeline and plot_coastlines() helper
    # both work unchanged).
    shp_path_borders = shpreader.natural_earth(
        resolution="110m", category="cultural", name="admin_0_boundary_lines_land")
    reader_borders = shpreader.Reader(shp_path_borders)

    borders = []
    for record in reader_borders.geometries():
        clipped = record.intersection(hemi_bbox)
        if clipped.is_empty:
            continue
        proj_geom = shp_transform(lambda x, y: fwd.transform(x, y), clipped)
        borders.append(proj_geom)

    return {
        "col_index": col_index,
        "X": X,
        "Y": Y,
        "boundary_x": boundary_x,
        "boundary_y": boundary_y,
        "coastlines": coastlines,
        "borders": borders,
    }


def plot_coastlines(ax, coastlines, color="black", lw=0.5):
    for geom in coastlines:
        if geom.geom_type == "LineString":
            xs, ys = geom.xy
            ax.plot(xs, ys, color=color, linewidth=lw)
        elif geom.geom_type == "MultiLineString":
            for part in geom.geoms:
                xs, ys = part.xy
                ax.plot(xs, ys, color=color, linewidth=lw)


def main():
    print(f"Opening {NC_PATH}")
    ds = xr.open_dataset(NC_PATH)
    var_name = find_precip_var(ds)
    print(f"Using variable: '{var_name}'")

    da = ds[var_name]  # already in cm, no unit conversion needed
    lat_name = "latitude" if "latitude" in da.coords else "lat"
    lon_name = "longitude" if "longitude" in da.coords else "lon"

    lats = da[lat_name].values
    lons = da[lon_name].values
    n_frames = da.sizes["frame"]
    print(f"{n_frames} frames, grid {lats.size} x {lons.size}")

    frame_labels = ds['frame_label'].values if 'frame_label' in ds else None
    if frame_labels is None:
        print("WARNING: no 'frame_label' variable found -- falling back to frame index in titles.")

    # VMIN/VMAX are fixed (matching the other three precipitation scripts'
    # 0/21 scale), not derived from this data's actual range -- print a
    # heads-up if real values would exceed it, since pcolormesh silently
    # clips out-of-range values to the end color rather than erroring.
    actual_min, actual_max = float(da.min()), float(da.max())
    print(f"Actual data range: {actual_min:.2f} to {actual_max:.2f} cm "
          f"(color scale fixed at {VMIN} to {VMAX})")
    if actual_min < VMIN or actual_max > VMAX:
        print("  NOTE: some values fall outside the fixed color scale and will clip.")

    print("Precomputing Western Hemisphere geometry...")
    west = build_hemisphere(lons, lats, WEST_LON)
    print("Precomputing Eastern Hemisphere geometry...")
    east = build_hemisphere(lons, lats, EAST_LON)

    def render_frame(i):
        frame_data = da.isel(frame=i).values
        label = str(frame_labels[i]) if frame_labels is not None else f"Frame {i}"

        fig, axes = plt.subplots(1, 2, figsize=(13, 7), constrained_layout=True)

        for ax, hemi in [(axes[0], west), (axes[1], east)]:
            data_hemi = frame_data[:, hemi["col_index"]]
            mesh = ax.pcolormesh(hemi["X"], hemi["Y"], data_hemi,
                                  shading="auto", cmap=CMAP,
                                  vmin=VMIN, vmax=VMAX)
            plot_coastlines(ax, hemi["coastlines"], color="black", lw=1.0)
            plot_coastlines(ax, hemi["borders"], color="black", lw=0.8)
            ax.plot(hemi["boundary_x"], hemi["boundary_y"], color="black", linewidth=1)
            ax.set_aspect("equal")
            ax.axis("off")

        # Header text matching the other precipitation scripts' ERA-5/
        # Climate Reanalyzer style. Using ax.set_title(loc=...) instead of
        # fig.text() here -- constrained_layout only reserves vertical
        # space for real axes artists like titles, not for arbitrary
        # figure-coordinate text, so fig.text() was overlapping the top of
        # the circles.
        #
        # y is pinned explicitly (not left as the default "auto") because
        # auto-placement measures each axes' own content bounding box --
        # with axis("off") there's no spine to anchor to, and the two
        # hemisphere circles aren't pixel-identical in extent, so "auto"
        # put the two titles at very slightly different heights.
        TITLE_Y = 1.02
        TITLE_FONTSIZE = 16
        axes[0].set_title('ERA-5 | Climate Reanalyzer', loc='left', y=TITLE_Y,
                            fontsize=TITLE_FONTSIZE, fontweight='bold')
        axes[1].set_title(f'{label}; 1979\u20132000 Mean Weekly Total', loc='right', y=TITLE_Y,
                            fontsize=TITLE_FONTSIZE, fontweight='bold')

        cbar = fig.colorbar(mesh, ax=axes, orientation='horizontal',
                             pad=0.08, fraction=0.12, shrink=0.5,
                             aspect=18, extend='neither')
        cbar.set_label('Precipitation (cm)', fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_ticks(BOUNDS)
        cbar.set_ticklabels([str(t) for t in BOUNDS])

        buf = io.BytesIO()
        # NOTE: no bbox_inches="tight" here -- it makes each frame a
        # different pixel size depending on content extents, which
        # breaks the video encoder (needs identical, even dimensions
        # every frame).
        fig.savefig(buf, format="png", dpi=DPI)
        plt.close(fig)
        buf.seek(0)
        img = imageio.imread(buf)

        # drop alpha channel (write_frames wants RGB24) and force even
        # width/height (required by libx264)
        img = img[:, :, :3]
        h, w = img.shape[:2]
        img = img[: h - (h % 2), : w - (w % 2)]
        return np.ascontiguousarray(img)

    # Render frame 0 first to lock in the output size, then start the
    # ffmpeg process directly via imageio_ffmpeg -- this bypasses
    # imageio.get_writer()'s plugin auto-resolution, which on some
    # installs picks the tifffile backend for .mp4 output instead of
    # ffmpeg (that backend doesn't accept fps=, hence the TypeError).
    print("Rendering frame 0...")
    first_frame = render_frame(0)
    height, width = first_frame.shape[:2]

    ffmpeg_writer = imageio_ffmpeg.write_frames(
        OUT_MP4_PATH, size=(width, height), fps=FPS,
        codec="libx264", quality=8, macro_block_size=1,
    )
    ffmpeg_writer.send(None)  # seed/initialize the generator
    ffmpeg_writer.send(first_frame.tobytes())
    print(f"  frame 1/{n_frames} done")

    for i in range(1, n_frames):
        frame_img = render_frame(i)
        ffmpeg_writer.send(frame_img.tobytes())
        if i % 8 == 0:
            print(f"  frame {i + 1}/{n_frames} done")

    ffmpeg_writer.close()
    print(f"Saved: {OUT_MP4_PATH}")


if __name__ == "__main__":
    main()