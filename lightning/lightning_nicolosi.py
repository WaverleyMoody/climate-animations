"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: lightning_nicolosi.py
Description: Generates the lightning climatology animation from WWLLN/WGLC data (2010-2025), rendered as a double-hemisphere Nicolosi Globular projection. 365 calendar-day frames (Feb 29 excluded), each the multi-year daily-climatological mean, annualized (x365) to match the reference product's strokes km-2 yr-1 units.
Note: For the Plate Carrée, Robinson, and Foucaut projections, see the other scripts in the lightning scripts folder.

--- Implementation notes ---
PROJ's `nicol` operation has no inverse transform, so unlike Robinson and
Foucaut this cannot use a Cartopy GeoAxes at all. Instead: pyproj forward-
projects the data grid and Natural Earth coastline/border geometries by
hand for each hemisphere (Western: lon_0=-90, Eastern: lon_0=90), rendered
on plain Matplotlib axes with imageio_ffmpeg handling video export
directly (bypassing imageio.get_writer()'s plugin auto-resolution, which
can silently pick the wrong backend for .mp4 output).

Because annual lightning stroke density spans several orders of magnitude
globally (near-zero over oceans/poles vs. 100+ over the most active
tropical land regions), this uses a LOG-scale colormap (LogNorm passed to
pcolormesh) rather than the linear vmin/vmax used for the other variables'
Nicolosi scripts — a linear scale would make almost the entire globe look
blank except for a few saturated hotspots. Color gradient: white -> Tiffany
blue -> yellow -> orange -> red, matching the Plate Carrée/Robinson/Foucaut
lightning scripts exactly.
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
NC_PATH = "/Volumes/CLIMATEDATA/lightning/climatology/climatology_lightning_365frame.nc"
OUTPUT_DIR = Path("/Users/waverleymoody/Downloads/lightning_animation")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
OUT_MP4_PATH = str(OUTPUT_DIR / "lightning_climatology_nicolosi.mp4")
WEST_LON = -90.0
EAST_LON = 90.0
FPS = 8
DPI = 130

# Same log-scale color scale/gradient as the Plate Carrée, Robinson, and
# Foucaut lightning scripts, so this animation matches the others on the
# site. White at 0.003 strokes km-2 yr-1 is just a thin sliver, Tiffany
# blue holds as a plateau from ~0.0045 to 0.03, then continues through
# yellow -> orange -> red up to 30.
VMIN, VMAX = 0.003, 30
TICK_VALUES = [0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]

log_range = np.log10(VMAX / VMIN)
def log_pos(value):
    return np.log10(value / VMIN) / log_range

color_stops = [
    (log_pos(0.003),  '#ffffff'),  # white          (0.003, thin sliver)
    (log_pos(0.0045), '#0ABAB5'),  # Tiffany blue starts (0.0045)
    (log_pos(0.03),   '#0ABAB5'),  # Tiffany blue plateau ends (0.03)
    (log_pos(0.3),    '#ffff00'),  # yellow         (0.3)
    (log_pos(3.0),    '#ff8c00'),  # orange         (3)
    (log_pos(30.0),   '#ff0000'),  # red            (30)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_lightning', color_stops, N=256)
NORM = mcolors.LogNorm(vmin=VMIN, vmax=VMAX)


def find_lightning_var(ds):
    candidates = [v for v in ds.data_vars if v != "frame_label"]
    if "lightning_density_annual" in candidates:
        return "lightning_density_annual"
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

    # US state boundaries -- admin_1_states_provinces_lines covers every
    # country's first-level admin boundaries, so filter to the US only via
    # .records() (which carries attributes) rather than .geometries() (which
    # doesn't), matching the "US only, not every country's provinces" rule
    # used in the Robinson/Foucaut R scripts.
    shp_path_states = shpreader.natural_earth(
        resolution="110m", category="cultural", name="admin_1_states_provinces_lines")
    reader_states = shpreader.Reader(shp_path_states)
    all_state_records = list(reader_states.records())

    states = []
    matched_any = False
    for record in all_state_records:
        if record.attributes.get("ADM0_NAME") != "United States of America":
            continue
        matched_any = True
        clipped = record.geometry.intersection(hemi_bbox)
        if clipped.is_empty:
            continue
        proj_geom = shp_transform(lambda x, y: fwd.transform(x, y), clipped)
        states.append(proj_geom)

    # Defensive check: if the 'ADM0_NAME' field doesn't exist or has a
    # different value than expected in some other Natural Earth version,
    # the filter above would silently produce zero states with no error.
    # Fail loudly instead, with the actual available field names/values,
    # so this is easy to fix rather than a mysteriously missing layer.
    if not matched_any and all_state_records:
        sample_attrs = all_state_records[0].attributes
        print("WARNING: no state records matched ADM0_NAME=='United States "
              "of America'. Available attribute keys on this shapefile: "
              f"{list(sample_attrs.keys())}")
        print("Sample record attributes:", sample_attrs)

    return {
        "col_index": col_index,
        "X": X,
        "Y": Y,
        "boundary_x": boundary_x,
        "boundary_y": boundary_y,
        "coastlines": coastlines,
        "borders": borders,
        "states": states,
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
        elif geom.geom_type == "GeometryCollection":
            # Clipping a line against the hemisphere box can produce a
            # mixed-type collection rather than a clean LineString/
            # MultiLineString -- recurse into it so those pieces don't get
            # silently dropped (this was the states-not-showing bug: the
            # original version had no branch for this case at all).
            plot_coastlines(ax, list(geom.geoms), color=color, lw=lw)
        # Polygon/MultiPolygon pieces (shouldn't occur when clipping line
        # geometries, but skip silently rather than crash if they do --
        # they're not something we want outlined anyway).


def main():
    print(f"Opening {NC_PATH}")
    ds = xr.open_dataset(NC_PATH)
    var_name = find_lightning_var(ds)
    print(f"Using variable: '{var_name}'")

    da = ds[var_name]  # already annualized strokes km-2 yr-1, no conversion needed
    lat_name = "latitude" if "latitude" in da.coords else "lat"
    lon_name = "longitude" if "longitude" in da.coords else "lon"

    lats = da[lat_name].values
    lons = da[lon_name].values
    n_frames = da.sizes["frame"]
    print(f"{n_frames} frames, grid {lats.size} x {lons.size}")

    frame_labels = ds['frame_label'].values if 'frame_label' in ds else None
    if frame_labels is None:
        print("WARNING: no 'frame_label' variable found -- falling back to frame index in titles.")

    actual_min, actual_max = float(da.min()), float(da.max())
    print(f"Actual data range: {actual_min:.4f} to {actual_max:.2f} strokes km-2 yr-1 "
          f"(color scale fixed at {VMIN} to {VMAX})")
    if actual_max > VMAX:
        print("  NOTE: some values exceed the fixed color scale and will clip to the top color (extend='max').")

    print("Precomputing Western Hemisphere geometry...")
    west = build_hemisphere(lons, lats, WEST_LON)
    print(f"  West hemisphere: {len(west['states'])} state boundary pieces "
          f"(geom types: {set(g.geom_type for g in west['states'])})")
    print("Precomputing Eastern Hemisphere geometry...")
    east = build_hemisphere(lons, lats, EAST_LON)
    print(f"  East hemisphere: {len(east['states'])} state boundary pieces "
          f"(geom types: {set(g.geom_type for g in east['states'])})")

    def render_frame(i):
        frame_data = da.isel(frame=i).values
        # Log scale can't render exact zero (common over oceans/poles) --
        # clip up to VMIN so those cells render as the lowest color instead
        # of erroring or dropping out.
        frame_data = np.clip(frame_data, VMIN, None)
        label = str(frame_labels[i]) if frame_labels is not None else f"Frame {i}"

        fig, axes = plt.subplots(1, 2, figsize=(13, 7), constrained_layout=True)

        for ax, hemi in [(axes[0], west), (axes[1], east)]:
            data_hemi = frame_data[:, hemi["col_index"]]
            mesh = ax.pcolormesh(hemi["X"], hemi["Y"], data_hemi,
                                  shading="auto", cmap=CMAP, norm=NORM)
            plot_coastlines(ax, hemi["coastlines"], color="black", lw=1.0)
            plot_coastlines(ax, hemi["borders"], color="black", lw=0.8)
            plot_coastlines(ax, hemi["states"], color="black", lw=0.5)
            ax.plot(hemi["boundary_x"], hemi["boundary_y"], color="black", linewidth=1)
            ax.set_aspect("equal")
            ax.axis("off")

        # Header text matching the other lightning scripts' WWLLN/WGLC
        # style. Using ax.set_title(loc=...) instead of fig.text() here --
        # constrained_layout only reserves vertical space for real axes
        # artists like titles, not for arbitrary figure-coordinate text, so
        # fig.text() was overlapping the top of the circles.
        #
        # y is pinned explicitly (not left as the default "auto") because
        # auto-placement measures each axes' own content bounding box --
        # with axis("off") there's no spine to anchor to, and the two
        # hemisphere circles aren't pixel-identical in extent, so "auto"
        # put the two titles at very slightly different heights.
        TITLE_Y = 1.02
        TITLE_FONTSIZE = 16
        axes[0].set_title('WWLLN | WGLC', loc='left', y=TITLE_Y,
                            fontsize=TITLE_FONTSIZE, fontweight='bold')
        axes[1].set_title(f'{label}; 2010\u20132025 Climatology', loc='right', y=TITLE_Y,
                            fontsize=TITLE_FONTSIZE, fontweight='bold')

        cbar = fig.colorbar(mesh, ax=axes, orientation='horizontal',
                             pad=0.08, fraction=0.12, shrink=0.5,
                             aspect=18, extend='max')
        cbar.set_label('Lightning stroke density (strokes km$^{-2}$ yr$^{-1}$)', fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        # LogNorm colorbars auto-generate minor ticks (2,3,4...9 between
        # each decade) by default -- switch them off since we only want
        # our explicit TICK_VALUES showing.
        cbar.ax.minorticks_off()
        cbar.set_ticks(TICK_VALUES)
        cbar.set_ticklabels([str(t) for t in TICK_VALUES])

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