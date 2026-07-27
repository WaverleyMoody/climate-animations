"""
animate_nicolosi_wind.py

Full 48-frame double-hemisphere Nicolosi Globular animation for 10m wind
speed, with mean wind vectors overlaid.

Same base pipeline as animate_nicolosi_2m_temp.py / animate_nicolosi_SLP.py
(hemisphere masking with the wraparound-ordering fix, Shapely clip-then-
reproject coastlines, imageio_ffmpeg direct-write export, pinned-title-y
header). This file adds one new piece: cartopy's ax.quiver(transform=...)
normally rotates wind vectors automatically when the display projection
differs from the data's native lon/lat grid. Since Nicolosi is rendered by
hand (no GeoAxes), that rotation has to be replicated manually here via a
local finite-difference Jacobian of the forward projection at each arrow
location -- see transform_vectors() below.
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

matplotlib.rcParams['font.family'] = 'Arial'

# ---- config ------------------------------------------------------------
NC_PATH = "/Users/waverleymoody/Downloads/climate_data_by_variable/wind_climatology.nc"
OUT_MP4_PATH = "/Users/waverleymoody/Desktop/nicolosi_wind.mp4"
WEST_LON = -90.0
EAST_LON = 90.0
FPS = 4
DPI = 130

# Same fixed color scale and custom colormap as animate_wind.py.
VMIN, VMAX = 0, 16
colors = [
    (0.00, '#ffffff'),  # white         (0 m/s)
    (0.20, '#d0e8f0'),  # light blue    (3 m/s)
    (0.35, "#469a46"),  # light green   (5.5 m/s)
    (0.44, '#ffff00'),  # yellow        (7 m/s)
    (0.50, '#ff6600'),  # orange        (8 m/s)
    (0.60, '#ff3300'),  # darker orange (9.5 m/s)
    (0.70, '#ff0000'),  # light red     (11 m/s)
    (0.80, '#cc0000'),  # red           (13 m/s)
    (0.90, '#990000'),  # dark red      (14.5 m/s)
    (1.00, '#660000'),  # darkest red   (16 m/s)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_wind', colors, N=35)

# Quiver settings. QUIVER_STEP subsamples the grid the same way
# animate_wind.py's step=30 does -- and needs to actually match: the real
# ERA5 grid is 0.25deg resolution (721 x 1440), same as the other scripts
# in this project, so QUIVER_STEP=30 gives the same ~7.5deg arrow spacing
# animate_wind.py uses.
#
# ARROW_SCALE is NOT copied from animate_wind.py's scale=300, because that
# value is tuned for PlateCarree's degree-based coordinate space --
# Nicolosi's X/Y here are in meters (radius ~6.37e6 m), a totally
# different numeric scale, so it needs its own value. With scale_units='xy'
# and scale=1/ARROW_SCALE, arrow length in meters = wind_speed * ARROW_SCALE.
#
# ARROW_WIDTH is the shaft width IN METERS (not axes-fraction). This
# matters: quiver has two independent unit systems -- scale_units controls
# arrow LENGTH, but width/headwidth/headlength are controlled separately by
# `units`, which defaults to 'width' (a fraction of the axes box, totally
# disconnected from the data coordinate system). With the hemisphere
# spanning ~2e7 meters, a 'width'-fraction-sized head ends up physically
# bigger than the entire arrow length in data units, swallowing the shaft
# completely -- no amount of increasing ARROW_SCALE fixes that, since it's
# a units mismatch, not a magnitude problem. Setting units='xy' below makes
# width/head sizing scale consistently with arrow length in meters, same as
# scale_units already did for length.
QUIVER_STEP = 24
ARROW_SCALE = 90000
ARROW_WIDTH = 22000  # shaft width in meters



def transform_vectors(fwd, lons, lats, u, v, eps=0.05):
    """Rotate wind vectors (u=eastward, v=northward, in m/s) into the
    projected map's local x/y directions, via a finite-difference Jacobian
    of the forward projection. This replicates what cartopy's
    ax.quiver(transform=...) does automatically for a real GeoAxes."""
    R = 6371000.0
    x0, y0 = fwd.transform(lons, lats)
    xe, ye = fwd.transform(lons + eps, lats)
    xn, yn = fwd.transform(lons, lats + eps)
    east_dist = eps * np.pi / 180 * R * np.cos(np.radians(lats))
    north_dist = eps * np.pi / 180 * R
    # guard against the pole, where east_dist -> 0
    east_dist = np.where(np.abs(east_dist) < 1e-6, 1e-6, east_dist)
    dxde = (xe - x0) / east_dist
    dyde = (ye - y0) / east_dist
    dxdn = (xn - x0) / north_dist
    dydn = (yn - y0) / north_dist
    u_map = dxde * u + dxdn * v
    v_map = dyde * u + dydn * v
    return u_map, v_map


def find_wind_vars(ds):
    required = ["wind_speed", "u10", "v10"]
    missing = [v for v in required if v not in ds.data_vars]
    if missing:
        raise ValueError(f"Missing expected wind variable(s): {missing}. "
                          f"Found: {list(ds.data_vars)}")
    return required


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

    # Quiver geometry: subsample the (already hemisphere-restricted,
    # correctly-ordered) row/col indices, project those points, and
    # precompute the local rotation Jacobian at each -- none of this
    # depends on frame data, only on the grid and the projection.
    quiver_rows = np.arange(0, len(lats), QUIVER_STEP)
    quiver_cols = col_index[::QUIVER_STEP]

    qlon2d, qlat2d = np.meshgrid(lons[quiver_cols], lats[quiver_rows])
    quiver_X, quiver_Y = fwd.transform(qlon2d, qlat2d)

    R = 6371000.0
    eps = 0.05
    qx0, qy0 = fwd.transform(qlon2d, qlat2d)
    qxe, qye = fwd.transform(qlon2d + eps, qlat2d)
    qxn, qyn = fwd.transform(qlon2d, qlat2d + eps)
    east_dist = eps * np.pi / 180 * R * np.cos(np.radians(qlat2d))
    north_dist = eps * np.pi / 180 * R
    east_dist = np.where(np.abs(east_dist) < 1e-6, 1e-6, east_dist)
    quiver_dxde = (qxe - qx0) / east_dist
    quiver_dyde = (qye - qy0) / east_dist
    quiver_dxdn = (qxn - qx0) / north_dist
    quiver_dydn = (qyn - qy0) / north_dist

    return {
        "col_index": col_index,
        "X": X,
        "Y": Y,
        "boundary_x": boundary_x,
        "boundary_y": boundary_y,
        "coastlines": coastlines,
        "borders": borders,
        "quiver_rows": quiver_rows,
        "quiver_cols": quiver_cols,
        "quiver_X": quiver_X,
        "quiver_Y": quiver_Y,
        "quiver_dxde": quiver_dxde,
        "quiver_dyde": quiver_dyde,
        "quiver_dxdn": quiver_dxdn,
        "quiver_dydn": quiver_dydn,
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
    find_wind_vars(ds)  # raises if wind_speed/u10/v10 are missing
    print("Using variables: 'wind_speed', 'u10', 'v10'")

    speed_da = ds["wind_speed"]  # already m/s, no unit conversion needed
    u_da = ds["u10"]
    v_da = ds["v10"]
    lat_name = "latitude" if "latitude" in speed_da.coords else "lat"
    lon_name = "longitude" if "longitude" in speed_da.coords else "lon"

    lats = speed_da[lat_name].values
    lons = speed_da[lon_name].values
    n_frames = speed_da.sizes["frame"]
    print(f"{n_frames} frames, grid {lats.size} x {lons.size}")

    frame_labels = ds['frame_label'].values if 'frame_label' in ds else None
    if frame_labels is None:
        print("WARNING: no 'frame_label' variable found -- falling back to frame index in titles.")

    # VMIN/VMAX are fixed (matching animate_wind.py's 0/16 scale), not
    # derived from this data's actual range -- print a heads-up if real
    # values would exceed it, since pcolormesh silently clips out-of-range
    # values to the end color rather than erroring.
    actual_min, actual_max = float(speed_da.min()), float(speed_da.max())
    print(f"Actual wind speed range: {actual_min:.1f} to {actual_max:.1f} m/s "
          f"(color scale fixed at {VMIN} to {VMAX})")
    if actual_min < VMIN or actual_max > VMAX:
        print("  NOTE: some values fall outside the fixed color scale and will clip.")

    print("Precomputing Western Hemisphere geometry...")
    west = build_hemisphere(lons, lats, WEST_LON)
    print("Precomputing Eastern Hemisphere geometry...")
    east = build_hemisphere(lons, lats, EAST_LON)

    def render_frame(i):
        speed_frame = speed_da.isel(frame=i).values
        u_frame = u_da.isel(frame=i).values
        v_frame = v_da.isel(frame=i).values
        label = str(frame_labels[i]) if frame_labels is not None else f"Frame {i}"

        fig, axes = plt.subplots(1, 2, figsize=(13, 7), constrained_layout=True)

        for ax, hemi in [(axes[0], west), (axes[1], east)]:
            data_hemi = speed_frame[:, hemi["col_index"]]
            mesh = ax.pcolormesh(hemi["X"], hemi["Y"], data_hemi,
                                  shading="auto", cmap=CMAP,
                                  vmin=VMIN, vmax=VMAX)
            plot_coastlines(ax, hemi["coastlines"], color="black", lw=1.0)
            plot_coastlines(ax, hemi["borders"], color="black", lw=0.8)
            ax.plot(hemi["boundary_x"], hemi["boundary_y"], color="black", linewidth=1)

            # Quiver overlay: pull u/v at the precomputed subsampled grid
            # points, then rotate them into map-local x/y using the
            # Jacobian computed once in build_hemisphere.
            u_q = u_frame[np.ix_(hemi["quiver_rows"], hemi["quiver_cols"])]
            v_q = v_frame[np.ix_(hemi["quiver_rows"], hemi["quiver_cols"])]
            u_map = hemi["quiver_dxde"] * u_q + hemi["quiver_dxdn"] * v_q
            v_map = hemi["quiver_dyde"] * u_q + hemi["quiver_dydn"] * v_q
            ax.quiver(hemi["quiver_X"], hemi["quiver_Y"], u_map, v_map,
                       color='black', scale=1 / ARROW_SCALE, scale_units='xy',
                       units='xy', width=ARROW_WIDTH,
                       headwidth=3, headlength=4, headaxislength=3.5)

            ax.set_aspect("equal")
            ax.axis("off")

        # Header text matching animate_2m_temp.py's ERA-5/Climate Reanalyzer
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
        TITLE_FONTSIZE = 16  # was 12
        axes[0].set_title('ERA-5 | Climate Reanalyzer', loc='left', y=TITLE_Y,
                            fontsize=TITLE_FONTSIZE, fontweight='bold')
        axes[1].set_title(f'{label}; 1979\u20132000 Weekly Mean', loc='right', y=TITLE_Y,
                            fontsize=TITLE_FONTSIZE, fontweight='bold')

        cbar = fig.colorbar(mesh, ax=axes, orientation='horizontal',
                             pad=0.08, fraction=0.12, shrink=0.5,
                             aspect=18, extend='neither')
        cbar.set_label('Wind Speed at 10 m (m/s)', fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_ticks(np.arange(0, 17, 2))
        cbar.set_ticklabels([str(t) for t in np.arange(0, 17, 2)])
        # NOTE: no divider-line loop here -- that's an SLP-specific style
        # (animate_SLP.py uses it, animate_wind.py doesn't), left over from
        # using the SLP script as a starting point. Confirmed by checking
        # animate_wind.py directly: no such loop in the original.

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