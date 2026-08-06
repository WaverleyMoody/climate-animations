"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: wind700_nicolosi.py
Description: Generates the 700mb wind speed climatology animation from
ERA5 reanalysis data (1979-2000), with mean wind vectors overlaid,
rendered as a double-hemisphere Nicolosi Globular projection.

Note: For the Plate Carrée, Robinson, and Foucaut projections, see the
other scripts in the wind_700mb scripts folder.
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

# ── Config ────────────────────────────────────────────────────────────────────
NC_PATH      = "/Volumes/CLIMATEDATA/wind700/climatology/climatology_wind700_48frame.nc"
OUT_MP4_PATH = "/Users/waverleymoody/Downloads/wind700_animation/wind700_climatology_nicolosi.mp4"
WEST_LON = -90.0
EAST_LON =  90.0
FPS = 4
DPI = 130

# ── Plot settings ─────────────────────────────────────────────────────────────
VMIN, VMAX = 0, 30

colors = [
    (0.00, '#ffffff'),  # white         (0 m/s)
    (0.20, '#d0e8f0'),  # light blue    (6 m/s)
    (0.35, '#469a46'),  # light green   (10.5 m/s)
    (0.44, '#ffff00'),  # yellow        (13 m/s)
    (0.50, '#ff6600'),  # orange        (15 m/s)
    (0.60, '#ff3300'),  # darker orange (18 m/s)
    (0.70, '#ff0000'),  # light red     (21 m/s)
    (0.80, '#cc0000'),  # red           (24 m/s)
    (0.90, '#990000'),  # dark red      (27 m/s)
    (1.00, '#660000'),  # darkest red   (30 m/s)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_wind700', colors, N=35)

# ── Quiver settings ───────────────────────────────────────────────────────────
# QUIVER_STEP subsamples the grid -- 28 gives a visually clean density for
# 700mb without overcrowding (tuned during 10m wind Nicolosi development).
# ARROW_SCALE is in projected meters per (m/s). 700mb winds are faster than
# 10m, so the same scale value produces longer arrows -- tune if needed.
# MAX_QUIVER_LAT excludes near-pole points where the Nicolosi projection is
# singular and the Jacobian blows up, causing starburst artifacts.
QUIVER_STEP    = 28
ARROW_SCALE    = 75000
ARROW_WIDTH    = 40000   # shaft width in projected meters (units='xy') -- thickened for visibility
MAX_QUIVER_LAT = 85


def build_hemisphere(lons, lats, central_lon):
    """Precompute everything that doesn't depend on frame data:
    the hemisphere mask, projected grid coords, boundary circle,
    clipped/projected coastlines, and the quiver Jacobian."""
    proj_str = f"+proj=nicol +lon_0={central_lon} +R=6371000"
    fwd = Transformer.from_crs("EPSG:4326", proj_str, always_xy=True)

    lon_diff = ((lons - central_lon + 180) % 360) - 180
    mask_1d  = np.abs(lon_diff) <= 90
    if mask_1d.sum() < 2:
        raise RuntimeError(f"Hemisphere mask kept almost no columns for lon_0={central_lon}")

    # Sort selected columns by angular offset -- fixes the wraparound bug.
    candidate_idx = np.where(mask_1d)[0]
    order         = np.argsort(lon_diff[candidate_idx])
    col_index     = candidate_idx[order]

    lons_hemi       = lons[col_index]
    lon2d, lat2d    = np.meshgrid(lons_hemi, lats)
    X, Y            = fwd.transform(lon2d, lat2d)

    edge_lats = np.linspace(-90, 90, 400)
    bx1, by1  = fwd.transform(np.full_like(edge_lats, central_lon + 90), edge_lats)
    bx2, by2  = fwd.transform(np.full_like(edge_lats, central_lon - 90), edge_lats)
    boundary_x = np.concatenate([bx1, bx2[::-1]])
    boundary_y = np.concatenate([by1, by2[::-1]])

    shp_path = shpreader.natural_earth(resolution="110m", category="physical", name="coastline")
    reader   = shpreader.Reader(shp_path)
    hemi_bbox = box(central_lon - 90, -90, central_lon + 90, 90)

    coastlines = []
    for record in reader.geometries():
        clipped = record.intersection(hemi_bbox)
        if clipped.is_empty:
            continue
        proj_geom = shp_transform(lambda x, y: fwd.transform(x, y), clipped)
        coastlines.append(proj_geom)

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

    # ── Quiver geometry ───────────────────────────────────────────────────────
    # Subsample grid and precompute the local rotation Jacobian at each point.
    # Exclude near-pole rows: Nicolosi is singular at the poles (all longitudes
    # collapse to one point), so the Jacobian blows up there and produces
    # starburst artifacts. MAX_QUIVER_LAT caps the latitude range.
    all_quiver_rows = np.arange(0, len(lats), QUIVER_STEP)
    quiver_rows     = all_quiver_rows[np.abs(lats[all_quiver_rows]) <= MAX_QUIVER_LAT]
    quiver_cols     = col_index[::QUIVER_STEP]

    qlon2d, qlat2d = np.meshgrid(lons[quiver_cols], lats[quiver_rows])
    quiver_X, quiver_Y = fwd.transform(qlon2d, qlat2d)

    R   = 6371000.0
    eps = 0.05
    qx0, qy0 = fwd.transform(qlon2d,       qlat2d)
    qxe, qye = fwd.transform(qlon2d + eps, qlat2d)
    qxn, qyn = fwd.transform(qlon2d,       qlat2d + eps)
    east_dist  = eps * np.pi / 180 * R * np.cos(np.radians(qlat2d))
    north_dist = eps * np.pi / 180 * R
    east_dist  = np.where(np.abs(east_dist) < 1e-6, 1e-6, east_dist)
    quiver_dxde = (qxe - qx0) / east_dist
    quiver_dyde = (qye - qy0) / east_dist
    quiver_dxdn = (qxn - qx0) / north_dist
    quiver_dydn = (qyn - qy0) / north_dist

    # Build a Shapely polygon from the boundary for arrow clipping.
    from shapely.geometry import Polygon
    boundary_poly = Polygon(zip(boundary_x, boundary_y))

    return {
        "col_index":     col_index,
        "X":             X,
        "Y":             Y,
        "boundary_x":    boundary_x,
        "boundary_y":    boundary_y,
        "boundary_poly": boundary_poly,
        "coastlines":    coastlines,
        "borders":       borders,
        "quiver_rows":   quiver_rows,
        "quiver_cols":   quiver_cols,
        "quiver_X":      quiver_X,
        "quiver_Y":      quiver_Y,
        "quiver_dxde":   quiver_dxde,
        "quiver_dyde":   quiver_dyde,
        "quiver_dxdn":   quiver_dxdn,
        "quiver_dydn":   quiver_dydn,
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
    print("Variables:", list(ds.data_vars))

    speed_da = ds["wind_speed_700mb"]
    u_da     = ds["u700"]
    v_da     = ds["v700"]

    lats = speed_da["latitude"].values
    lons = speed_da["longitude"].values
    n_frames = speed_da.sizes["frame"]
    print(f"{n_frames} frames, grid {lats.size} x {lons.size}")

    # No frame_label in this dataset -- generate week labels from frame index.
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    week_starts = [1, 8, 15, 22]
    frame_labels = [f'{month_names[m]} {d}'
                    for m in range(12) for d in week_starts]

    actual_min = float(speed_da.min())
    actual_max = float(speed_da.max())
    print(f"Actual wind speed range: {actual_min:.1f} to {actual_max:.1f} m/s "
          f"(color scale fixed at {VMIN} to {VMAX})")
    if actual_min < VMIN or actual_max > VMAX:
        print("  NOTE: some values fall outside the fixed color scale and will clip.")

    print("Precomputing Western Hemisphere geometry...")
    west = build_hemisphere(lons, lats, WEST_LON)
    print("Precomputing Eastern Hemisphere geometry...")
    east = build_hemisphere(lons, lats, EAST_LON)

    def render_frame(i):
        # squeeze() drops the pressure_level dimension (size 1).
        speed_frame = speed_da.isel(frame=i).squeeze().values
        u_frame     = u_da.isel(frame=i).squeeze().values
        v_frame     = v_da.isel(frame=i).squeeze().values
        label       = frame_labels[i]

        fig, axes = plt.subplots(1, 2, figsize=(13, 7), constrained_layout=True)

        for ax, hemi in [(axes[0], west), (axes[1], east)]:
            data_hemi = speed_frame[:, hemi["col_index"]]
            mesh = ax.pcolormesh(hemi["X"], hemi["Y"], data_hemi,
                                 shading="auto", cmap=CMAP,
                                 vmin=VMIN, vmax=VMAX)
            plot_coastlines(ax, hemi["coastlines"], color="black", lw=1.0)
            plot_coastlines(ax, hemi["borders"],    color="black", lw=0.8)
            ax.plot(hemi["boundary_x"], hemi["boundary_y"], color="black", linewidth=1)

            # Quiver overlay: pull u/v at the precomputed subsampled grid
            # points, then rotate them into map-local x/y using the
            # Jacobian computed once in build_hemisphere.
            u_q = u_frame[np.ix_(hemi["quiver_rows"], hemi["quiver_cols"])]
            v_q = v_frame[np.ix_(hemi["quiver_rows"], hemi["quiver_cols"])]

            # Normalize: preserve wind direction from the Jacobian but scale
            # arrow length by true wind speed magnitude. Prevents projection
            # distortion from producing misleading arrow lengths.
            dx_raw     = hemi["quiver_dxde"] * u_q + hemi["quiver_dxdn"] * v_q
            dy_raw     = hemi["quiver_dyde"] * u_q + hemi["quiver_dydn"] * v_q
            wind_speed = np.sqrt(u_q**2 + v_q**2)
            proj_norm  = np.sqrt(dx_raw**2 + dy_raw**2)
            valid      = (proj_norm > 1e-10) & (wind_speed > 0.05)
            u_map = np.where(valid, dx_raw / proj_norm * wind_speed, 0.0)
            v_map = np.where(valid, dy_raw / proj_norm * wind_speed, 0.0)

            # Clip arrows whose tip extends outside the hemisphere boundary.
            # Near-edge grid points have valid start positions inside the
            # boundary, but their arrow tips can land just outside the circle,
            # producing arrows that visually poke through the outline.
            from shapely.geometry import Point
            qX_flat = hemi["quiver_X"].ravel()
            qY_flat = hemi["quiver_Y"].ravel()
            u_flat  = u_map.ravel()
            v_flat  = v_map.ravel()
            tip_x   = qX_flat + u_flat * ARROW_SCALE
            tip_y   = qY_flat + v_flat * ARROW_SCALE
            poly    = hemi["boundary_poly"]
            in_bounds = np.array([
                poly.contains(Point(tx, ty))
                for tx, ty in zip(tip_x, tip_y)
            ])
            in_bounds = in_bounds.reshape(hemi["quiver_X"].shape)
            u_map = np.where(in_bounds, u_map, 0.0)
            v_map = np.where(in_bounds, v_map, 0.0)

            ax.quiver(hemi["quiver_X"], hemi["quiver_Y"], u_map, v_map,
                      color='black', scale=1 / ARROW_SCALE, scale_units='xy',
                      units='xy', width=ARROW_WIDTH,
                      headwidth=3, headlength=4, headaxislength=3.5)

            ax.set_aspect("equal")
            ax.axis("off")

        TITLE_Y        = 1.02
        TITLE_FONTSIZE = 16
        axes[0].set_title('ERA-5 | Climate Reanalyzer', loc='left', y=TITLE_Y,
                           fontsize=TITLE_FONTSIZE, fontweight='bold')
        axes[1].set_title(f'{label}; 1979\u20132000 Weekly Mean', loc='right', y=TITLE_Y,
                           fontsize=TITLE_FONTSIZE, fontweight='bold')

        cbar = fig.colorbar(mesh, ax=axes, orientation='horizontal',
                            pad=0.08, fraction=0.12, shrink=0.5,
                            aspect=18, extend='neither')
        cbar.set_label('Wind Speed at 700 mb (m/s)', fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_ticks(np.arange(0, 31, 5))
        cbar.set_ticklabels([str(t) for t in np.arange(0, 31, 5)])

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=DPI)
        plt.close(fig)
        buf.seek(0)
        img = imageio.imread(buf)

        # Drop alpha channel and force even dimensions (required by libx264).
        img = img[:, :, :3]
        h, w = img.shape[:2]
        img = img[: h - (h % 2), : w - (w % 2)]
        return np.ascontiguousarray(img)

    print("Rendering frame 0...")
    first_frame = render_frame(0)
    height, width = first_frame.shape[:2]

    ffmpeg_writer = imageio_ffmpeg.write_frames(
        OUT_MP4_PATH, size=(width, height), fps=FPS,
        codec="libx264", quality=8, macro_block_size=1,
    )
    ffmpeg_writer.send(None)
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