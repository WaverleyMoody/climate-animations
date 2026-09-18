"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: animate_gpp_nicolosi.py
Description: Generates the Gross Primary Productivity climatology
    animation from MOD17A2HGF (2000-2009), rendered as a double-
    hemisphere Nicolosi Globular projection.
Note: For the Plate Carrée, Robinson, and Foucaut projections, see
    animate_gpp_platecarree.py, animate_gpp_robinson.py, and
    animate_gpp_foucaut.py in the gpp scripts folder. Data excludes
    2000-01-01 through 2000-02-10 (early Terra commissioning gap) -
    the animation effectively starts 2000-02-18, 454 frames total.

PROJ's `nicol` operation has no inverse transform, so Cartopy's
GeoAxes (which needs one to draw gridlines/features and to auto-detect
wrapped cells) can't be used here - unlike the Plate Carrée, Robinson,
and Foucaut siblings. This instead forward-projects the data and
Natural Earth land/coastline geometries by hand with `pyproj`, on a
pair of plain (non-cartopy) Matplotlib axes - one per hemisphere, the
same manual pipeline already established in this project's
2m_temp_nicolosi.py et al.
"""

import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FuncAnimation
from matplotlib.cm import ScalarMappable
from matplotlib.colors import PowerNorm
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from pyproj import Transformer
from shapely.geometry import box
from shapely.ops import transform as shp_transform
from tqdm import tqdm

# ---- Config -------------------------------------------------------------
DATA_PATH = "/Volumes/CLIMATEDATA/gpp_2000_2009.nc"
OUTPUT_PATH = "/Volumes/CLIMATEDATA/gpp_2000_2009_nicolosi.mp4"
WEST_LON = -90.0   # hemisphere centers, matching the sibling Nicolosi scripts
EAST_LON = 90.0
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
# for its tick labels and axis label, which are reserved via the hemisphere
# rects below rather than left to matplotlib's default margins.
COLORBAR_RECT = [0.28, 0.085, 0.44, 0.035]
COLORBAR_LABEL = "GPP (kg C/m$^2$/8-day)"
COLORBAR_TICKS = [0, 0.03, 0.06, 0.09, 0.12]
COLORBAR_TEXT_COLOR = "black"

# Hemisphere axes placement, figure-fraction [left, bottom, width, height].
# Same bottom/height as the single-map siblings' MAP_RECT, split into two
# side-by-side squares with a small gutter between them.
WEST_RECT = [0.03, 0.13, 0.44, 0.80]
EAST_RECT = [0.53, 0.13, 0.44, 0.80]
TITLE_Y = 1.02  # axis("off") leaves no spine to anchor "auto" title placement to


def load_data():
    ds = xr.open_dataset(DATA_PATH, chunks={"time": 1})
    return ds["gpp"]


def _polygons_to_patch(geom, **kwargs):
    """
    Converts a Shapely Polygon/MultiPolygon into a single Matplotlib
    PathPatch, holes (interior rings) included via Path's MOVETO/LINETO
    codes - the standard way to carry Shapely geometry into Matplotlib
    without a third-party bridge like descartes.
    """
    polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    vertices, codes = [], []
    for poly in polys:
        for ring in (poly.exterior, *poly.interiors):
            coords = np.asarray(ring.coords)
            ring_codes = np.full(len(coords), MplPath.LINETO, dtype=MplPath.code_type)
            ring_codes[0] = MplPath.MOVETO
            vertices.append(coords)
            codes.append(ring_codes)
    if not vertices:
        return None
    return PathPatch(MplPath(np.concatenate(vertices), np.concatenate(codes)), **kwargs)


def build_hemisphere(lons, lats, central_lon):
    """
    Precomputes everything that doesn't depend on frame data: the
    forward-projected grid coordinates and column ordering, the boundary
    circle, and the clipped/reprojected ocean fill, land fill, and
    coastline for this hemisphere.
    """
    proj_str = f"+proj=nicol +lon_0={central_lon} +R=6371000"
    fwd = Transformer.from_crs("EPSG:4326", proj_str, always_xy=True)

    lon_diff = ((lons - central_lon + 180) % 360) - 180  # angular offset from center
    mask_1d = np.abs(lon_diff) <= 90
    if mask_1d.sum() < 2:
        raise RuntimeError(f"Hemisphere mask kept almost no columns for lon_0={central_lon}")

    # Sort the selected columns by angular offset from the hemisphere
    # center rather than trusting the original array's storage order -
    # lon_diff has no discontinuity within a single +/-90deg window (its
    # only jump is at +/-180 from center, always excluded by the mask
    # above), so this ordering is safe for any input longitude convention.
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
    boundary_path = MplPath(np.column_stack([boundary_x, boundary_y]))

    hemi_bbox = box(central_lon - 90, -90, central_lon + 90, 90)

    def clipped_reprojected(category, name):
        shp_path = shpreader.natural_earth(resolution="110m", category=category, name=name)
        reader = shpreader.Reader(shp_path)
        geoms = []
        for record in reader.geometries():
            clipped = record.intersection(hemi_bbox)
            if clipped.is_empty:
                continue
            geoms.append(shp_transform(lambda x, y: fwd.transform(x, y), clipped))
        return geoms

    land_polys = clipped_reprojected("physical", "land")
    coastlines = clipped_reprojected("physical", "coastline")

    ocean_patch = PathPatch(boundary_path, facecolor=OCEAN_COLOR, edgecolor="none", zorder=0)
    land_patches = [
        p for geom in land_polys
        if (p := _polygons_to_patch(geom, facecolor=LAND_BASE_COLOR, edgecolor="none", zorder=1))
    ]
    boundary_patch = PathPatch(
        boundary_path, facecolor="none", edgecolor="gray", linewidth=0.6, zorder=3,
    )

    return {
        "col_index": col_index,
        "X": X,
        "Y": Y,
        "ocean_patch": ocean_patch,
        "land_patches": land_patches,
        "coastlines": coastlines,
        "boundary_patch": boundary_patch,
    }


def setup_hemisphere_axes(ax, hemi):
    """
    Draws everything static for this hemisphere - ocean fill, land fill,
    coastline, boundary circle - once, before the animation starts. None
    of it depends on frame data, so it's never touched again per frame;
    only the GPP pcolormesh gets added/removed each frame (see
    render_frame), which is both simpler and far cheaper than clearing
    and rebuilding the whole axes 454 times.
    """
    ax.add_patch(hemi["ocean_patch"])
    for patch in hemi["land_patches"]:
        ax.add_patch(patch)
    for geom in hemi["coastlines"]:
        lines = geom.geoms if geom.geom_type == "MultiLineString" else [geom]
        for line in lines:
            xs, ys = line.xy
            ax.plot(xs, ys, color="gray", linewidth=0.3, zorder=3)
    ax.add_patch(hemi["boundary_patch"])

    ax.set_xlim(hemi["X"].min(), hemi["X"].max())
    ax.set_ylim(hemi["Y"].min(), hemi["Y"].max())
    ax.set_aspect("equal")
    ax.axis("off")


def render_frame(ax, hemi, da_frame, prev_mesh):
    if prev_mesh is not None:
        prev_mesh.remove()

    data_hemi = da_frame.values[:, hemi["col_index"]]
    # NaN-transparent overlay: ocean/no-data pixels show the land/ocean
    # fill underneath rather than a solid color, matching the original
    # clip's look of GPP-colored land against a plain navy ocean.
    mesh = ax.pcolormesh(
        hemi["X"], hemi["Y"], np.ma.masked_invalid(data_hemi),
        shading="auto", cmap=CMAP, norm=PowerNorm(gamma=GAMMA, vmin=VMIN, vmax=VMAX),
        zorder=2,
    )
    return mesh


def _add_gpp_colorbar(fig):
    """
    Unlike the single-map siblings, the Nicolosi axes are never cleared
    per frame (see setup_hemisphere_axes), and VMIN/VMAX/GAMMA/CMAP never
    change between frames - so the colorbar's mapping is identical every
    frame and only needs to be built once, from a standalone
    ScalarMappable rather than any particular frame's mesh.
    """
    sm = ScalarMappable(norm=PowerNorm(gamma=GAMMA, vmin=VMIN, vmax=VMAX), cmap=CMAP)
    cax = fig.add_axes(COLORBAR_RECT)
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")

    cbar.set_ticks(COLORBAR_TICKS)
    cbar.set_label(COLORBAR_LABEL, color=COLORBAR_TEXT_COLOR, fontsize=11, labelpad=6)
    cbar.ax.tick_params(labelsize=9, colors=COLORBAR_TEXT_COLOR)
    cbar.outline.set_edgecolor(COLORBAR_TEXT_COLOR)
    cbar.outline.set_linewidth(0.6)


def main():
    da = load_data()
    n_frames = da.sizes["time"]
    print(f"Rendering {n_frames} frames at {FPS} fps...")

    lons, lats = da.x.values, da.y.values
    print("Precomputing Western Hemisphere geometry...")
    west = build_hemisphere(lons, lats, WEST_LON)
    print("Precomputing Eastern Hemisphere geometry...")
    east = build_hemisphere(lons, lats, EAST_LON)

    fig = plt.figure(figsize=(12, 6))
    ax_west = fig.add_axes(WEST_RECT)
    ax_east = fig.add_axes(EAST_RECT)
    setup_hemisphere_axes(ax_west, west)
    setup_hemisphere_axes(ax_east, east)
    _add_gpp_colorbar(fig)

    mesh_state = {"west": None, "east": None}
    progress = tqdm(total=n_frames)

    def update(i):
        frame = da.isel(time=i).load()
        mesh_state["west"] = render_frame(ax_west, west, frame, mesh_state["west"])
        mesh_state["east"] = render_frame(ax_east, east, frame, mesh_state["east"])

        date_str = pd.Timestamp(frame.time.values).strftime("%b %-d, %Y")
        ax_west.set_title("Gross Primary Productivity", loc="left", y=TITLE_Y, fontsize=12)
        ax_east.set_title(date_str, loc="right", y=TITLE_Y, fontsize=12)

        progress.update(1)
        return [mesh_state["west"], mesh_state["east"]]

    # cache_frame_data=False - OOM discipline on the 8GB Air; don't let
    # FuncAnimation hold references to every frame's data simultaneously
    anim = FuncAnimation(
        fig, update, frames=n_frames, cache_frame_data=False,
    )

    anim.save(OUTPUT_PATH, fps=FPS, writer="ffmpeg", dpi=150)
    progress.close()
    print(f"\nSaved animation to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()