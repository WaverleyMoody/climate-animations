"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: south_stereo.py
Description: Generates the Southern Hemisphere seasonal-cycle stereographic animation from NASA 
Blue Marble Next Generation monthly imagery (2004) with an ERA5 sea_ice_cover overlay 
(24 biweekly frames, 1st and 15th of each month), rendered in the South Polar Stereographic projection.
Note: For the Northern Hemisphere, see north_stereo.py in the stereographic scripts folder.
"""

import datetime as dt
from pathlib import Path

import cartopy.crs as ccrs
import imageio.v2 as imageio
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from PIL import Image

# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------

BLUE_MARBLE_DIR = Path("/Volumes/CLIMATEDATA/blue_marble/2004")
ICE_FILE = Path("/Volumes/CLIMATEDATA/era5_sea_ice/2004/era5_sea_ice_cover_2004_biweekly.nc")

FRAMES_DIR = Path("/Volumes/CLIMATEDATA/stereographic/south/frames")
VIDEO_OUT = Path("/Volumes/CLIMATEDATA/stereographic/south/south_seasonal_cycle.mp4")

# Target frame dates: 1st and 15th of each month, 2004
TARGET_DATES = [
    dt.date(2004, m, d) for m in range(1, 13) for d in (1, 15)
]

ICE_THRESHOLD = 0.15   # minimum concentration to render as visible ice edge
SOUTH_MAX_LAT = 0      # northern boundary of the hemisphere view (equator)
SOUTH_CENTRAL_LON = 300  # meridian pointing to the top of the map -- SouthPolarStereo's
                          # rotation direction differs from NorthPolarStereo, so this is
                          # a starting guess carried over from the north script; retune
                          # by eye rather than assuming it gives the same orientation
FPS = 6

BASE_EXTENT = (-180, 180, -90, 90)  # Blue Marble global equirectangular extent


# ----------------------------------------------------------------------------
# Blue Marble loading and biweekly blending
# ----------------------------------------------------------------------------

def month_file(month: int) -> Path:
    return BLUE_MARBLE_DIR / f"world.2004{month:02d}.3x5400x2700.jpg"


def load_month_image(month: int) -> np.ndarray:
    """Load a monthly Blue Marble composite as a float32 array in [0, 1]."""
    img = Image.open(month_file(month)).convert("RGB")
    return np.asarray(img, dtype=np.float32) / 255.0


def blend_base_image(target_date: dt.date) -> np.ndarray:
    """
    Blue Marble is monthly, not biweekly. Blend the two nearest monthly
    composites (each anchored to the 15th of its month) to build a base
    image for the target biweekly date. Dates outside the first/last
    anchor fall back to the nearest single month with no blending.
    """
    centers = {m: dt.date(2004, m, 15) for m in range(1, 13)}

    if target_date <= centers[1]:
        return load_month_image(1)
    if target_date >= centers[12]:
        return load_month_image(12)

    prev_month = max(m for m, c in centers.items() if c <= target_date)
    next_month = prev_month + 1
    prev_center, next_center = centers[prev_month], centers[next_month]

    if target_date == prev_center:
        return load_month_image(prev_month)

    span = (next_center - prev_center).days
    weight = (target_date - prev_center).days / span

    prev_img = load_month_image(prev_month)
    next_img = load_month_image(next_month)
    return (1 - weight) * prev_img + weight * next_img


# ----------------------------------------------------------------------------
# ERA5 sea ice loading
# ----------------------------------------------------------------------------

def get_time_dim(ds: xr.Dataset) -> str:
    return "valid_time" if "valid_time" in ds.dims else "time"


def load_ice_dataset() -> xr.Dataset:
    ds = xr.open_dataset(ICE_FILE)
    # ERA5 uses 0-360 longitude; normalize to -180/180 and re-sort before reprojecting
    ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
    ds = ds.sortby("longitude")
    return ds


def get_ice_frame(ds: xr.Dataset, target_date: dt.date) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (lon, lat, concentration) for the nearest available date."""
    time_dim = get_time_dim(ds)
    target = np.datetime64(target_date.isoformat())
    frame = ds.sel({time_dim: target}, method="nearest")

    var_name = list(ds.data_vars)[0]  # robust to siconc/ci naming differences
    sic = frame[var_name].values.astype(np.float32)
    lon = frame["longitude"].values
    lat = frame["latitude"].values
    return lon, lat, sic


def build_ice_rgba(sic: np.ndarray, threshold: float = ICE_THRESHOLD) -> np.ndarray:
    """
    Build an RGBA overlay: white where ice concentration exceeds the
    threshold, transparent elsewhere. ERA5 sea_ice_cover is NaN over land,
    which doubles as a free ocean mask -- land pixels get alpha=0.
    """
    ny, nx = sic.shape
    rgba = np.zeros((ny, nx, 4), dtype=np.float32)
    rgba[..., 0:3] = 1.0  # white

    alpha = np.where(sic >= threshold, np.clip(sic, 0, 1), 0.0)
    alpha = np.nan_to_num(alpha, nan=0.0)  # NaN (land) -> fully transparent
    rgba[..., 3] = alpha
    return rgba


# ----------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------

def polar_boundary_circle():
    """Standard Cartopy recipe for a round clip boundary on a polar stereographic axes."""
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)


def render_frame(target_date: dt.date, base_img: np.ndarray,
                  lon: np.ndarray, lat: np.ndarray, ice_rgba: np.ndarray,
                  output_path: Path) -> None:
    fig = plt.figure(figsize=(8, 8), dpi=150, facecolor="white")
    proj = ccrs.SouthPolarStereo(central_longitude=SOUTH_CENTRAL_LON)
    # Fixed, symmetric rect (left, bottom, width, height) in figure fraction --
    # equal left/right AND top/bottom margins keep the globe centered both
    # horizontally and vertically, and using an explicit rect (rather than
    # default placement + a tight-crop save) guarantees every frame renders
    # at the same size and position.
    ax = fig.add_axes([0.05, 0.09, 0.90, 0.82], projection=proj)
    ax.set_facecolor("white")

    ax.set_extent([-180, 180, -90, SOUTH_MAX_LAT], ccrs.PlateCarree())
    ax.set_boundary(polar_boundary_circle(), transform=ax.transAxes)

    # Base layer: reprojected Blue Marble imagery
    ax.imshow(
        base_img, origin="upper", extent=BASE_EXTENT,
        transform=ccrs.PlateCarree(), interpolation="bilinear",
    )

    # Overlay: reprojected sea ice concentration
    ice_extent = (lon.min(), lon.max(), lat.min(), lat.max())
    ice_origin = "upper" if lat[0] > lat[-1] else "lower"
    ax.imshow(
        ice_rgba, origin=ice_origin, extent=ice_extent,
        transform=ccrs.PlateCarree(), interpolation="nearest",
    )

    fig.suptitle("Southern Hemisphere Seasonal Cycle", fontsize=16, y=0.97, x=0.5, ha="center")

    fig.text(
        0.98, 0.02, target_date.strftime("%b %d, %Y"),
        ha="right", va="bottom", fontsize=12, color="black",
        path_effects=None,
    )

    fig.savefig(output_path, facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    FRAMES_DIR.mkdir(parents=True, exist_ok=True)
    ice_ds = load_ice_dataset()

    frame_paths = []
    for target_date in TARGET_DATES:
        print(f"Rendering frame: {target_date.isoformat()}")
        base_img = blend_base_image(target_date)
        lon, lat, sic = get_ice_frame(ice_ds, target_date)
        ice_rgba = build_ice_rgba(sic)

        out_path = FRAMES_DIR / f"south_{target_date.isoformat()}.png"
        render_frame(target_date, base_img, lon, lat, ice_rgba, out_path)
        frame_paths.append(out_path)

    ice_ds.close()

    print(f"Assembling {len(frame_paths)} frames into {VIDEO_OUT}...")
    VIDEO_OUT.parent.mkdir(parents=True, exist_ok=True)
    with imageio.get_writer(VIDEO_OUT, fps=FPS) as writer:
        for frame_path in frame_paths:
            writer.append_data(imageio.imread(frame_path))

    print("Done.")


if __name__ == "__main__":
    main()