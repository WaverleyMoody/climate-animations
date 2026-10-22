"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: arctic_smoke_orthographic.py
Description: Generates the Arctic boreal forest wildfire smoke animation from
    NASA MERRA-2 hourly total aerosol optical depth TOTEXTTAU (2019-06-01 to
    2019-08-31), rendered in the Orthographic projection.
Note: For the data download, see download_merra2_aod_2019_jja.py in the
    arctic_smoke scripts folder.

Run:
    conda activate climate
    # 1) Set PREVIEW_ONLY = True, check the PNG
    # 2) Set PREVIEW_ONLY = False for the full render
    caffeinate -i python arctic_smoke_orthographic.py
"""

import shutil
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import dask
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
VARIABLE = "TOTEXTTAU"
DATA_PATH = Path("/Volumes/CLIMATEDATA/smoke/merra2_2019_jja/MERRA2_TOTEXTTAU_hourly_2019-06-08.nc")

LOCAL_DIR = Path.home() / "CLIMATE ANIMATIONS" / "arctic_smoke" / "_render"   # local APFS
FINAL_DIR = Path("/Volumes/CLIMATEDATA/smoke/animations")                     # exFAT
OUT_NAME = "arctic_smoke_jja2019_orthographic.mp4"

PREVIEW_ONLY = False                       # True: save one PNG and stop
PREVIEW_TIME = "2019-06-13T03:30"         # matches the reference frame (04z)

# Globe view, measured from the reference frame
CENTRAL_LON = -115.0
CENTRAL_LAT = 70.0
LAT_CIRCLES = [63, 80]                    # dashed parallels in the reference

# Resolution of the precomputed projection grid (pixels across the globe)
GRID_N = 1400

# Colorbar levels from the original (log-spaced, extended both ends)
LEVELS = [0.10, 0.15, 0.21, 0.31, 0.45, 0.66, 0.97, 1.40, 2.10, 3.00]

BG_COLOR = "#fdf5e8"                      # cream globe background
CMAP_ANCHORS = ["#fbe9cf", "#f6cf9c", "#eeac72", "#e2864f",
                "#cd5f34", "#ad3c21", "#842312", "#5a0e06"]
OVER_COLOR = "#3d0703"

# Line styles: light coastlines worldwide, dark outlines for US/Canada + states
LIGHT_LINE = "#9a9a9a"
DARK_LINE = "#3a3a3a"
DARK_COUNTRIES = ["USA", "CAN"]

SHOW_US_CITIES = True
CITY_MIN_POP = 50_000

TITLE = "Smoke from Arctic Boreal Forest Wildfires, June\u2013August 2019"
SUBTITLE = "NASA MERRA-2 hourly aerosol optical depth (550 nm)"
CBAR_LABEL = "Dust optical depth (550nm)"

FPS = 35                                  # 2208 frames -> 63 s, same as original
FIGSIZE = (12, 10.4)                      # 1800 x 1560 px at DPI 150
DPI = 150
BITRATE = 12000

dask.config.set(scheduler="synchronous")
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_cmap():
    smooth = LinearSegmentedColormap.from_list("aod", CMAP_ANCHORS)
    n_bins = len(LEVELS) - 1
    cmap = ListedColormap(smooth(np.linspace(0, 1, n_bins)), name="aod_binned")
    cmap.set_under(BG_COLOR)
    cmap.set_over(OVER_COLOR)
    cmap.set_bad(BG_COLOR)
    return cmap, BoundaryNorm(LEVELS, n_bins)


def _records(resolution, category, name):
    path = shpreader.natural_earth(resolution=resolution, category=category, name=name)
    for rec in shpreader.Reader(path).records():
        attrs = {k.upper(): v for k, v in rec.attributes.items()}
        yield attrs, rec.geometry


def _us_city_points():
    try:
        lons, lats = [], []
        for a, geom in _records("10m", "cultural", "populated_places"):
            if a.get("ADM0_A3") == "USA" and (a.get("POP_MAX") or 0) >= CITY_MIN_POP:
                lons.append(geom.x)
                lats.append(geom.y)
        return np.array(lons), np.array(lats)
    except Exception as exc:
        print(f"City layer skipped ({exc})")
        return np.array([]), np.array([])


def _stamp(t):
    return pd.Timestamp(t).strftime("%Y-%m-%d %H:%Mz")


class OrthoResampler:
    """Precompute bilinear lookups from the lat-lon grid onto the globe image.

    Reprojecting with imshow(transform=PlateCarree) every frame is slow; this
    does the projection math once, so each frame is just an array gather.
    """

    def __init__(self, proj, lon, lat, n):
        r = proj.x_limits[1]
        self.extent = [-r, r, -r, r]
        self.n = n

        centers = -r + (np.arange(n) + 0.5) * (2 * r / n)
        xx, yy = np.meshgrid(centers, centers)
        inside = (xx ** 2 + yy ** 2) < r ** 2

        pts = ccrs.PlateCarree().transform_points(proj, xx[inside], yy[inside])
        plon, plat = pts[:, 0], pts[:, 1]

        nlon, nlat = lon.size, lat.size
        dlon, dlat = lon[1] - lon[0], lat[1] - lat[0]

        fi = np.mod((plon - lon[0]) / dlon, nlon)
        i0 = np.floor(fi).astype(np.int64)
        wx = (fi - i0).astype(np.float32)
        i0 = i0 % nlon
        i1 = (i0 + 1) % nlon

        fj = (plat - lat[0]) / dlat
        valid = np.isfinite(fj) & (fj >= 0) & (fj <= nlat - 1)
        j0 = np.clip(np.floor(np.nan_to_num(fj)), 0, nlat - 2).astype(np.int64)
        wy = np.clip(np.nan_to_num(fj) - j0, 0, 1).astype(np.float32)
        j1 = j0 + 1

        flat_inside = np.flatnonzero(inside.ravel())
        self.target = flat_inside[valid]
        self.i00 = (j0 * nlon + i0)[valid]
        self.i01 = (j0 * nlon + i1)[valid]
        self.i10 = (j1 * nlon + i0)[valid]
        self.i11 = (j1 * nlon + i1)[valid]
        self.wx = wx[valid]
        self.wy = wy[valid]

    def __call__(self, field):
        f = np.asarray(field, dtype=np.float32).ravel()
        top = (1 - self.wx) * f[self.i00] + self.wx * f[self.i01]
        bot = (1 - self.wx) * f[self.i10] + self.wx * f[self.i11]
        out = np.full(self.n * self.n, np.nan, dtype=np.float32)
        out[self.target] = (1 - self.wy) * top + self.wy * bot
        return out.reshape(self.n, self.n)


def _build_figure(resampler, proj):
    cmap, norm = _build_cmap()
    pc = ccrs.PlateCarree()

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI, facecolor="white")
    ax = fig.add_axes([0.01, 0.03, 0.85, 0.86], projection=proj)
    ax.set_global()
    ax.set_facecolor(BG_COLOR)
    ax.spines["geo"].set_visible(False)

    im = ax.imshow(np.full((resampler.n, resampler.n), np.nan), origin="lower",
                   extent=resampler.extent, transform=proj, cmap=cmap, norm=norm,
                   interpolation="bilinear", zorder=1)

    # Light layers worldwide
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.45,
                   edgecolor=LIGHT_LINE, zorder=2)
    ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="none",
                   linewidth=0.35, edgecolor=LIGHT_LINE, zorder=2)

    # Dark outlines for the US and Canada, plus US states
    countries = [g for a, g in _records("50m", "cultural", "admin_0_countries")
                 if a.get("ADM0_A3") in DARK_COUNTRIES]
    ax.add_geometries(countries, crs=pc, facecolor="none", edgecolor=DARK_LINE,
                      linewidth=0.6, zorder=3)
    states = [g for a, g in _records("50m", "cultural", "admin_1_states_provinces_lakes")
              if a.get("ADM0_A3") == "USA"]
    ax.add_geometries(states, crs=pc, facecolor="none", edgecolor=DARK_LINE,
                      linewidth=0.45, zorder=3)

    # Dashed latitude circles
    circle_lons = np.linspace(-180, 180, 721)
    for lat in LAT_CIRCLES:
        ax.plot(circle_lons, np.full_like(circle_lons, lat), transform=pc,
                color="#777777", linewidth=0.7, linestyle=(0, (4, 4)), zorder=3)

    if SHOW_US_CITIES:
        clon, clat = _us_city_points()
        if clon.size:
            ax.scatter(clon, clat, s=1.5, c="black", linewidths=0,
                       transform=pc, zorder=4)

    fig.text(0.435, 0.955, TITLE, ha="center", va="bottom",
             fontsize=17, fontweight="bold")
    fig.text(0.435, 0.925, SUBTITLE, ha="center", va="bottom",
             fontsize=12, color="#444444")

    cax = fig.add_axes([0.895, 0.27, 0.018, 0.46])
    cb = fig.colorbar(im, cax=cax, extend="both", ticks=LEVELS)
    cb.ax.set_yticklabels([f"{v:.2f}" for v in LEVELS], fontsize=13)
    cb.set_label(CBAR_LABEL, fontsize = 14, labelpad = 8)
    cb.ax.yaxis.set_label_position("left")

    stamp = fig.text(0.93, 0.03, "", ha="right", va="bottom", fontsize=13)
    return fig, im, stamp


def _verify_mp4(path, expected_frames):
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"{path.name} is missing or empty")
    n_frames, secs = imageio_ffmpeg.count_frames_and_secs(str(path))
    if n_frames != expected_frames:
        raise ValueError(f"{path.name}: {n_frames} frames, expected {expected_frames}")
    print(f"Verified {path.name}: {n_frames} frames, {secs:.1f} s")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if not DATA_PATH.exists():
        sys.exit(f"{DATA_PATH} not found. Is /Volumes/CLIMATEDATA mounted?")
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)

    ds = xr.open_dataset(DATA_PATH, chunks={"time": 1})
    da = ds[VARIABLE].transpose("time", "lat", "lon")
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.isel(lat=slice(None, None, -1))
    if not np.all(np.diff(da.lon.values) > 0):
        sys.exit("Longitude axis is not ascending; check the input file")

    times = da.time.values
    n = times.size
    print(f"Loaded {VARIABLE}: {n} hourly frames, {da.sizes['lat']} x {da.sizes['lon']}")

    proj = ccrs.Orthographic(central_longitude=CENTRAL_LON, central_latitude=CENTRAL_LAT)
    print("Precomputing orthographic lookup grid...")
    resampler = OrthoResampler(proj, da.lon.values, da.lat.values, GRID_N)

    fig, im, stamp = _build_figure(resampler, proj)

    if PREVIEW_ONLY:
        t = da.sel(time=PREVIEW_TIME, method="nearest")
        im.set_data(resampler(t.values))
        stamp.set_text(_stamp(t.time.values))
        out = LOCAL_DIR / "preview_frame.png"
        fig.savefig(out, dpi=DPI)
        print(f"Preview saved: {out}")
        return

    local_mp4 = LOCAL_DIR / OUT_NAME
    writer = FFMpegWriter(fps=FPS, codec="libx264", bitrate=BITRATE,
                          extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, str(local_mp4), dpi=DPI):
        for i in range(n):
            im.set_data(resampler(da.isel(time=i).values))
            stamp.set_text(_stamp(times[i]))
            writer.grab_frame()
            if i % 24 == 0:
                print(f"  frame {i + 1}/{n}  {_stamp(times[i])}")
    plt.close(fig)
    ds.close()

    _verify_mp4(local_mp4, n)

    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    final_mp4 = FINAL_DIR / OUT_NAME
    if final_mp4.exists():
        final_mp4.unlink()
    shutil.move(str(local_mp4), str(final_mp4))
    print(f"Saved {final_mp4}")


if __name__ == "__main__":
    main()