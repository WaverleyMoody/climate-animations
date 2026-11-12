"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: arctic_co_orthographic.py
Description: Generates the Arctic wildfire carbon monoxide animation from NASA
    MERRA-2 hourly CO total column COCL (2019-06-01 to 2019-08-31), rendered in
    the Orthographic projection.
Note: For the data download, see download_merra2_co_2019_jja.py in the
    C02_arctic scripts folder.

Background removal:
    The original (GEOS-5 FP) plots CO tagged as biomass burning, so clean air
    is ~0 DU. MERRA-2 only has total CO, where clean high-latitude air is
    ~40-55 DU. To isolate the fire signal, each pixel's background is the
    BG_PERCENTILE-th percentile of its daily means within +/- BG_HALF_WINDOW
    days, interpolated to each hour and subtracted. The background is computed
    once and cached as a .npy file in LOCAL_DIR.

Run:
    conda activate climate
    # 1) Set PREVIEW_ONLY = True, check the PNG and the printed percentiles
    # 2) Set PREVIEW_ONLY = False for the full render
    caffeinate -i python arctic_co_orthographic.py
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
VARIABLE = "CO_column"                    # DU, from download_merra2_co_2019_jja.py
DATA_PATH = Path("/Volumes/CLIMATEDATA/co/merra2_2019_jja/MERRA2_COCL_hourly_2019-06-08.nc")

LOCAL_DIR = Path.home() / "CLIMATE ANIMATIONS" / "C02_arctic" / "_render"   # local APFS
FINAL_DIR = Path("/Volumes/CLIMATEDATA/co/animations")                      # exFAT
OUT_NAME = "arctic_co_jja2019_orthographic.mp4"

PREVIEW_ONLY = False                       # True: save one PNG and stop
PREVIEW_TIME = "2019-07-08T21:30"         # matches the reference frame (21z)

# Background removal (see docstring). SUBTRACT_BACKGROUND = False plots total CO.
SUBTRACT_BACKGROUND = True
BG_PERCENTILE = 10
BG_HALF_WINDOW = 15                       # days -> 31-day window
BG_CACHE = LOCAL_DIR / f"co_background_p{BG_PERCENTILE}_w{2 * BG_HALF_WINDOW + 1}.npy"

# Globe view (same as the smoke animation)
CENTRAL_LON = -115.0
CENTRAL_LAT = 70.0
LAT_CIRCLES = [63, 80]

GRID_N = 1400

# Colorbar levels from the original (extended both ends), Dobson units
LEVELS = [5, 7, 11, 16, 23, 33, 48, 70, 100, 150]

BG_COLOR = "#f7f9fd"                      # near-white globe background
CMAP_ANCHORS = ["#e6eef9", "#c6daf0", "#9ec5e6", "#72a8d9",
                "#4a89c8", "#2e6db4", "#1b5299", "#0f3a78"]
OVER_COLOR = "#08245a"

LIGHT_LINE = "#9a9a9a"
DARK_LINE = "#3a3a3a"
DARK_COUNTRIES = ["USA", "CAN"]

SHOW_US_CITIES = True
CITY_MIN_POP = 50_000

TITLE = "Carbon Monoxide from Arctic Wildfires, June–August 2019"
SUBTITLE = ("NASA MERRA-2 hourly CO column above background"
            if SUBTRACT_BACKGROUND else "NASA MERRA-2 hourly total CO column")
UNITS_LABEL = "Dobson units"

FPS = 35                                  # 2208 frames -> 63 s
FIGSIZE = (12, 10.4)                      # 1800 x 1560 px at DPI 150
DPI = 150
BITRATE = 12000

dask.config.set(scheduler="synchronous")
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_cmap():
    smooth = LinearSegmentedColormap.from_list("co", CMAP_ANCHORS)
    n_bins = len(LEVELS) - 1
    cmap = ListedColormap(smooth(np.linspace(0, 1, n_bins)), name="co_binned")
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
    """Precompute bilinear lookups from the lat-lon grid onto the globe image."""

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


class Background:
    """Rolling low-percentile background, one field per day, interpolated hourly.

    Daily means are built one day (24 frames, ~20 MB) at a time, so peak memory
    stays near 0.2 GB. The result (92 x 361 x 576 float32, ~76 MB) is cached.
    """

    def __init__(self, da):
        days = pd.DatetimeIndex(da.time.values).normalize().unique()
        self.day_centers = (days + pd.Timedelta(hours=12)).values

        if BG_CACHE.exists():
            self.fields = np.load(BG_CACHE)
            if self.fields.shape != (days.size, da.sizes["lat"], da.sizes["lon"]):
                raise ValueError(f"{BG_CACHE.name} has the wrong shape; delete it and rerun")
            print(f"Loaded background cache {BG_CACHE.name}")
            return

        print(f"Building daily means for {days.size} days...")
        day_index = pd.DatetimeIndex(da.time.values).normalize()
        daily = np.empty((days.size, da.sizes["lat"], da.sizes["lon"]), dtype=np.float32)
        for k, d in enumerate(days):
            idx = np.flatnonzero(day_index == d)
            daily[k] = da.isel(time=idx).mean("time").values
            if k % 10 == 0:
                print(f"  day {k + 1}/{days.size}")

        print(f"Rolling {BG_PERCENTILE}th percentile, +/-{BG_HALF_WINDOW} days...")
        self.fields = np.empty_like(daily)
        for k in range(days.size):
            lo, hi = max(0, k - BG_HALF_WINDOW), min(days.size, k + BG_HALF_WINDOW + 1)
            self.fields[k] = np.percentile(daily[lo:hi], BG_PERCENTILE, axis=0)

        np.save(BG_CACHE, self.fields)
        np.load(BG_CACHE, mmap_mode="r")[-1].sum()          # read-back check
        print(f"Saved background cache {BG_CACHE.name}")

    def at(self, t):
        x = (np.datetime64(t) - self.day_centers[0]) / np.timedelta64(1, "D")
        x = float(np.clip(x, 0, len(self.fields) - 1))
        k = min(int(np.floor(x)), len(self.fields) - 2)
        w = np.float32(x - k)
        return (1 - w) * self.fields[k] + w * self.fields[k + 1]


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

    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.45,
                   edgecolor=LIGHT_LINE, zorder=2)
    ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="none",
                   linewidth=0.35, edgecolor=LIGHT_LINE, zorder=2)

    countries = [g for a, g in _records("50m", "cultural", "admin_0_countries")
                 if a.get("ADM0_A3") in DARK_COUNTRIES]
    ax.add_geometries(countries, crs=pc, facecolor="none", edgecolor=DARK_LINE,
                      linewidth=0.6, zorder=3)
    states = [g for a, g in _records("50m", "cultural", "admin_1_states_provinces_lakes")
              if a.get("ADM0_A3") == "USA"]
    ax.add_geometries(states, crs=pc, facecolor="none", edgecolor=DARK_LINE,
                      linewidth=0.45, zorder=3)

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
    cb.ax.set_yticklabels([f"{v:g}" for v in LEVELS], fontsize=13)
    fig.text(0.904, 0.235, UNITS_LABEL, ha="center", va="top",
             fontsize=13, style="italic")

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

    ds = xr.open_dataset(DATA_PATH, chunks={"time": 24})
    da = ds[VARIABLE].transpose("time", "lat", "lon")
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.isel(lat=slice(None, None, -1))
    if not np.all(np.diff(da.lon.values) > 0):
        sys.exit("Longitude axis is not ascending; check the input file")

    times = da.time.values
    n = times.size
    print(f"Loaded {VARIABLE}: {n} hourly frames, {da.sizes['lat']} x {da.sizes['lon']}")

    bg = Background(da) if SUBTRACT_BACKGROUND else None

    def frame(i):
        f = da.isel(time=i).values
        return f - bg.at(times[i]) if bg is not None else f

    proj = ccrs.Orthographic(central_longitude=CENTRAL_LON, central_latitude=CENTRAL_LAT)
    print("Precomputing orthographic lookup grid...")
    resampler = OrthoResampler(proj, da.lon.values, da.lat.values, GRID_N)

    fig, im, stamp = _build_figure(resampler, proj)

    if PREVIEW_ONLY:
        i = int(np.argmin(np.abs(times - np.datetime64(PREVIEW_TIME))))
        f = frame(i)
        nh = f[da.lat.values >= 30]
        p = np.percentile(nh, [5, 50, 90, 99, 99.9])
        print(f"30-90N at {_stamp(times[i])}: p5={p[0]:.1f}  median={p[1]:.1f}  "
              f"p90={p[2]:.1f}  p99={p[3]:.1f}  p99.9={p[4]:.1f}  max={nh.max():.1f} DU")
        im.set_data(resampler(f))
        stamp.set_text(_stamp(times[i]))
        out = LOCAL_DIR / "preview_frame.png"
        fig.savefig(out, dpi=DPI)
        print(f"Preview saved: {out}")
        return

    local_mp4 = LOCAL_DIR / OUT_NAME
    writer = FFMpegWriter(fps=FPS, codec="libx264", bitrate=BITRATE,
                          extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, str(local_mp4), dpi=DPI):
        for i in range(n):
            im.set_data(resampler(frame(i)))
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
