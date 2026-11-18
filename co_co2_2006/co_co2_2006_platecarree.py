"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: co_co2_2006_platecarree.py
Description: Generates the "CO and CO2, 2006" animation from NASA GMAO MERRA-2
    column CO (COCL, hourly) and CAMS EGG4 column-mean CO2 (XCO2, 3-hourly)
    (2006-01-01 to 2006-12-31), rendered in the Plate Carree projection.
    Column CO2 is drawn in color (377-395 ppmv) and column CO as white plumes
    (0-6 x 1e18 molecules cm-2) over a base map of blue ocean, grey land and
    white Antarctica, following NASA GMAO's GEOS-5 Nature Run visualization
    "Simulated Atmospheric Carbon Concentrations" (SVS 30515).
Note: This is currently the only projection for this variable; the input
    files come from download_co_co2_2006.py in the co_co2_2006 scripts folder.

Frames:
    One frame per 3-hourly CO2 time step (2,920 frames). The hourly CO fields
    whose centers fall within 45 minutes of that time (MERRA-2 hours centered
    at HH:30, so normally the two hours either side) are averaged to match it.
    Both fields are bilinearly interpolated to a 0.25 deg display grid.

Compositing (back to front):
    1) Base map: flat blue ocean; land from the January 2004 Blue Marble,
       desaturated to grey; Antarctica blended to white. Land mask from Natural
       Earth 50m land (cached as _landmask_0.1deg.npy).
    2) CO2: colormap with opacity ramping up from CO2_ALPHA_LO to CO2_ALPHA_HI,
       so low-CO2 air (mostly Southern Hemisphere) lets the base map show.
    3) CO: white with opacity w = clip((CO - CO_FADE_START) / (6 - CO_FADE_START))
       ^ CO_GAMMA, so fire plumes are white and clean air is invisible.

Run:
    conda activate climate
    caffeinate -i python co_co2_2006_platecarree.py
    (Set PREVIEW_ONLY = True first to write 4 test stills in a few seconds.)
"""

import shutil
import sys
import time
from pathlib import Path

import cartopy.crs as ccrs
import dask
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap, Normalize
from tqdm import tqdm

dask.config.set(scheduler="synchronous")
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
YEAR = 2006
WORK_DIR = Path.home() / "CLIMATE ANIMATIONS" / "co_co2_2006"     # local APFS
DRIVE_ROOT = Path("/Volumes/CLIMATEDATA")
DATA_DIR = DRIVE_ROOT / "co_co2" / str(YEAR)                      # exFAT

CO_FILE = f"MERRA2_COCL_hourly_{YEAR}.nc"
CO_VAR = "CO_column"                     # 1e18 molecules cm-2
CO2_GLOB = f"CAMS_EGG4_XCO2*{YEAR}*.nc"
CO2_VAR = "CO2_column"                   # ppmv

BLUE_MARBLE = DRIVE_ROOT / "blue_marble" / "2004" / "world.200401.3x5400x2700.jpg"
LANDMASK_CACHE = WORK_DIR / "_landmask_0.1deg.npy"

OUT_NAME = f"co_co2_{YEAR}_platecarree.mp4"
PREVIEW_ONLY = False                     # True: write 4 PNG stills and stop
PREVIEW_DATES = ["2006-01-20T00", "2006-04-15T12", "2006-08-20T12",
                 "2006-10-10T12"]

# Color scales (the original's colorbars)
CO2_MIN, CO2_MAX = 377.0, 395.0          # ppmv
CO_MIN, CO_MAX = 0.0, 6.0                # 1e18 molecules cm-2
CO_FADE_START = 1.0                      # CO below this (clean air) is invisible
CO_GAMMA = 0.8                           # <1 brightens moderate plumes
CO2_ALPHA_LO, CO2_ALPHA_HI = 379.0, 382.0  # ppmv: transparent -> fully opaque

# Base map
OCEAN_RGB = (0.11, 0.17, 0.42)           # deep blue, close to the 378-379 color
LAND_BRIGHTNESS = 0.80                   # scales Blue Marble luminance
LAND_KEEP_COLOR = 0.15                   # 0 = pure grey, 1 = original colors
ANTARCTIC_WHITE = 0.85                   # blend of Antarctica toward white
BASE_RES = 0.1                           # deg

# Display grid and video
DISPLAY_RES = 0.25                       # deg
FPS = 30                                 # 2,920 frames -> ~97 s
DPI = 100
FIG_W, FIG_H = 19.2, 10.8                # 1920 x 1080
CRF = "18"

# CO2 (the original's colorbar): black 377 -> navy -> teal -> yellow 383
#      -> red 385 -> dark red -> purple -> pink 395
def _stop(v):
    return (v - CO2_MIN) / (CO2_MAX - CO2_MIN)


CO2_CMAP = LinearSegmentedColormap.from_list("gmao_co2", [
    (_stop(377.0), "#14142a"), (_stop(378.5), "#1c2b6c"),
    (_stop(380.0), "#2a6fb3"), (_stop(381.0), "#62c6cf"),
    (_stop(382.0), "#b9e18f"), (_stop(383.0), "#f3e15c"),
    (_stop(384.0), "#f4a443"), (_stop(385.0), "#df4429"),
    (_stop(386.0), "#8f1b1b"), (_stop(387.0), "#5b1f5f"),
    (_stop(388.5), "#9a72c2"), (_stop(390.0), "#d9b6e9"),
    (_stop(392.0), "#efb5df"), (_stop(395.0), "#f8ccee")])
CO_CMAP = LinearSegmentedColormap.from_list("co_white", ["#000000", "#ffffff"])

CO_TOLERANCE = pd.Timedelta(minutes=45)


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------
def _locate(name_or_glob):
    for d in (DATA_DIR, WORK_DIR):
        hits = sorted(p for p in d.glob(name_or_glob)
                      if not p.name.startswith("._")) if d.exists() else []
        if hits:
            return hits[0]
    sys.exit(f"Cannot find {name_or_glob} in {DATA_DIR} or {WORK_DIR}")


def _axis_weights(src, dst, periodic):
    """Bilinear indices/weights mapping a 1-D source axis onto dst."""
    src = np.asarray(src, dtype="float64")
    if periodic:
        src = np.concatenate([[src[-1] - 360.0], src, [src[0] + 360.0]])
        offset = -1                              # padded index -> real index
    else:
        offset = 0
    dst = np.clip(dst, src[0], src[-1])
    i1 = np.clip(np.searchsorted(src, dst, side="right"), 1, len(src) - 1)
    i0 = i1 - 1
    w = (dst - src[i0]) / (src[i1] - src[i0])
    if periodic:
        n = len(src) - 2
        i0 = (i0 + offset) % n
        i1 = (i1 + offset) % n
    return i0, i1, w.astype("float32")


def _regrid(field, lat_w, lon_w):
    """Separable bilinear interpolation of a (lat, lon) array."""
    j0, j1, wx = lon_w
    a = field[:, j0] * (1 - wx) + field[:, j1] * wx
    i0, i1, wy = lat_w
    return a[i0] * (1 - wy)[:, None] + a[i1] * wy[:, None]


class DayCache:
    """Reads a variable one day-sized block at a time (matches file chunks)."""

    def __init__(self, da, per_day):
        self.da, self.per_day, self.blocks = da, per_day, {}

    def get(self, idx):
        day = idx // self.per_day
        if day not in self.blocks:
            if len(self.blocks) >= 2:
                self.blocks.pop(min(self.blocks))
            s = day * self.per_day
            self.blocks[day] = self.da.isel(
                time=slice(s, s + self.per_day)).values.astype("float32")
        return self.blocks[day][idx - day * self.per_day]


# ---------------------------------------------------------------------------
# Compositing
# ---------------------------------------------------------------------------
_co2_norm = Normalize(CO2_MIN, CO2_MAX, clip=True)


def _co2_alpha(co2):
    a = np.clip((co2 - CO2_ALPHA_LO) / (CO2_ALPHA_HI - CO2_ALPHA_LO), 0, 1)
    return a * a * (3 - 2 * a)                     # smoothstep


def _co_weight(co):
    return np.clip((co - CO_FADE_START) / (CO_MAX - CO_FADE_START), 0, 1) \
        ** CO_GAMMA


def composite(co2, co):
    """RGBA overlay: white CO drawn over the semi-transparent CO2 colors."""
    rgb = CO2_CMAP(_co2_norm(co2))[..., :3].astype("float32")
    a2 = _co2_alpha(co2).astype("float32")
    w = _co_weight(co).astype("float32")
    alpha = w + a2 * (1 - w)
    out = np.empty(co2.shape + (4,), dtype="float32")
    out[..., :3] = (w[..., None] + rgb * (a2 * (1 - w))[..., None]) / \
        np.maximum(alpha, 1e-6)[..., None]
    out[..., 3] = alpha
    return out


# ---------------------------------------------------------------------------
# Base map: blue ocean, grey land, white Antarctica
# ---------------------------------------------------------------------------
def _land_mask(nlat, nlon):
    """Land mask on the base grid (row 0 = 90N), cached after the first run."""
    if LANDMASK_CACHE.exists():
        m = np.load(LANDMASK_CACHE)
        if m.shape == (nlat, nlon):
            return m
    lat = 90 - (np.arange(nlat) + 0.5) * 180 / nlat
    lon = -180 + (np.arange(nlon) + 0.5) * 360 / nlon
    LON, LAT = np.meshgrid(lon, lat)
    try:
        import shapely
        from cartopy.io import shapereader
        shp = shapereader.natural_earth("50m", "physical", "land")
        land = shapely.union_all(list(shapereader.Reader(shp).geometries()))
        shapely.prepare(land)
        print("Building land mask from Natural Earth 50m (first run only)...")
        m = shapely.contains_xy(land, LON, LAT)
    except Exception as e:
        print(f"  Natural Earth unavailable ({e}); using a colour-based mask")
        return None
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    np.save(LANDMASK_CACHE, m)
    return m


def build_base():
    from PIL import Image
    nlat, nlon = int(round(180 / BASE_RES)), int(round(360 / BASE_RES))
    bm = np.asarray(Image.open(BLUE_MARBLE).convert("RGB")
                    .resize((nlon, nlat), Image.LANCZOS), dtype="float32") / 255
    lum = bm @ np.array([0.299, 0.587, 0.114], dtype="float32")
    land = _land_mask(nlat, nlon)
    if land is None:                               # Blue Marble ocean is dark blue
        land = ~((bm[..., 2] > bm[..., 0] + 0.03) & (lum < 0.35))
    grey = lum[..., None] * LAND_BRIGHTNESS
    land_rgb = grey * (1 - LAND_KEEP_COLOR) + bm * LAND_BRIGHTNESS * LAND_KEEP_COLOR
    lat = 90 - (np.arange(nlat) + 0.5) * 180 / nlat
    ant = (land & (lat[:, None] < -60))[..., None]
    land_rgb = np.where(ant, land_rgb * (1 - ANTARCTIC_WHITE) + ANTARCTIC_WHITE,
                        land_rgb)
    base = np.where(land[..., None], land_rgb,
                    np.array(OCEAN_RGB, dtype="float32"))
    return np.clip(base, 0, 1)


# ---------------------------------------------------------------------------
def main():
    co_path, co2_path = _locate(CO_FILE), _locate(CO2_GLOB)
    print(f"CO : {co_path}\nCO2: {co2_path}")
    co_ds = xr.open_dataset(co_path, chunks={"time": 24})
    co2_ds = xr.open_dataset(co2_path, chunks={"time": 8})
    co_da, co2_da = co_ds[CO_VAR], co2_ds[CO2_VAR]

    co_t = pd.to_datetime(co_da["time"].values)
    co2_t = pd.to_datetime(co2_da["time"].values)
    print(f"CO : {len(co_t)} hourly steps, {co_t[0]} to {co_t[-1]}")
    print(f"CO2: {len(co2_t)} 3-hourly steps, {co2_t[0]} to {co2_t[-1]}")
    co_per_day = int(round(pd.Timedelta(days=1) / (co_t[1] - co_t[0])))
    co2_per_day = int(round(pd.Timedelta(days=1) / (co2_t[1] - co2_t[0])))

    # Match each CO2 time to the CO hours centered within 45 minutes of it
    matches = []
    for t in co2_t:
        lo = np.searchsorted(co_t, t - CO_TOLERANCE, side="left")
        hi = np.searchsorted(co_t, t + CO_TOLERANCE, side="right")
        if hi <= lo:
            sys.exit(f"No CO hour within 45 min of {t}")
        matches.append(list(range(lo, hi)))
    print(f"Matched {len(matches)} frames "
          f"({sum(len(m) == 2 for m in matches)} averaged from 2 CO hours)")

    # Scale check against the original's colorbars
    s = co2_da.isel(time=slice(0, None, co2_per_day * 7)).values
    p = np.percentile(s, [1, 50, 99])
    print(f"CO2 weekly sample: p1={p[0]:.1f} median={p[1]:.1f} "
          f"p99={p[2]:.1f} ppmv (scale {CO2_MIN:.0f}-{CO2_MAX:.0f})")

    # Display grid and interpolation weights
    lat = np.arange(-90, 90 + DISPLAY_RES / 2, DISPLAY_RES)
    lon = np.arange(-180, 180 + DISPLAY_RES / 2, DISPLAY_RES)
    co_w = (_axis_weights(co_da["lat"].values, lat, False),
            _axis_weights(co_da["lon"].values, lon, True))
    co2_w = (_axis_weights(co2_da["lat"].values, lat, False),
             _axis_weights(co2_da["lon"].values, lon, True))
    co_cache, co2_cache = DayCache(co_da, co_per_day), DayCache(co2_da,
                                                                co2_per_day)

    def frame_rgb(k):
        co = np.mean([co_cache.get(i) for i in matches[k]], axis=0)
        return composite(_regrid(co2_cache.get(k), *co2_w),
                         _regrid(co, *co_w))

    # ---------------- Figure ----------------
    # Full-width 2:1 map (1920 x 960) over a 120 px black legend band
    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI, facecolor="black")
    band = 1 - (FIG_W / 2) / FIG_H
    ax = fig.add_axes([0, band, 1, 1 - band], projection=ccrs.PlateCarree())
    ax.set_global()
    ax.set_facecolor("black")
    ax.spines["geo"].set_visible(False)
    ax.imshow(build_base(), origin="upper", extent=[-180, 180, -90, 90],
              transform=ccrs.PlateCarree(), interpolation="bilinear", zorder=0)
    img = ax.imshow(frame_rgb(0), origin="lower", extent=[-180, 180, -90, 90],
                    transform=ccrs.PlateCarree(), interpolation="bilinear",
                    zorder=1)

    # Date block (left)
    stamp = fig.text(0.125, band * 0.60, "", color="white", fontsize=40,
                     fontweight="bold", ha="left", va="center")
    hour = fig.text(0.0, band * 0.60, "", color="#bbbbbb", fontsize=20,
                    ha="left", va="center")
    fig.text(0.125, band * 0.15, "SDSU Climate Informatics Lab", color="white",
             fontsize=13, fontweight="bold", ha="left", va="center")

    # Colorbars (right): CO grey, then CO2 as it appears over black
    def _band_cbar(x0, width, cmap, norm, ticks, label, fmt):
        cax = fig.add_axes([x0, band * 0.38, width, band * 0.20])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm, cmap), cax=cax,
                          orientation="horizontal")
        cb.set_ticks(ticks)
        cb.ax.set_xticklabels([fmt.format(t) for t in ticks])
        cb.ax.tick_params(colors="white", labelsize=9, length=0, pad=3)
        cb.outline.set_visible(False)
        fig.text(x0, band * 0.80, label, color="white", fontsize=11,
                 ha="left", va="center")

    gamma_cmap = LinearSegmentedColormap.from_list(
        "co_gamma", CO_CMAP(_co_weight(np.linspace(CO_MIN, CO_MAX, 256))))
    _band_cbar(0.385, 0.29, gamma_cmap, Normalize(CO_MIN, CO_MAX),
               np.arange(CO_MIN, CO_MAX + 0.01, 0.6),
               "Carbon Monoxide Column Abundance [10$^{18}$ molec cm$^{-2}$]",
               "{:.1f}")
    v = np.linspace(CO2_MIN, CO2_MAX, 256)
    co2_bar = CO2_CMAP(_co2_norm(v))
    co2_bar[:, :3] *= _co2_alpha(v)[:, None]       # as seen over black
    _band_cbar(0.695, 0.29, LinearSegmentedColormap.from_list("co2_bar",
                                                              co2_bar),
               _co2_norm, np.arange(CO2_MIN, CO2_MAX + 0.1, 2),
               "Carbon Dioxide Column Concentration [ppmv]", "{:.0f}")

    # Place the hour label just right of the widest date of the year
    renderer = fig.canvas.get_renderer()
    x1 = 0.0
    for d in sorted({f"{t:%Y / %m / %d}" for t in co2_t}):
        stamp.set_text(d)
        bb = stamp.get_window_extent(renderer).transformed(
            fig.transFigure.inverted())
        x1 = max(x1, bb.x1)
    hour.set_x(x1 + 0.008)

    def draw(k):
        img.set_data(frame_rgb(k))
        stamp.set_text(f"{co2_t[k]:%Y / %m / %d}")
        hour.set_text(f"{co2_t[k]:%H}z")

    # ---------------- Preview stills ----------------
    if PREVIEW_ONLY:
        for d in PREVIEW_DATES:
            k = int(np.argmin(np.abs(co2_t - pd.Timestamp(d))))
            draw(k)
            out = WORK_DIR / f"preview_{co2_t[k]:%Y%m%d_%H}.png"
            fig.savefig(out, dpi=DPI, facecolor="black")
            print(f"  wrote {out}")
        return

    # ---------------- Video ----------------
    local_out = WORK_DIR / OUT_NAME
    tmp = local_out.with_name(local_out.stem + "_tmp.mp4")
    writer = FFMpegWriter(fps=FPS, codec="libx264", bitrate=-1,
                          extra_args=["-pix_fmt", "yuv420p", "-crf", CRF,
                                      "-preset", "veryfast"])
    t0 = time.time()
    with writer.saving(fig, str(tmp), DPI):
        for k in tqdm(range(len(co2_t)), desc="Rendering", unit="frame"):
            draw(k)
            writer.grab_frame(facecolor="black")
    plt.close(fig)
    co_ds.close()
    co2_ds.close()
    tmp.replace(local_out)
    print(f"Rendered {len(co2_t)} frames in {(time.time() - t0) / 60:.1f} min "
          f"-> {local_out} ({local_out.stat().st_size / 1e6:.0f} MB)")

    if DRIVE_ROOT.exists():
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        final = DATA_DIR / OUT_NAME
        if final.exists():
            final.unlink()
        shutil.move(str(local_out), str(final))
        print(f"Moved to {final}")
    else:
        print(f"{DRIVE_ROOT} not mounted - video left at {local_out}")


if __name__ == "__main__":
    main()
