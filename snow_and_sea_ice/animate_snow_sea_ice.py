"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations Library
by Professor John Michael Wallace.

Script: animate_snow_sea_ice.py

Description: Generates the snow cover and sea ice animation from IMS 4 km NH snow/ice
             (NSIDC G02156) and NOAA/NSIDC Sea Ice Concentration CDR V6 south (G02202_V6)
             (2019-2021), rendered as side-by-side north and south polar orthographic globes
             over Blue Marble Next Generation (2004) imagery.

Note: Adapted from NASA SVS 4995 ("Global Snow Cover and Sea Ice Cycle at Both Poles").

Approach: every output pixel of each globe is mapped once (at startup) to its source pixel in
the Blue Marble image, the IMS grid, and the G02202 grid. Each frame is then pure numpy
indexing, so memory stays low and rendering is fast.

"""

import gc
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from PIL import Image
from pyproj import CRS, Transformer
from tqdm import tqdm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation, font_manager
import imageio_ffmpeg

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
START, END = "2019-01-01", "2021-12-31"
FRAME_STEP_DAYS = 5            # 5 -> ~220 frames; 1 -> daily (1,096 frames)
FPS = 6                        # 6 fps at 5-day steps ~ 37 s; for daily frames use ~18
PREVIEW_ONLY = False            # True = write one PNG and stop
PREVIEW_DATE = "2019-08-05"

DRIVE = Path("/Volumes/CLIMATEDATA")
IMS_DIR = DRIVE / "ims_4km"
SIC_DIR = DRIVE / "sea_ice_concentration_south"
BM_DIR = DRIVE / "blue_marble" / "2004"

SCRIPT_DIR = Path.home() / "CLIMATE ANIMATIONS" / "snow_and_sea_ice"
OUT_NAME = "snow_sea_ice_2019_2021.mp4"
LOCAL_OUT = SCRIPT_DIR / OUT_NAME                      # APFS
FINAL_OUT = DRIVE / "snow_and_sea_ice" / OUT_NAME      # exFAT

# Blue Marble month per globe: each hemisphere's snow-minimum month, so all snow
# and ice on screen comes from the 2019-2021 data, not from 2004 imagery.
NORTH_BM_MONTH = 8
SOUTH_BM_MONTH = 2

# Frame layout (1920x1080)
W, H = 1920, 1080
R = 420                                                # globe radius, px
CENTERS = {"north": (510, 555), "south": (1410, 555)}
DATE_XY = (1880, 1045)                                 # lower-right date label (right-aligned)

# Colours / blending
SNOW_RGB = np.array([248, 248, 252], np.float32)
ICE_RGB = np.array([232, 240, 250], np.float32)
SNOW_ALPHA = 0.90
ICE_ALPHA = 0.95
SIC_FULL = 60.0            # concentration (%) at which southern ice is fully opaque
SIC_GAMMA = 0.7
GRID_RGB = np.array([45, 45, 50], np.float32)
GRID_ALPHA = 0.35
HAZE_RGB = np.array([170, 195, 230], np.float32)
HAZE_POWER, HAZE_STRENGTH = 8, 0.55

IMS_VAR = "IMS_Surface_Values"      # 1 water, 2 snow-free land, 3 sea/lake ice, 4 snow on land
SIC_VAR = "cdr_seaice_conc"         # uint8 0-100 %, 255 = land / fill


# ---------------------------------------------------------------------------
# File indexing
# ---------------------------------------------------------------------------
def index_files(root, regex, fmt):
    pat = re.compile(regex)
    files = {}
    for p in root.rglob("*.nc"):
        if p.name.startswith("._"):
            continue
        m = pat.search(p.name)
        if m:
            files[pd.to_datetime(m.group(1), format=fmt)] = p
    if not files:
        sys.exit(f"No files found under {root}")
    return files


def nearest_file(files, d, max_days=3):
    for k in range(max_days + 1):
        for cand in (d - pd.Timedelta(days=k), d + pd.Timedelta(days=k)):
            if cand in files:
                return files[cand]
    raise FileNotFoundError(f"No file within {max_days} days of {d.date()}")


def read_raw(path, var):
    with xr.open_dataset(path, mask_and_scale=False) as ds:
        return ds[var].values[0]


# ---------------------------------------------------------------------------
# Projection helpers
# ---------------------------------------------------------------------------
def grid_crs(ds, crs_var):
    attrs = ds[crs_var].attrs
    try:
        return CRS.from_cf(attrs)
    except Exception:
        pass
    for key in ("crs_wkt", "spatial_ref", "proj4text", "proj4_string", "proj4"):
        if key in attrs:
            return CRS.from_user_input(attrs[key])
    raise ValueError(f"Could not read a CRS from '{crs_var}' attributes: {list(attrs)}")


def grid_fraction(sample_path, crs_var, lat, lon):
    """Fractional (col, row) position of each lat/lon inside the file's x/y grid."""
    with xr.open_dataset(sample_path, mask_and_scale=False) as ds:
        crs = grid_crs(ds, crs_var)
        xs = ds["x"].values.astype(np.float64)
        ys = ds["y"].values.astype(np.float64)
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = tr.transform(lon, lat)
    fx = (np.asarray(x) - xs[0]) / (xs[1] - xs[0])
    fy = (np.asarray(y) - ys[0]) / (ys[1] - ys[0])
    fx = np.where(np.isfinite(fx), fx, -9.0)
    fy = np.where(np.isfinite(fy), fy, -9.0)
    return fx, fy, len(xs), len(ys)


# ---------------------------------------------------------------------------
# Globe geometry (computed once)
# ---------------------------------------------------------------------------
def load_blue_marble(month):
    paths = sorted(p for p in BM_DIR.glob(f"world.2004{month:02d}*.jpg")
                   if not p.name.startswith("._"))
    if not paths:
        sys.exit(f"No Blue Marble image for month {month:02d} in {BM_DIR}")
    return np.asarray(Image.open(paths[0]).convert("RGB"))


def graticule(lat2d, rc, lon2d, hemi):
    """Anti-aliased lat circles (every 10 deg) and meridians (every 30 deg), in exact px distance."""
    sign = 1 if hemi == "north" else -1
    d_lat = np.full(rc.shape, np.inf, np.float32)
    for k in range(10, 90, 10):
        d_lat = np.minimum(d_lat, R * np.abs(rc - np.cos(np.radians(k))))
    dlon = ((lon2d + 15.0) % 30.0) - 15.0
    d_lon = R * rc * np.abs(np.sin(np.radians(dlon)))
    d_lon[sign * lat2d > 80] = np.inf                  # stop meridians before the pole
    d = np.minimum(d_lat, d_lon)
    return np.clip(1.0 - d / 0.8, 0.0, 1.0).astype(np.float32)


def build_globe(hemi, bm_month, ims_sample, sic_sample):
    cx, cy = CENTERS[hemi]
    ys, xs = np.mgrid[cy - R:cy + R, cx - R:cx + R]
    u = (xs + 0.5 - cx) / R
    v = (cy - (ys + 0.5)) / R
    r2d = np.hypot(u, v)
    inside = r2d < 1.0
    rc = np.clip(r2d, 0.0, 1.0)

    if hemi == "north":            # Greenwich meridian pointing down, 90E to the right
        lat2d = np.degrees(np.arccos(rc))
        lon2d = np.degrees(np.arctan2(u, -v))
    else:                          # Greenwich meridian pointing up, 90E to the right
        lat2d = -np.degrees(np.arccos(rc))
        lon2d = np.degrees(np.arctan2(u, v))

    grid2d = graticule(lat2d, rc, lon2d, hemi)

    g = {
        "rows": ys[inside],
        "cols": xs[inside],
        "lat": lat2d[inside],
        "lon": lon2d[inside],
        "grid": grid2d[inside] * GRID_ALPHA,
    }
    r = r2d[inside]
    g["edge"] = np.clip((1.0 - r) * R, 0.0, 1.0).astype(np.float32)
    g["haze"] = (r ** HAZE_POWER * HAZE_STRENGTH).astype(np.float32)

    img = load_blue_marble(bm_month)
    h, w, _ = img.shape
    col = np.clip(((g["lon"] + 180.0) / 360.0 * w).astype(np.int64), 0, w - 1)
    row = np.clip(((90.0 - g["lat"]) / 180.0 * h).astype(np.int64), 0, h - 1)
    g["base"] = img[row, col].astype(np.float32)
    del img
    gc.collect()

    if hemi == "north":
        fx, fy, nx, ny = grid_fraction(ims_sample, "projection", g["lat"], g["lon"])
        ix, iy = np.rint(fx).astype(np.int64), np.rint(fy).astype(np.int64)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        g["ims_sel"] = np.nonzero(ok)[0]
        g["ims_ix"], g["ims_iy"] = ix[ok], iy[ok]
    else:
        sel = np.nonzero(g["lat"] < -30.0)[0]
        fx, fy, nx, ny = grid_fraction(sic_sample, "crs", g["lat"][sel], g["lon"][sel])
        x0, y0 = np.floor(fx).astype(np.int64), np.floor(fy).astype(np.int64)
        ok = (x0 >= 0) & (x0 < nx - 1) & (y0 >= 0) & (y0 < ny - 1)
        g["sic_sel"] = sel[ok]
        g["sic_x0"], g["sic_y0"] = x0[ok], y0[ok]
        g["sic_wx"] = (fx[ok] - x0[ok]).astype(np.float32)
        g["sic_wy"] = (fy[ok] - y0[ok]).astype(np.float32)
    return g


# ---------------------------------------------------------------------------
# Per-frame compositing
# ---------------------------------------------------------------------------
def blend(rgb, color, alpha):
    return rgb * (1.0 - alpha[:, None]) + color * alpha[:, None]


def composite(g, canvas, ice_alpha=None, snow_alpha=None):
    rgb = g["base"].copy()
    if snow_alpha is not None:
        rgb = blend(rgb, SNOW_RGB, snow_alpha)
    if ice_alpha is not None:
        rgb = blend(rgb, ICE_RGB, ice_alpha)
    rgb = blend(rgb, GRID_RGB, g["grid"])
    rgb = blend(rgb, HAZE_RGB, g["haze"])
    rgb *= g["edge"][:, None]
    canvas[g["rows"], g["cols"]] = np.clip(rgb, 0, 255).astype(np.uint8)


def north_layers(g, ims):
    n = g["lat"].size
    vals = ims[g["ims_iy"], g["ims_ix"]]
    snow = np.zeros(n, np.float32)
    ice = np.zeros(n, np.float32)
    snow[g["ims_sel"][vals == 4]] = SNOW_ALPHA
    ice[g["ims_sel"][vals == 3]] = ICE_ALPHA
    return ice, snow


def south_layer(g, sic):
    valid = (sic <= 100).astype(np.float32)             # 255 = land / fill
    c = np.where(sic <= 100, sic, 0).astype(np.float32)
    A = np.clip((c - 15.0) / (SIC_FULL - 15.0), 0.0, 1.0) ** SIC_GAMMA
    x0, y0, wx, wy = g["sic_x0"], g["sic_y0"], g["sic_wx"], g["sic_wy"]
    # Bilinear over ocean cells only, so land cells don't thin the ice along the coast
    num = np.zeros(x0.size, np.float32)
    den = np.zeros(x0.size, np.float32)
    for dy, dx, w in ((0, 0, (1 - wx) * (1 - wy)), (0, 1, wx * (1 - wy)),
                      (1, 0, (1 - wx) * wy), (1, 1, wx * wy)):
        v = valid[y0 + dy, x0 + dx] * w
        num += A[y0 + dy, x0 + dx] * v
        den += v
    a = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    ice = np.zeros(g["lat"].size, np.float32)
    ice[g["sic_sel"]] = a * ICE_ALPHA
    return ice


# ---------------------------------------------------------------------------
# Figure / overlays
# ---------------------------------------------------------------------------
def pick_font():
    names = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"):
        if name in names:
            return name
    return "sans-serif"


def format_date(d):
    return f"{d.strftime('%b')} {d.day}, {d.year}"      # e.g. Aug 5, 2019


def build_figure(canvas):
    font = pick_font()
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor="black")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    im = ax.imshow(canvas, extent=(0, W, H, 0), interpolation="nearest")

    txt = dict(color="white", fontfamily=font, va="center")
    ax.text(W / 2, 58, "Snow Cover & Sea Ice", fontsize=30, fontweight="bold", ha="center", **txt)
    ax.text(CENTERS["north"][0], 112, "Northern Hemisphere", fontsize=20, ha="center", **txt)
    ax.text(CENTERS["south"][0], 112, "Southern Hemisphere", fontsize=20, ha="center", **txt)
    date_text = ax.text(*DATE_XY, "", fontsize=24, ha="right", **txt)
    return fig, im, date_text


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if not DRIVE.exists():
        sys.exit(f"{DRIVE} is not mounted.")

    ims_files = index_files(IMS_DIR, r"ims(\d{7})_", "%Y%j")
    sic_files = index_files(SIC_DIR, r"sic_pss25_(\d{8})_", "%Y%m%d")
    print(f"Indexed {len(ims_files)} IMS files, {len(sic_files)} G02202 south files")

    print("Precomputing globe geometry ...")
    north = build_globe("north", NORTH_BM_MONTH, next(iter(ims_files.values())), None)
    south = build_globe("south", SOUTH_BM_MONTH, None, next(iter(sic_files.values())))

    canvas = np.zeros((H, W, 3), np.uint8)
    fig, im, date_text = build_figure(canvas)

    def draw(d):
        ice_n, snow_n = north_layers(north, read_raw(nearest_file(ims_files, d), IMS_VAR))
        composite(north, canvas, ice_alpha=ice_n, snow_alpha=snow_n)
        ice_s = south_layer(south, read_raw(nearest_file(sic_files, d), SIC_VAR))
        composite(south, canvas, ice_alpha=ice_s)
        im.set_data(canvas)
        date_text.set_text(format_date(d))

    if PREVIEW_ONLY:
        d = pd.Timestamp(PREVIEW_DATE)
        draw(d)
        out = SCRIPT_DIR / f"preview_{d.date()}.png"
        fig.savefig(out, dpi=100, facecolor="black")
        print(f"Preview written: {out}")
        return

    dates = pd.date_range(START, END, freq=f"{FRAME_STEP_DAYS}D")
    print(f"Rendering {len(dates)} frames at {FPS} fps (~{len(dates) / FPS:.0f} s)")

    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    writer = animation.FFMpegWriter(
        fps=FPS, codec="libx264", bitrate=-1,
        extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"],
    )
    anim = animation.FuncAnimation(
        fig, lambda i: draw(dates[i]), frames=len(dates),
        blit=False, cache_frame_data=False,
    )
    pbar = tqdm(total=len(dates), unit="frame")
    anim.save(str(LOCAL_OUT), writer=writer, dpi=100,
              savefig_kwargs={"facecolor": "black"},
              progress_callback=lambda i, n: pbar.update(1))
    pbar.close()
    plt.close(fig)

    n_frames, secs = imageio_ffmpeg.count_frames_and_secs(str(LOCAL_OUT))
    print(f"Verified {LOCAL_OUT.name}: {n_frames} frames, {secs:.1f} s")
    if n_frames < len(dates):
        sys.exit("Frame count lower than expected — not moving to drive.")

    FINAL_OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(LOCAL_OUT), str(FINAL_OUT))
    print(f"Moved to {FINAL_OUT}")


if __name__ == "__main__":
    main()