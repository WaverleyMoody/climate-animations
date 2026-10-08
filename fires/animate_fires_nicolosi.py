"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation Animations
Library by Professor John Michael Wallace.

Script: animate_fires_nicolosi.py

Description: Generates the fire carbon emissions time-series animation from
             GFED5.1 weekly totals (March 2003 - January 2022), rendered as two
             Nicolosi globular hemispheres (west centered on 90W, east on 90E) in
             the style of the NASA SVS "Fires - a global perspective" animation.

Note: For the Plate Carree, Robinson, and Foucaut projections, see the other
      scripts in the fires scripts folder.

Run inside the climate env:  conda activate climate
  1. RENDER_MODE = "test"  -> writes a few PNG frames for review (fast)
  2. RENDER_MODE = "full"  -> caffeinate -i python animate_fires_nicolosi.py
If you change any basemap colors or the layout, delete BG_CACHE.
"""

import io
import shutil
import sys
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import cartopy.feature as cfeature          # used only to fetch Natural Earth geometries
import imageio_ffmpeg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import requests
import shapely
import xarray as xr
from matplotlib import patheffects
from matplotlib.animation import FFMpegWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
from PIL import Image
from shapely.geometry.polygon import orient
from tqdm import tqdm

# ---------------------------------------------------------------- settings
RENDER_MODE = "full"            # "test" or "full"
TEST_DATES = ["2004-03-06", "2010-08-14", "2015-09-26", "2019-12-28"]

WORK = Path.home() / "CLIMATE ANIMATIONS" / "fires"                    # local APFS
WEEKLY = Path("/Volumes/CLIMATEDATA/fires/GFED5.1_weekly")
BASEMAP_DIR = Path("/Volumes/CLIMATEDATA/fires/basemap")
STAGING = WORK / "_staging"
TEST_DIR = WORK / "test_frames"
BG_CACHE = WORK / "fires_background_nicolosi.png"
FINAL_VIDEO = Path("/Volumes/CLIMATEDATA/fires_2003_2022_nicolosi.mp4")

WIDTH, HEIGHT, DPI = 1920, 1080, 100
FPS = 12                        # 987 weeks -> ~82 s (NASA original runs 1:24)
EARTH_RADIUS = 6371007.2        # m (sphere)
HEMI_CENTERS = (-90.0, 90.0)    # central longitudes: western, eastern hemisphere
HEMI_GAP = 24                   # px between the two circles
MAP_MARGIN_X = 24               # min px on each side
MAP_MARGIN_Y = 24               # min px top and bottom
SUPERSAMPLE = 4                 # lon/lat samples per data cell, per axis, for the pixel table
LANDMASK_RES = 0.05             # degrees; resolution of the rasterized land mask

# Basemap colors (NASA SVS look)
FRAME_COLOR = "#1c1e22"             # outside the map outline
OCEAN_COLOR = "#4d5868"
LAND_DARK, LAND_LIGHT = 0.10, 0.36   # gray range for shaded relief (0 = black, 1 = white)
LAND_FLAT = "#2e2e2e"               # used only if the relief raster is unavailable
BORDER_COLOR = "#d9d9d9"
OUTLINE_COLOR = "#8b919b"
RELIEF_URLS = [
    "https://naciscdn.org/naturalearth/50m/raster/GRAY_50M_SR.zip",
    "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/GRAY_50M_SR.zip",
]

# Data colors: green -> orange -> white, 0-20 g C m-2 week-1, saturating above 20
VMIN, VMAX = 0.0, 20.0
CMAP = LinearSegmentedColormap.from_list("nasa_fire", [
    (0.00, "#5e8c3a"),
    (0.25, "#a3b257"),
    (0.50, "#e5b260"),
    (0.75, "#f0d5a4"),
    (1.00, "#f3f1ec"),
])
CMAP.set_over("#f3f1ec")
NORM = Normalize(VMIN, VMAX)
ALPHA_MIN, ALPHA_FULL = 0.05, 1.0   # transparent below ALPHA_MIN, opaque above ALPHA_FULL
TEXT_COLOR = "#e6e6e6"
TEXT_HALO = [patheffects.withStroke(linewidth=3, foreground="#1c1e22")]  # legibility over land
CREDIT = "Data: GFED5.1 (van der Werf et al.)  |  SDSU Climate Informatics Lab"



# ---------------------------------------------------------------- projection geometry
class MapGeometry:
    """
    Two Nicolosi hemispheres side by side. Forward projection from lon/lat to
    continuous pixel coordinates (col, row from the top-left corner).
    """

    def __init__(self):
        self.fwd = [pyproj.Transformer.from_crs(
            "EPSG:4326", f"+proj=nicol +lon_0={lon0} +R={EARTH_RADIUS} +type=crs", always_xy=True)
            for lon0 in HEMI_CENTERS]
        r_m, _ = self.fwd[1].transform(HEMI_CENTERS[1] + 90.0, 0.0)    # hemisphere radius (m)
        d = min((WIDTH - 2 * MAP_MARGIN_X - HEMI_GAP) / 2, HEIGHT - 2 * MAP_MARGIN_Y)
        self.r_px = d / 2
        self.scale = self.r_px / r_m
        x_start = (WIDTH - (2 * d + HEMI_GAP)) / 2
        self.centers = [(x_start + d / 2, HEIGHT / 2), (x_start + 1.5 * d + HEMI_GAP, HEIGHT / 2)]

    @staticmethod
    def hemisphere(lon) -> np.ndarray:
        """0 = western (lon < 0), 1 = eastern (lon >= 0)."""
        return (np.asarray(lon) >= 0).astype(np.int8)

    def lonlat_to_px(self, lon, lat):
        lon = np.asarray(lon, dtype=np.float64)
        lat = np.asarray(lat, dtype=np.float64)
        hemi = self.hemisphere(lon)
        col = np.full(lon.shape, np.nan)
        row = np.full(lon.shape, np.nan)
        for k, (cx, cy) in enumerate(self.centers):
            m = hemi == k
            if m.any():
                x, y = self.fwd[k].transform(lon[m], lat[m])
                col[m] = cx + np.asarray(x) * self.scale
                row[m] = cy - np.asarray(y) * self.scale
        return col, row

    def inside_mask(self) -> np.ndarray:
        cc, rr = np.meshgrid(np.arange(WIDTH) + 0.5, np.arange(HEIGHT) + 0.5)
        inside = np.zeros((HEIGHT, WIDTH), dtype=bool)
        for cx, cy in self.centers:
            inside |= (cc - cx) ** 2 + (rr - cy) ** 2 <= self.r_px ** 2
        return inside


def build_pixel_lonlat(geom: MapGeometry, cell_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """
    (HEIGHT, WIDTH) float32 arrays of lon and lat for every pixel inside the map
    (NaN outside), built by forward-projecting a dense lon/lat sample grid.
    """
    step = cell_deg / SUPERSAMPLE
    lons = -180 + (np.arange(int(round(360 / step))) + 0.5) * step
    lats = -90 + (np.arange(int(round(180 / step))) + 0.5) * step
    pix_lon = np.full(HEIGHT * WIDTH, np.nan, dtype=np.float32)
    pix_lat = np.full(HEIGHT * WIDTH, np.nan, dtype=np.float32)

    for i0 in tqdm(range(0, len(lats), 160), desc="Pixel table", unit="band", leave=False):
        LO, LA = np.meshgrid(lons, lats[i0:i0 + 160])
        LO, LA = LO.ravel(), LA.ravel()
        col, row = geom.lonlat_to_px(LO, LA)
        ok = np.isfinite(col) & np.isfinite(row)
        c = np.floor(col[ok]).astype(np.int64)
        r = np.floor(row[ok]).astype(np.int64)
        inb = (c >= 0) & (c < WIDTH) & (r >= 0) & (r < HEIGHT)
        idx = r[inb] * WIDTH + c[inb]
        pix_lon[idx] = LO[ok][inb]
        pix_lat[idx] = LA[ok][inb]
    pix_lon = pix_lon.reshape(HEIGHT, WIDTH)
    pix_lat = pix_lat.reshape(HEIGHT, WIDTH)

    # Fill the few pixels inside the circles that no sample landed in
    inside = geom.inside_mask()
    holes0 = int((inside & np.isnan(pix_lon)).sum())
    for _ in range(20):
        holes = inside & np.isnan(pix_lon)
        if not holes.any():
            break
        for dr, dc in ((0, 1), (0, -1), (1, 0), (-1, 0)):
            src_lon = np.roll(pix_lon, (dr, dc), axis=(0, 1))
            src_lat = np.roll(pix_lat, (dr, dc), axis=(0, 1))
            take = holes & np.isfinite(src_lon)
            pix_lon[take], pix_lat[take] = src_lon[take], src_lat[take]
            holes &= ~take
    left = int((inside & np.isnan(pix_lon)).sum())
    print(f"Pixel table: {int(np.isfinite(pix_lon).sum()):,} map pixels "
          f"({holes0:,} gap pixels filled, {left} left)")
    return pix_lon, pix_lat


# ---------------------------------------------------------------- basemap
def load_relief() -> np.ndarray | None:
    """Natural Earth 1:50m gray shaded relief as uint8 (lat, lon), north up; None if unavailable."""
    BASEMAP_DIR.mkdir(parents=True, exist_ok=True)
    tif = BASEMAP_DIR / "GRAY_50M_SR.tif"
    if not tif.exists():
        for url in RELIEF_URLS:
            try:
                print(f"Downloading shaded relief: {url}")
                r = requests.get(url, timeout=120)
                r.raise_for_status()
                with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
                    member = next(m for m in zf.namelist()
                                  if m.lower().endswith(".tif") and not Path(m).name.startswith("._"))
                    tif.write_bytes(zf.read(member))
                break
            except Exception as e:  # noqa: BLE001
                print(f"  failed ({type(e).__name__}: {e})")
    if not tif.exists():
        print("Shaded relief unavailable; using flat land color.")
        return None
    Image.MAX_IMAGE_PIXELS = None
    arr = np.asarray(Image.open(tif))
    if arr.ndim == 3:
        arr = arr[..., :3].mean(axis=2)
    return arr[::2, ::2].astype(np.uint8)          # 10800x5400 -> 5400x2700


def _polygon_path(geoms) -> MplPath | None:
    """Shapely (multi)polygons in lon/lat -> one compound Matplotlib path (holes kept)."""
    verts, codes = [], []
    stack = list(geoms)
    while stack:
        g = stack.pop()
        if g.geom_type in ("MultiPolygon", "GeometryCollection"):
            stack.extend(g.geoms)
            continue
        if g.geom_type != "Polygon" or g.is_empty:
            continue
        g = orient(g, 1.0)
        for ring in (g.exterior, *g.interiors):
            xy = np.asarray(ring.coords)[:, :2]
            verts.append(xy)
            codes += [MplPath.MOVETO] + [MplPath.LINETO] * (len(xy) - 2) + [MplPath.CLOSEPOLY]
    return MplPath(np.concatenate(verts), codes) if verts else None


def land_fraction_raster() -> np.ndarray:
    """Antialiased land mask (1 = land, 0 = ocean/lake) on a regular lon/lat grid, north up."""
    nx, ny = int(round(360 / LANDMASK_RES)), int(round(180 / LANDMASK_RES))
    fig = plt.figure(figsize=(nx / DPI, ny / DPI), dpi=DPI, facecolor="black")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    for name, color in (("land", "white"), ("lakes", "black")):
        path = _polygon_path(cfeature.NaturalEarthFeature("physical", name, "50m").geometries())
        if path is not None:
            ax.add_patch(PathPatch(path, facecolor=color, edgecolor="none", antialiased=True))
    fig.canvas.draw()
    mask = np.asarray(fig.canvas.buffer_rgba())[..., 0].astype(np.float32) / 255.0
    plt.close(fig)
    return mask


def sample_bilinear(grid: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Bilinear sample of a north-up global lon/lat grid at the given points."""
    ny, nx = grid.shape
    fx = (lon + 180.0) / 360.0 * nx - 0.5
    fy = (90.0 - lat) / 180.0 * ny - 0.5
    x0 = np.floor(fx).astype(np.int64)
    y0 = np.floor(fy).astype(np.int64)
    tx, ty = fx - x0, fy - y0
    x0c, x1c = np.mod(x0, nx), np.mod(x0 + 1, nx)               # wrap in longitude
    y0c, y1c = np.clip(y0, 0, ny - 1), np.clip(y0 + 1, 0, ny - 1)
    top = grid[y0c, x0c] * (1 - tx) + grid[y0c, x1c] * tx
    bot = grid[y1c, x0c] * (1 - tx) + grid[y1c, x1c] * tx
    return top * (1 - ty) + bot * ty


def border_segments(geom: MapGeometry) -> list[np.ndarray]:
    """Country borders forward-projected to pixel coordinates, split where they change hemisphere."""
    segs = []
    for g in cfeature.NaturalEarthFeature("cultural", "admin_0_boundary_lines_land", "50m").geometries():
        g = shapely.segmentize(g, 0.25)
        lines = g.geoms if g.geom_type.startswith("Multi") else [g]
        for line in lines:
            xy = np.asarray(line.coords)[:, :2]
            col, row = geom.lonlat_to_px(xy[:, 0], xy[:, 1])
            pts = np.column_stack([col, row])
            cuts = np.flatnonzero(np.diff(geom.hemisphere(xy[:, 0])) != 0) + 1
            for part in np.split(pts, cuts):
                if len(part) >= 2 and np.isfinite(part).all():
                    segs.append(part)
    return segs


def build_background(geom: MapGeometry, pix_lon: np.ndarray, pix_lat: np.ndarray) -> np.ndarray:
    """Render the static basemap once; returns a (HEIGHT, WIDTH, 3) float array."""
    if BG_CACHE.exists():
        return plt.imread(BG_CACHE)[..., :3]

    print("Building background ...")
    valid = np.isfinite(pix_lon)
    lon, lat = pix_lon[valid].astype(np.float64), pix_lat[valid].astype(np.float64)

    land = sample_bilinear(land_fraction_raster(), lon, lat)[:, None]
    relief = load_relief()
    if relief is not None:
        ny, nx = relief.shape
        ir = np.clip(((90.0 - lat) / 180.0 * ny).astype(np.int64), 0, ny - 1)
        ic = np.clip(((lon + 180.0) / 360.0 * nx).astype(np.int64), 0, nx - 1)
        gray = LAND_DARK + (LAND_LIGHT - LAND_DARK) * (relief[ir, ic].astype(np.float32) / 255.0)
        land_rgb = np.repeat(gray[:, None], 3, axis=1)
    else:
        land_rgb = np.broadcast_to(np.array(to_rgb(LAND_FLAT)), (lon.size, 3))

    img = np.empty((HEIGHT, WIDTH, 3), dtype=np.float32)
    img[:] = to_rgb(FRAME_COLOR)
    img[valid] = np.array(to_rgb(OCEAN_COLOR)) * (1 - land) + land_rgb * land

    # Vector layers: borders and the two hemisphere outlines
    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI, facecolor=FRAME_COLOR)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.imshow(img, extent=(0, WIDTH, HEIGHT, 0), interpolation="nearest")
    ax.add_collection(LineCollection(border_segments(geom), colors=BORDER_COLOR,
                                     linewidths=0.45, alpha=0.6))
    for cx, cy in geom.centers:
        ax.add_patch(plt.Circle((cx, cy), geom.r_px, fill=False,
                                edgecolor=OUTLINE_COLOR, linewidth=0.8))
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(HEIGHT, 0)
    fig.canvas.draw()
    bg = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)

    WORK.mkdir(parents=True, exist_ok=True)
    plt.imsave(BG_CACHE, bg)
    return bg.astype(np.float32) / 255.0


# ---------------------------------------------------------------- frames
def build_lookup(pix_lon: np.ndarray, pix_lat: np.ndarray,
                 lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Flat pixel indices inside the map, and the flat source-cell index for each."""
    pix = np.flatnonzero(np.isfinite(pix_lon.ravel()))
    plon, plat = pix_lon.ravel()[pix], pix_lat.ravel()[pix]
    dlat, dlon = float(lat[1] - lat[0]), float(lon[1] - lon[0])
    ilat = np.clip(np.rint((plat - lat[0]) / dlat), 0, len(lat) - 1).astype(np.int64)
    ilon = np.clip(np.rint((plon - lon[0]) / dlon), 0, len(lon) - 1).astype(np.int64)
    return pix, ilat * len(lon) + ilon


def fire_rgba_u8(field: np.ndarray) -> np.ndarray:
    """Weekly field -> flat (cells, 4) uint8 RGBA with an alpha ramp for near-zero cells."""
    rgba = CMAP(NORM(field), bytes=True)
    alpha = np.clip((field - ALPHA_MIN) / (ALPHA_FULL - ALPHA_MIN), 0.0, 1.0)
    rgba[..., 3] = (alpha * 255).astype(np.uint8)
    return rgba.reshape(-1, 4)


class Projector:
    """Turns a weekly field into a full-frame RGBA overlay using the pixel lookup."""

    def __init__(self, pix: np.ndarray, src: np.ndarray):
        self.pix, self.src = pix, src
        self.buf = np.zeros((HEIGHT * WIDTH, 4), dtype=np.uint8)

    def __call__(self, field: np.ndarray) -> np.ndarray:
        self.buf[self.pix] = fire_rgba_u8(field)[self.src]
        return self.buf.reshape(HEIGHT, WIDTH, 4)


def setup_figure(bg: np.ndarray):
    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI, facecolor=FRAME_COLOR)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.imshow(bg, origin="upper", aspect="auto", interpolation="nearest", zorder=0)
    fire_im = ax.imshow(np.zeros((HEIGHT, WIDTH, 4), dtype=np.uint8), origin="upper",
                        aspect="auto", interpolation="nearest", zorder=1)
    ax.set_xlim(-0.5, WIDTH - 0.5)
    ax.set_ylim(HEIGHT - 0.5, -0.5)

    # Legend block, bottom center, in the open space between the two hemispheres
    bar_w = 0.20
    cax = fig.add_axes([0.5 - bar_w / 2, 0.075, bar_w, 0.022])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                      orientation="horizontal", extend="max", ticks=[0, 5, 10, 15, 20])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=TEXT_COLOR, labelsize=14, length=4)
    cb.set_label("grams of carbon per square meter per week", color=TEXT_COLOR,
                 fontsize=14, labelpad=6, path_effects=TEXT_HALO)
    for lbl in cb.ax.get_xticklabels():
        lbl.set_path_effects(TEXT_HALO)
    xc = 0.5
    date_txt = fig.text(xc, 0.150, "", ha="center", va="bottom", color=TEXT_COLOR,
                        fontsize=22, fontweight="bold", path_effects=TEXT_HALO)
    fig.text(xc, 0.112, "Carbon Emissions", ha="center", va="bottom",
             color=TEXT_COLOR, fontsize=18, path_effects=TEXT_HALO)
    fig.text(0.99, 0.012, CREDIT, ha="right", va="bottom", color="#bdbdbd", fontsize=11)
    return fig, fire_im, date_txt


def week_index() -> list[tuple[Path, int, pd.Timestamp]]:
    files = sorted(p for p in WEEKLY.glob("GFED5.1_C_weekly_*.nc") if not p.name.startswith("._"))
    if not files:
        sys.exit(f"No weekly files in {WEEKLY}")
    index = []
    for p in files:
        with xr.open_dataset(p) as ds:
            index += [(p, k, t) for k, t in enumerate(pd.DatetimeIndex(ds["time"].values))]
    return index


# ---------------------------------------------------------------- main
def main() -> None:
    if not Path("/Volumes/CLIMATEDATA").is_dir():
        sys.exit("External drive /Volumes/CLIMATEDATA is not mounted.")
    WORK.mkdir(parents=True, exist_ok=True)
    STAGING.mkdir(parents=True, exist_ok=True)

    index = week_index()
    print(f"{len(index)} weekly frames: {index[0][2].date()} to {index[-1][2].date()}")
    with xr.open_dataset(index[0][0]) as ds:
        lat, lon = ds["lat"].values, ds["lon"].values

    geom = MapGeometry()
    pix_lon, pix_lat = build_pixel_lonlat(geom, abs(float(lon[1] - lon[0])))
    bg = build_background(geom, pix_lon, pix_lat)
    project = Projector(*build_lookup(pix_lon, pix_lat, lat, lon))
    del pix_lon, pix_lat
    fig, fire_im, date_txt = setup_figure(bg)

    if RENDER_MODE == "test":
        TEST_DIR.mkdir(parents=True, exist_ok=True)
        times = pd.DatetimeIndex([t for _, _, t in index])
        for d in TEST_DATES:
            i = int(np.abs(times - pd.Timestamp(d)).argmin())
            path, k, t = index[i]
            with xr.open_dataset(path) as ds:
                field = ds["C_weekly"].isel(time=k).values
            fire_im.set_data(project(field))
            date_txt.set_text(t.strftime("%b  %Y"))
            out = TEST_DIR / f"fires_nicolosi_{t.date()}.png"
            fig.savefig(out, dpi=DPI, facecolor=fig.get_facecolor())
            print(f"  wrote {out}")
        plt.close(fig)
        print("Review the PNGs, then set RENDER_MODE = 'full'.")
        return

    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    writer = FFMpegWriter(fps=FPS, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    tmp = STAGING / FINAL_VIDEO.name
    current_path, ds = None, None
    with writer.saving(fig, str(tmp), dpi=DPI):
        for path, k, t in tqdm(index, desc="Rendering", unit="frame"):
            if path != current_path:
                if ds is not None:
                    ds.close()
                ds, current_path = xr.open_dataset(path), path
            fire_im.set_data(project(ds["C_weekly"].isel(time=k).values))
            date_txt.set_text(t.strftime("%b  %Y"))
            writer.grab_frame()
    ds.close()
    plt.close(fig)

    nframes, secs = imageio_ffmpeg.count_frames_and_secs(str(tmp))
    if nframes != len(index):
        sys.exit(f"Video has {nframes} frames, expected {len(index)}; left at {tmp}")
    shutil.move(str(tmp), str(FINAL_VIDEO))
    print(f"Done: {FINAL_VIDEO}  ({nframes} frames, {secs:.1f} s)")


if __name__ == "__main__":
    main()