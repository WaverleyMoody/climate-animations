"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace

Script: lightning_platecarree.py
Description: Generates the lightning climatology animation from WWLLN/WGLC
data (2010-2025), rendered in the Plate Carrée projection. 365 calendar-day
frames (Feb 29 excluded), each the multi-year daily-climatological mean,
annualized (x365) to match the reference product's strokes km-2 yr-1 units.

Note: For the Robinson projection, see lightning_robinson.py.
For the Foucaut and Nicolosi projections, visit the lightning
scripts folder.
"""

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import imageio
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Volumes/CLIMATEDATA/lightning/climatology/climatology_lightning_365frame.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/lightning_animation')
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ── Projection to render ──────────────────────────────────────────────────────
PROJ_NAME = 'platecarree'
PROJ_CRS = ccrs.PlateCarree()

# ── Plot settings (log scale) ────────────────────────────────────────────────
VMIN, VMAX = 0.003, 30  # strokes km-2 yr-1, matching the reference legend
TICK_VALUES = [0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]

# Color stops: white (near-zero) -> Tiffany blue -> yellow -> orange -> red
# (highest), on a LOG scale (fractions of log10(VMAX/VMIN)) so spacing is
# even in log-space. White is now just a thin sliver near VMIN, and Tiffany
# blue is repeated at two anchor values (0.0045 and 0.03) to create a solid
# blue plateau spanning that range, rather than an instant transition —
# this pushes a lot more of the map into blue before yellow/orange/red kick in.
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

# ── Load pre-computed climatology ──────────────────────────────────────────────
ds_clim = xr.open_dataset(DATA_FILE)
lats = ds_clim['latitude'].values
lons = ds_clim['longitude'].values
labels = ds_clim['frame_label'].values

frames_data = []
for i in range(ds_clim.sizes['frame']):
    field = ds_clim['lightning_density_annual'].isel(frame=i).values
    # Log scale can't render exact zero (common over oceans/poles) — clip
    # up to VMIN so those cells render as the lowest color instead of
    # erroring or dropping out.
    field = np.clip(field, VMIN, None)
    label = str(labels[i])
    frames_data.append((label, field))
    print(f'Loaded: {label}')

# ── Generate frames + MP4 ───────────────────────────────────────────────────────
print(f'=== Rendering projection: {PROJ_NAME} ===')
frames_dir = OUTPUT_DIR / PROJ_NAME / 'frames'
frames_dir.mkdir(parents=True, exist_ok=True)
frame_paths = []

log_levels = np.logspace(np.log10(VMIN), np.log10(VMAX), 121)

for i, (label, field) in enumerate(frames_data):
    fig, ax = plt.subplots(
        figsize=(12, 6),
        subplot_kw={'projection': PROJ_CRS}
    )

    im = ax.contourf(
        lons, lats, field,
        levels=log_levels,
        cmap=CMAP,
        norm=NORM,
        extend='max',
        transform=ccrs.PlateCarree()
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.8, edgecolor='black')
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')
    ax.set_global()

    # ── Lat/lon axis labels (matching the ERA-5 | Climate Reanalyzer style) ──
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
    gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    # ── Colorbar (log scale, ticks matching the reference legend) ────────────
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal',
                         pad=0.08, fraction=0.05, aspect=20, extend='max')
    cbar.set_label('Lightning stroke density (strokes km$^{-2}$ yr$^{-1}$)', fontsize=11)
    cbar.set_ticks(TICK_VALUES)
    cbar.set_ticklabels([str(t) for t in TICK_VALUES])
    cbar.ax.tick_params(labelsize=9)

    # ── Titles ───────────────────────────────────────────────────────────────
    ax.text(0.0, 1.02, 'WWLLN | WGLC',
             transform=ax.transAxes, fontsize=10,
             fontweight='bold', ha='left', va='bottom')
    ax.text(1.0, 1.02, f'{label}; 2010–2025 Climatology',
             transform=ax.transAxes, fontsize=10,
             fontweight='bold', ha='right', va='bottom')

    fpath = frames_dir / f'frame_{i:03d}.png'
    plt.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close()
    frame_paths.append(str(fpath))
    print(f'Saved frame {i+1}/{len(frames_data)}')

# ── Assemble into MP4 ─────────────────────────────────────────────────────────
with imageio.get_writer(OUTPUT_DIR / 'lightning_climatology.mp4',
                         fps=8, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print('Done! Animation saved to:', OUTPUT_DIR / 'lightning_climatology.mp4')
