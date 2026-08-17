"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: precipitation_platecarree.py
Description: Generates the precipitation climatology animation from ERA5 reanalysis data (1979-2000), rendered in the Plate Carrée projection, recentered on the Pacific. Weekly totals (cm water equivalent) computed server-side by the derived-era5-single-levels-daily-statistics CDS dataset, averaged into a 22-year climatological mean per 48-frame weekly cycle.
Note: For the Robinson, Foucaut, and Nicolosi projections, see the other scripts in the precipitation scripts folder.
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
from cartopy.util import add_cyclic_point
import imageio
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Volumes/CLIMATEDATA/precip/climatology/climatology_precip_48frame.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/precipitation_animation')
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ── Projection to render ──────────────────────────────────────────────────────
PROJ_NAME = 'platecarree'
CENTRAL_LON = 180  # centers the Pacific instead of the Atlantic (default is 0)
PROJ_CRS = ccrs.PlateCarree(central_longitude=CENTRAL_LON)

# ── Plot settings ─────────────────────────────────────────────────────────────
VMIN, VMAX = 0, 21
TICK_STEP = 3

# Smooth continuous gradient (not discrete flat blocks): white at 0 cm
# gradually blends into blue over just the first 0.3 cm (shrunk from 1 cm
# to cut down on how much of the colorbar reads as blank/white), then
# continues through the standard blue → green → yellow → orange → red →
# purple wet gradient up to 20 cm. Positions are fractions of VMAX.
color_stops = [
    (0 / VMAX, '#ffffff'),    # white           (0 cm)
    (0.3 / VMAX, '#d0e8f0'),  # pale blue       (0.3 cm)
    (2 / VMAX, '#7ab8e0'),    # light blue      (2 cm)
    (3 / VMAX, '#3060c0'),    # blue            (3 cm)
    (6 / VMAX, '#00994d'),    # green           (6 cm)
    (8 / VMAX, '#66cc00'),    # yellow-green    (8 cm)
    (9 / VMAX, '#ffff00'),    # yellow          (9 cm)
    (10 / VMAX, '#ff9900'),   # orange          (10 cm)
    (11 / VMAX, '#ff3300'),   # red-orange      (11 cm)
    (12 / VMAX, '#990000'),   # dark red        (12 cm)
    (21 / VMAX, '#660066'),   # purple          (21 cm)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_precip', color_stops, N=256)
BOUNDS = np.arange(VMIN, VMAX + 1, TICK_STEP)  # 0, 3, 6, ..., up to VMAX
if BOUNDS[-1] != VMAX:
    BOUNDS = np.append(BOUNDS, VMAX)  # always label the true top of the range

# ── Load pre-computed climatology ──────────────────────────────────────────────
ds_clim = xr.open_dataset(DATA_FILE)

# ERA5 stores longitude as 0-360; normalize to -180-180 and sort so Cartopy's
# PlateCarree doesn't draw a seam/duplicate wrap at the prime meridian.
ds_clim = ds_clim.assign_coords(
    longitude=(((ds_clim['longitude'] + 180) % 360) - 180)
).sortby('longitude')

lats = ds_clim['latitude'].values
lons = ds_clim['longitude'].values
labels = ds_clim['frame_label'].values

frames_data = []
for i in range(ds_clim.sizes['frame']):
    field = ds_clim['tp_weekly_total_cm'].isel(frame=i).values
    label = str(labels[i])
    frames_data.append((label, field))
    print(f'Loaded: {label}')

# ── Generate frames + MP4 ───────────────────────────────────────────────────────
print(f'=== Rendering projection: {PROJ_NAME} (centered on {CENTRAL_LON}) ===')
frames_dir = OUTPUT_DIR / PROJ_NAME / 'frames'
frames_dir.mkdir(parents=True, exist_ok=True)
frame_paths = []

for i, (label, field) in enumerate(frames_data):
    # The data's antimeridian seam (where lon wraps from -180 back to 180)
    # used to sit harmlessly at the map's left/right edge when centered on
    # the Atlantic. Recentered on the Pacific, that seam now falls at the
    # middle of the view -- add_cyclic_point closes the loop so there's no
    # visible gap/artifact right where the Pacific is supposed to be clean.
    field_cyclic, lons_cyclic = add_cyclic_point(field, coord=lons)

    fig, ax = plt.subplots(
        figsize=(12, 6),
        subplot_kw={'projection': PROJ_CRS}
    )

    im = ax.contourf(
        lons_cyclic, lats, field_cyclic,
        levels=np.linspace(VMIN, VMAX, 61),
        cmap=CMAP,
        extend='neither',
        transform=ccrs.PlateCarree()  # data stays in standard (unrotated) lon/lat
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.8, edgecolor='black')
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')
    ax.set_global()

    # ── Lat/lon axis labels (matching the ERA-5 | Climate Reanalyzer style:
    # 90°N/45°N/0°/45°S/90°S and 180°/90°W/0°/90°E/180°) ────────────────────
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
    gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    # ── Colorbar (rectangular, no pointed extend arrows; ticks 0-13) ────────
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal',
                         pad=0.08, fraction=0.05, aspect=20, extend='neither')
    cbar.set_label('Precipitation (cm)', fontsize=11)
    cbar.set_ticks(BOUNDS)
    cbar.set_ticklabels([str(t) for t in BOUNDS])
    cbar.ax.tick_params(labelsize=9)

    # ── Titles ───────────────────────────────────────────────────────────────
    ax.text(0.0, 1.02, 'ERA-5 | Climate Reanalyzer',
             transform=ax.transAxes, fontsize=10,
             fontweight='bold', ha='left', va='bottom')
    ax.text(1.0, 1.02, f'{label}; 1979–2000 Mean Weekly Total',
             transform=ax.transAxes, fontsize=10,
             fontweight='bold', ha='right', va='bottom')

    fpath = frames_dir / f'frame_{i:03d}.png'
    plt.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close()
    frame_paths.append(str(fpath))
    print(f'Saved frame {i+1}/{len(frames_data)}')

# ── Assemble into MP4 ─────────────────────────────────────────────────────────
with imageio.get_writer(OUTPUT_DIR / 'precipitation_climatology.mp4',
                         fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print('Done! Animation saved to:', OUTPUT_DIR / 'precipitation_climatology.mp4')