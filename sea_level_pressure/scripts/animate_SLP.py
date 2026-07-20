"""
SDSU Climate Informatics
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
San Diego State University

A reproduction of the University of Washington General Circulation
Animations Library, originally created by Professor John Michael Wallace.

Script: animate_SLP.py
Description: Generates the sea level pressure climatology animation from
ERA5 reanalysis data (1979-2000), rendered across multiple map projections.
"""

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Users/waverleymoody/Downloads/climate_data_by_variable/slp_climatology.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/slp_animation')
OUTPUT_DIR.mkdir(exist_ok=True)


# ── Projections to render ─────────────────────────────────────────────────────
PROJECTIONS = {
    'platecarree': ccrs.PlateCarree(),
    'robinson': ccrs.Robinson(),
}

# ── Plot settings ─────────────────────────────────────────────────────────────
VMIN, VMAX = 970, 1040

colors = [
    '#2d004b',  # dark purple  (970 mb)
    '#800080',  # purple
    '#8b008b',  # magenta
    '#ff00ff',  # bright pink
    '#003153',  # dark blue
    '#b0c4de',  # blue
    '#ffffff',  # white
    '#ff6600',  # orange
    '#ff0000',  # red
    '#8b0000',  # dark red     (1040 mb)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_slp', colors, N=60)

# ── Load pre-computed climatology ──────────
ds_clim = xr.open_dataset(DATA_FILE)
lats = ds_clim['latitude'].values
lons = ds_clim['longitude'].values
labels = ds_clim['frame_label'].values

frames_data = []
for i in range(ds_clim.sizes['frame']):
    field = (ds_clim['msl'].isel(frame=i) / 100).values  # Pa → hPa/mb
    label = str(labels[i])
    frames_data.append((label, field))
    print(f'Loaded: {label}')

frame_paths = []

# ── Generate one frame per date ───────────────────────────────────────────────
for proj_name, proj_crs in PROJECTIONS.items():
    print(f'=== Rendering projection: {proj_name} ===')
    frames_dir = OUTPUT_DIR / proj_name / 'frames'
    frames_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []

    for i, (label, field) in enumerate(frames_data):
        fig, ax = plt.subplots(
            figsize=(12, 6),
            subplot_kw={'projection': proj_crs}
        )

        im = ax.pcolormesh(
            lons, lats, field,
            vmin=VMIN, vmax=VMAX,
            cmap=CMAP,
            transform=ccrs.PlateCarree(),
            shading='auto'
        )

        ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='black')
        ax.add_feature(cfeature.BORDERS, linewidth=0.8, edgecolor='black')
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')
        ax.set_global()

        gl = ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', alpha=0.5)
        gl.xlocator = matplotlib.ticker.FixedLocator([-180, -90, 0, 90, 180])
        gl.ylocator = matplotlib.ticker.FixedLocator([-90, -45, 0, 45, 90])

        # Add longitude labels on bottom (axes-fraction placement — a data
        # Add latitude labels on left
        for lon, lon_label in [(-180, '180°'), (-90, '90°W'), (0, '0°'), (90, '90°E'), (180, '180°')]:
            x, y = proj_crs.transform_point(lon, -90, ccrs.PlateCarree())
            ax.annotate(lon_label, xy=(x, y), xycoords='data',
                        xytext=(0, -6), textcoords='offset points',
                        ha='center', va='top', fontsize=9, annotation_clip=False)

        # Add latitude labels on left
        for lat, lat_label in [(-90, '90°S'), (-45, '45°S'), (0, '0°'), (45, '45°N'), (90, '90°N')]:
            x, y = proj_crs.transform_point(-180, lat, ccrs.PlateCarree())
            if abs(lat) == 90 and proj_name == 'robinson':
                x_offset = -18
            else:
                x_offset = -8
            ax.annotate(lat_label, xy=(x, y), xycoords='data',
                        xytext=(x_offset, 0), textcoords='offset points',
                        ha='right', va='center', fontsize=9, annotation_clip=False)

        cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                            pad=0.06, fraction=0.12, shrink=0.5,
                            aspect=20, extend='neither')
        cbar.set_label('Pressure at Mean Sea Level (mb)', fontsize=11)
        cbar.ax.tick_params(labelsize=9)
        cbar.set_ticks(np.arange(970, 1041, 10))
        cbar.set_ticklabels([str(t) for t in np.arange(970, 1041, 10)])

        for x in np.linspace(0, 1, 61)[1:-1]:
            cbar.ax.axvline(x, color='black', linewidth=0.5)

        ax.text(0.0, 1.02, 'ERA-5 | Climate Reanalyzer',
                transform=ax.transAxes, fontsize=10,
                fontweight='bold', ha='left', va='bottom')
        ax.text(1.0, 1.02, f'{label}; 1979–2000 Weekly Mean',
                transform=ax.transAxes, fontsize=10,
                fontweight='bold', ha='right', va='bottom')

        fpath = frames_dir / f'frame_{i:03d}.png'
        plt.savefig(fpath, dpi=150, bbox_inches='tight')
        plt.close()
        frame_paths.append(str(fpath))
        print(f'Saved frame {i+1}/{len(frames_data)}')

    video_name = ('slp_climatology.mp4' if proj_name == 'platecarree'
                  else f'slp_climatology_{proj_name}.mp4')
    with imageio.get_writer(OUTPUT_DIR / video_name, fps=4, codec='libx264') as writer:
        for fp in frame_paths:
            writer.append_data(imageio.imread(fp))

    print(f'Done! Animation saved to: {OUTPUT_DIR / video_name}')