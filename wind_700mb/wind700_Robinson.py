"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: wind700_Robinson.py
Description: Generates the 700mb wind speed climatology animation from
ERA5 reanalysis data (1979-2000), with mean wind vectors overlaid,
rendered in the Robinson projection.

Note: For the Plate Carrée, Foucaut, and Nicolosi projections, see the
other scripts in the wind_700mb scripts folder.
"""

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio.v2 as imageio
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Volumes/CLIMATEDATA/wind700/climatology/climatology_wind700_48frame.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/wind700_animation')
OUTPUT_DIR.mkdir(exist_ok=True)

# ── Projection to render ──────────────────────────────────────────────────────
PROJ_NAME = 'robinson'
PROJ_CRS  = ccrs.Robinson()

# ── Plot settings ─────────────────────────────────────────────────────────────
VMIN, VMAX = 0, 30

colors = [
    (0.00, '#ffffff'),  # white         (0 m/s)
    (0.20, '#d0e8f0'),  # light blue    (6 m/s)
    (0.35, '#469a46'),  # light green   (10.5 m/s)
    (0.44, '#ffff00'),  # yellow        (13 m/s)
    (0.50, '#ff6600'),  # orange        (15 m/s)
    (0.60, '#ff3300'),  # darker orange (18 m/s)
    (0.70, '#ff0000'),  # light red     (21 m/s)
    (0.80, '#cc0000'),  # red           (24 m/s)
    (0.90, '#990000'),  # dark red      (27 m/s)
    (1.00, '#660000'),  # darkest red   (30 m/s)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_wind700', colors, N=35)

# ── Load pre-computed climatology ─────────────────────────────────────────────
ds_clim = xr.open_dataset(DATA_FILE)
print("Variables:", list(ds_clim.data_vars))
print("Dims:", dict(ds_clim.sizes))

lats = ds_clim['latitude'].values
lons = ds_clim['longitude'].values

# No frame_label in this dataset -- generate week labels from frame index.
# 48 frames = 12 months x 4 weeks (starting on 1st, 8th, 15th, 22nd).
month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
week_starts = [1, 8, 15, 22]
labels = [f'{month_names[m]} {d}'
          for m in range(12) for d in week_starts]

frames_data = []
for i in range(ds_clim.sizes['frame']):
    mean_speed = ds_clim['wind_speed_700mb'].isel(frame=i).squeeze().values
    frame_u    = ds_clim['u700'].isel(frame=i).squeeze().values
    frame_v    = ds_clim['v700'].isel(frame=i).squeeze().values
    label      = labels[i]
    frames_data.append((label, mean_speed, frame_u, frame_v))
    print(f'Loaded: {label}')

# ── Generate frames + MP4 ─────────────────────────────────────────────────────
print(f'=== Rendering projection: {PROJ_NAME} ===')
frames_dir = OUTPUT_DIR / PROJ_NAME / 'frames'
frames_dir.mkdir(parents=True, exist_ok=True)
frame_paths = []

for i, (label, mean_speed, frame_u, frame_v) in enumerate(frames_data):
    fig, ax = plt.subplots(
        figsize=(12, 6),
        subplot_kw={'projection': PROJ_CRS}
    )

    im = ax.pcolormesh(
        lons, lats, mean_speed,
        vmin=VMIN, vmax=VMAX,
        cmap=CMAP,
        transform=ccrs.PlateCarree(),
        shading='auto'
    )

    # Subsample for quiver arrows.
    # Convert longitudes from 0–360 to -180/180 if needed.
    step = 30
    lons_sub = lons[::step]
    lats_sub = lats[::step]
    u_sub    = frame_u[::step, ::step]
    v_sub    = frame_v[::step, ::step]

    lons_sub_180 = np.where(lons_sub > 180, lons_sub - 360, lons_sub)
    sort_idx     = np.argsort(lons_sub_180)
    lons_sub_180 = lons_sub_180[sort_idx]
    u_sub        = u_sub[:, sort_idx]
    v_sub        = v_sub[:, sort_idx]

    lons_grid, lats_grid = np.meshgrid(lons_sub_180, lats_sub)

    ax.quiver(
        lons_grid, lats_grid,
        u_sub, v_sub,
        transform=ccrs.PlateCarree(),
        color='black',
        scale=450,
        width=0.001,
        headwidth=3,
        headlength=3,
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='black')
    ax.add_feature(cfeature.BORDERS,   linewidth=0.8, edgecolor='black')
    ax.add_feature(cfeature.STATES,    linewidth=0.5, edgecolor='black')
    ax.set_global()

    gl = ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', alpha=0.5)
    gl.xlocator = matplotlib.ticker.FixedLocator([-180, -90, 0, 90, 180])
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -45, 0, 45, 90])

    # Longitude labels along the bottom edge (projection-aware)
    for lon, lon_label in [(-180, '180°'), (-90, '90°W'), (0, '0°'),
                            (90, '90°E'), (180, '180°')]:
        x, y = PROJ_CRS.transform_point(lon, -90, ccrs.PlateCarree())
        ax.annotate(lon_label, xy=(x, y), xycoords='data',
                    xytext=(0, -6), textcoords='offset points',
                    ha='center', va='top', fontsize=9, annotation_clip=False)

    # Latitude labels along the left edge (projection-aware)
    # Robinson needs extra clearance at the poles so the label doesn't
    # collide with the curved edge of the map.
    for lat, lat_label in [(-90, '90°S'), (-45, '45°S'), (0, '0°'),
                             (45, '45°N'), (90, '90°N')]:
        x, y = PROJ_CRS.transform_point(-180, lat, ccrs.PlateCarree())
        x_offset = -18 if abs(lat) == 90 else -8
        ax.annotate(lat_label, xy=(x, y), xycoords='data',
                    xytext=(x_offset, 0), textcoords='offset points',
                    ha='right', va='center', fontsize=9, annotation_clip=False)

    cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                        pad=0.06, fraction=0.12, shrink=0.5,
                        aspect=20, extend='neither')
    cbar.set_label('Wind Speed at 700 mb (m/s)', fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_ticks(np.arange(0, 31, 5))
    cbar.set_ticklabels([str(t) for t in np.arange(0, 31, 5)])

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

video_name = 'wind700_climatology_robinson.mp4'
with imageio.get_writer(OUTPUT_DIR / video_name, fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print(f'Done! Animation saved to: {OUTPUT_DIR / video_name}')