import xarray as xr
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio
from pathlib import Path
import re
import pandas as pd

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR = Path('/Users/waverleymoody/Downloads/climate_daily_means')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/wind_animation')
FRAMES_DIR = OUTPUT_DIR / 'frames'
OUTPUT_DIR.mkdir(exist_ok=True)
FRAMES_DIR.mkdir(exist_ok=True)

# ── Plot settings ─────────────────────────────────────────────────────────────
VMIN, VMAX = 0, 16

colors = [
    (0.00,  '#ffffff'),  # white        (0 m/s)
    (0.20,  '#d0e8f0'),  # light blue   (3 m/s)
    (0.35,  "#469a46"),  # light green        (5.5 m/s)
    (0.44,  '#ffff00'),  # yellow       (7 m/s)
    (0.50,  '#ff6600'),  # orange       (8 m/s)
    (0.60,  '#ff3300'),  # darker orange(9.5 m/s)
    (0.70,  '#ff0000'),  # light red    (11 m/s)
    (0.80,  '#cc0000'),  # red          (13 m/s)
    (0.90,  '#990000'),  # dark red     (14.5 m/s)
    (1.00,  '#660000'),  # darkest red  (16 m/s)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_wind', colors, N=35)

# Load all pipeline files 
target_days = [1, 8, 15, 22]
YEARS = list(range(1979, 2001))
MONTHS = list(range(1, 13))
MONTH_NAMES = ['Jan','Feb','Mar','Apr','May','Jun',
               'Jul','Aug','Sep','Oct','Nov','Dec']

all_files = sorted(DATA_DIR.glob('daily_mean_*.nc'))
ds_all = xr.open_mfdataset(all_files, combine='nested', concat_dim='time',
                            drop_variables=['tp', 'step', 'surface', 'number'],
                            coords='minimal',
                            compat='override')

lats = ds_all['latitude'].values
lons = ds_all['longitude'].values

# Parse dates from filenames
dates = []
for f in all_files:
    match = re.search(r'daily_mean_(\d{4})_(\d{2})_(\d{2})', f.name)
    year, month, day = int(match.group(1)), int(match.group(2)), int(match.group(3))
    dates.append((year, month, day))

time_index = pd.DatetimeIndex([pd.Timestamp(y, m, d) for y, m, d in dates])
ds_all = ds_all.assign_coords(time=time_index)

u10_all = ds_all['u10']
v10_all = ds_all['v10']
speed_all = np.sqrt(u10_all**2 + v10_all**2)

static_u = u10_all.mean(dim='time').values
static_v = v10_all.mean(dim='time').values

frames_data = []

for month in MONTHS:
    for day in target_days:
        mask = ((speed_all['time'].dt.month == month) &
                (speed_all['time'].dt.day == day))
        mean_speed = speed_all.sel(time=mask).mean(dim='time').values
        label = f'{MONTH_NAMES[month - 1]} {day}'
        frames_data.append((label, mean_speed))
        print(f'Processed: {label}')

frame_paths = []

# ── Generate one frame per date ───────────────────────────────────────────────
for i, (label, mean_speed) in enumerate(frames_data):
    fig, ax = plt.subplots(
        figsize=(12, 6),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )

    
    im = ax.pcolormesh(
        lons, lats, mean_speed,
        vmin=VMIN, vmax=VMAX,
        cmap=CMAP,
        transform=ccrs.PlateCarree(),
        shading='auto'
    )

    # Subsample for streamlines to reduce memory usage
    # Convert longitudes from 0-360 to -180-180
    step = 30
    lons_sub = lons[::step]
    lats_sub = lats[::step]
    u_sub = static_u[::step, ::step]
    v_sub = static_v[::step, ::step]

    u_sub = static_u[::step, ::step]
    v_sub = static_v[::step, ::step]
    lons_sub_180 = np.where(lons_sub > 180, lons_sub - 360, lons_sub)
    sort_idx = np.argsort(lons_sub_180)
    lons_sub_180 = lons_sub_180[sort_idx]
    u_sub = u_sub[:, sort_idx]
    v_sub = v_sub[:, sort_idx]

    lons_grid, lats_grid = np.meshgrid(lons_sub_180, lats_sub)

    ax.quiver(
        lons_grid, lats_grid,
        u_sub, v_sub,
        transform=ccrs.PlateCarree(),
        color='black',
        scale=300,
        width=0.001,
        headwidth=3,
        headlength=3,
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.8, edgecolor='black')
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')
    ax.set_global()

    gl = ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', alpha=0.5)
    gl.xlocator = matplotlib.ticker.FixedLocator([-180, -90, 0, 90, 180])
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -45, 0, 45, 90])

    # Add longitude labels on bottom
    for lon, lon_label in [(-180, '180°'), (-90, '90°W'), (0, '0°'), (90, '90°E'), (180, '180°')]:
        ax.text(lon, -95, lon_label, transform=ccrs.PlateCarree(),
                ha='center', va='top', fontsize=9)

    # Add latitude labels on left
    for lat, lat_label in [(-90, '90°S'), (-45, '45°S'), (0, '0°'), (45, '45°N'), (90, '90°N')]:
        ax.text(-0.01, (lat + 90) / 180, lat_label,
                transform=ax.transAxes,
                ha='right', va='center', fontsize=9)

    cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                        pad=0.06, fraction=0.12, shrink=0.5,
                        aspect=20, extend='neither')
    cbar.set_label('Wind Speed at 10 m (m/s)', fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_ticks(np.arange(0, 17, 2))
    cbar.set_ticklabels([str(t) for t in np.arange(0, 17, 2)])


    ax.text(0.0, 1.02, 'ERA-5 | Climate Reanalyzer',
            transform=ax.transAxes, fontsize=10,
            fontweight='bold', ha='left', va='bottom')
    ax.text(1.0, 1.02, f'{label}; 1979–2000 Weekly Mean',
            transform=ax.transAxes, fontsize=10,
            fontweight='bold', ha='right', va='bottom')

    fpath = FRAMES_DIR / f'frame_{i:03d}.png'
    plt.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close()
    frame_paths.append(str(fpath))
    print(f'Saved frame {i+1}/{len(frames_data)}')

# ── Assemble into MP4 ─────────────────────────────────────────────────────────
with imageio.get_writer(OUTPUT_DIR / 'wind_climatology.mp4',
                        fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print('Done! Animation saved to:', OUTPUT_DIR / 'wind_climatology.mp4')