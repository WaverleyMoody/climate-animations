import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio
from pathlib import Path
import matplotlib.colors as mcolors

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Users/waverleymoody/Downloads/climate_data_by_variable/2m_temp_data.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/2m_temp_animation')
FRAMES_DIR = OUTPUT_DIR / 'frames'
OUTPUT_DIR.mkdir(exist_ok=True)
FRAMES_DIR.mkdir(exist_ok=True)

# ── Plot settings ─────────────────────────────────────────────────────────────
VMIN, VMAX = -60, 50

# Each tuple is (position between 0-1, color)
# 0.0 = -60°C, 1.0 = 50°C
# 20°C is at position (20 - (-60)) / (50 - (-60)) = 80/110 = 0.727
colors = [
    '#2d004b',  # dark purple   (-60°C)
    '#800080',  # purple
    '#8b008b',  # magenta
    '#ff00ff',  # bright pink
    '#0000cd',  # blue
    '#6699cc',  # gray blue
    '#87ceeb',  # sky blue
    '#00bfff',  # cyan
    '#39ff14',  # neon green
    '#00cc00',  # green
    '#006400',  # dark green
    '#ffff00',  # yellow
    '#ff6600',  # orange
    '#ff0000',  # red
    '#8b0000',  # dark red
    '#ffffff',  # white         (+50°C)
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_temp', colors, N=60)

# ── Extract the 1st, 8th, 15th, 22nd of each month and average across years ──
target_days = [1, 8, 15, 22]
MONTHS = list(range(1, 13))
MONTH_NAMES = ['Jan','Feb','Mar','Apr','May','Jun',
               'Jul','Aug','Sep','Oct','Nov','Dec']

frames_data = []

# Load the single consolidated file — it already has a proper time coordinate.
# chunks= keeps this dask-backed/lazy so we don't pull the whole 22-year global
# grid into memory at once (without it, open_dataset loads eagerly and can crash).
ds_all = xr.open_dataset(DATA_FILE, chunks={'time': 50})
t2m_all = ds_all['t2m'] - 273.15

lats = ds_all['latitude'].values
lons = ds_all['longitude'].values

for month in MONTHS:
    for day in target_days:
        subset = t2m_all.sel(time=((t2m_all['time'].dt.month == month) &
                                   (t2m_all['time'].dt.day == day)))
        mean_field = subset.mean(dim='time').compute().values
        label = f'{MONTH_NAMES[month - 1]} {day}'
        frames_data.append((label, mean_field, lats, lons))
        print(f'Processed: {label}')

frame_paths = []

# ── Generate one frame per date ───────────────────────────────────────────────
for i, (label, field, lats, lons) in enumerate(frames_data):
    fig, ax = plt.subplots(
        figsize=(12, 6),
        subplot_kw={'projection': ccrs.PlateCarree()}
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
    cbar.set_label('Temperature at 2 m (°C)', fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_ticks(np.arange(-60, 51, 10))
    cbar.set_ticklabels([str(t) for t in np.arange(-60, 51, 10)])

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
with imageio.get_writer(OUTPUT_DIR / '2m_temperature_climatology.mp4',
                        fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print('Done! Animation saved to:', OUTPUT_DIR / '2m_temperature_climatology.mp4')