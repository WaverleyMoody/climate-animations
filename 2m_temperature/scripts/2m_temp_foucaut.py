"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: 2m_temp_foucaut.py
Description: Generates the 2m temperature climatology animation from ERA5 reanalysis data (1979-2000), rendered in the Foucaut projection.
Note: For the Plate Carrée, Robinson, and Nicolosi projections, see the other scripts in the 2m_temperature scripts folder.
"""

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shapely.geometry as sgeom
import imageio
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Users/waverleymoody/Downloads/climate_data_by_variable/2m_temp_climatology.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/2m_temp_animation')
OUTPUT_DIR.mkdir(exist_ok=True)


# ── Custom Foucaut projection ──────────────────────────────────────────────────
class Foucaut(ccrs.Projection):
    """PROJ's plain `fouc` operation (NOT `fouc_s`/Foucaut Sinusoidal --
    that's a different, parameterized projection family). Confirmed to
    produce the correct pointed-pole Foucaut shape with no tuning needed.
    """

    def __init__(self):
        super().__init__(proj4_params=[('proj', 'fouc')])

        # Trace the outer boundary: up the antimeridian at lon=-180, then
        # back down at lon=180, forming a closed loop. Vectorized via
        # transform_points rather than a per-point loop.
        lats = np.linspace(-90, 90, 200)
        lons_left = np.full_like(lats, -180)
        lons_right = np.full_like(lats, 180)
        all_lons = np.concatenate([lons_left, lons_right[::-1]])
        all_lats = np.concatenate([lats, lats[::-1]])

        xyz = self.transform_points(ccrs.Geodetic(), all_lons, all_lats)
        x, y = xyz[:, 0], xyz[:, 1]
        valid = np.isfinite(x) & np.isfinite(y)
        x, y = x[valid], y[valid]

        self._boundary = sgeom.LinearRing(zip(x, y))
        # Computed from the actual traced boundary rather than a generic
        # hardcoded guess, so this stays correct even if the true extent
        # doesn't match some other projection's typical bounds.
        self._x_limits = (x.min(), x.max())
        self._y_limits = (y.min(), y.max())

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits

    @property
    def threshold(self):
        return 1e4


PROJ_CRS = Foucaut()
PROJ_NAME = 'foucaut'

# ── Plot settings (identical to animate_2m_temp.py) ───────────────────────────
VMIN, VMAX = -60, 50
colors = [
    '#2d004b', '#800080', '#8b008b', '#ff00ff', '#0000cd', '#6699cc',
    '#87ceeb', '#00bfff', '#39ff14', '#00cc00', '#006400', '#ffff00',
    '#ff6600', '#ff0000', '#8b0000', '#ffffff',
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_temp', colors, N=60)

# ── Load pre-computed climatology ──────────
ds_clim = xr.open_dataset(DATA_FILE)

# Normalize longitude to [-180, 180) and sort ascending -- see module
# docstring for why this matters here.
lon_vals = ds_clim['longitude'].values
lon_norm = ((lon_vals + 180) % 360) - 180
ds_clim = ds_clim.assign_coords(longitude=lon_norm).sortby('longitude')

lats = ds_clim['latitude'].values
lons = ds_clim['longitude'].values
labels = ds_clim['frame_label'].values

frames_data = []
for i in range(ds_clim.sizes['frame']):
    field = (ds_clim['t2m'].isel(frame=i) - 273.15).values
    label = str(labels[i])
    frames_data.append((label, field))
    print(f'Loaded: {label}')

# ── Generate frames + MP4 ──────────────────────────────────────────────────────
frames_dir = OUTPUT_DIR / PROJ_NAME / 'frames'
frames_dir.mkdir(parents=True, exist_ok=True)
frame_paths = []

for i, (label, field) in enumerate(frames_data):
    fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': PROJ_CRS})

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

    # Single pole labels, since the whole pole is one point here (no
    # longitude labels -- removed per request).
    x, y = PROJ_CRS.transform_point(0, 90, ccrs.PlateCarree())
    ax.annotate('90°N', xy=(x, y), xycoords='data',
                xytext=(0, 6), textcoords='offset points',
                ha='center', va='bottom', fontsize=9, annotation_clip=False)
    x, y = PROJ_CRS.transform_point(0, -90, ccrs.PlateCarree())
    ax.annotate('90°S', xy=(x, y), xycoords='data',
                xytext=(0, -6), textcoords='offset points',
                ha='center', va='top', fontsize=9, annotation_clip=False)

    # Latitude labels along the left edge, excluding the poles (handled above).
    for lat, lat_label in [(-45, '45°S'), (0, '0°'), (45, '45°N')]:
        x, y = PROJ_CRS.transform_point(-180, lat, ccrs.PlateCarree())
        ax.annotate(lat_label, xy=(x, y), xycoords='data',
                    xytext=(-8, 0), textcoords='offset points',
                    ha='right', va='center', fontsize=9, annotation_clip=False)

    cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                         pad=0.08, fraction=0.12, shrink=0.5,
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

    fpath = frames_dir / f'frame_{i:03d}.png'
    plt.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close()
    frame_paths.append(str(fpath))
    print(f'Saved frame {i+1}/{len(frames_data)}')

video_name = '2m_temperature_climatology_foucaut.mp4'
with imageio.get_writer(OUTPUT_DIR / video_name, fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print(f'Done! Animation saved to: {OUTPUT_DIR / video_name}')