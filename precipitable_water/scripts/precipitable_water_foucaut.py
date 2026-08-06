"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: precipitable_water_foucaut.py
Description: Generates the precipitable water (total column water vapour) climatology animation from ERA5 reanalysis data (1979-2000), rendered in the Foucaut projection.
Note: For the Plate Carrée, Robinson, and Nicolosi projections, see the other scripts in the precipitable_water scripts folder.
"""

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shapely.geometry as sgeom
import imageio
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Users/waverleymoody/Downloads/climate_data_by_variable/precip_water_climatology.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/precipitable_water_animation')
OUTPUT_DIR.mkdir(exist_ok=True)


# ── Custom Foucaut projection ──────────────────────────────────────────────────
class Foucaut(ccrs.Projection):
    """PROJ's plain `fouc` operation (NOT `fouc_s`/Foucaut Sinusoidal).
    Confirmed working via animate_foucaut_2m_temp.py."""

    def __init__(self):
        super().__init__(proj4_params=[('proj', 'fouc')])

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

# ── Plot settings (identical to animate_precipitable_water.py) ───────────────
VMIN, VMAX = 0, 80
colors = [
    '#d3d3d3',  # light gray      (0 kg/m2)
    '#808080',  # gray            (6 kg/m2)
    "#4b2922",  # brown           (12 kg/m2)
    '#d2b48c',  # tan/light brown (20 kg/m2)
    "#adcfe6",  # light blue      (24 kg/m2)
    "#30306a",  # blue            (30 kg/m2)
    "#004A00",  # dark green      (35 kg/m2)
    '#00cc00',  # green           (40 kg/m2)
    '#ffff00',  # yellow          (50 kg/m2)
    '#ff6600',  # orange          (56 kg/m2)
    '#ff0000',  # red             (60 kg/m2)
    "#690909",  # dark red
    "#67013e",  # dark pink
    '#ff00ff',  # magenta         (80 kg/m2)
    "#38014b",  # dark purple
]
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_pw', colors, N=40)

# ── Load pre-computed climatology ──────────
ds_clim = xr.open_dataset(DATA_FILE)

# Normalize longitude to [-180, 180) and sort ascending -- fixes the
# antimeridian/reprojection distortion, same as temp/SLP/wind.
lon_vals = ds_clim['longitude'].values
lon_norm = ((lon_vals + 180) % 360) - 180
ds_clim = ds_clim.assign_coords(longitude=lon_norm).sortby('longitude')

lats = ds_clim['latitude'].values
lons = ds_clim['longitude'].values
labels = ds_clim['frame_label'].values

frames_data = []
for i in range(ds_clim.sizes['frame']):
    field = ds_clim['tcwv'].isel(frame=i).values
    label = str(labels[i])
    frames_data.append((label, field))
    print(f'Loaded: {label}')

# ── Generate frames + MP4 ──────────────────────────────────────────────────────
frames_dir = OUTPUT_DIR / PROJ_NAME / 'frames'
frames_dir.mkdir(parents=True, exist_ok=True)
frame_paths = []

for i, (label, field) in enumerate(frames_data):
    fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': PROJ_CRS})

    levels = np.linspace(VMIN, VMAX, 61)
    im = ax.contourf(
        lons, lats, field,
        levels=levels,
        cmap=CMAP,
        transform=ccrs.PlateCarree(),
        extend='neither'
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.8, edgecolor='black')
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')
    ax.set_global()

    gl = ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', alpha=0.5)
    gl.xlocator = matplotlib.ticker.FixedLocator([-180, -90, 0, 90, 180])
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -45, 0, 45, 90])

    # No longitude labels (removed per preference). Single pole labels,
    # since the whole pole is one point in this projection.
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
                        pad=0.06, fraction=0.12, shrink=0.5,
                        aspect=20)
    cbar.set_label('Precipitable Water (kg m-2)', fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_ticks(np.arange(0, 81, 10))
    cbar.set_ticklabels([str(t) for t in np.arange(0, 81, 10)])

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

video_name = 'precipitable_water_climatology_foucaut.mp4'
with imageio.get_writer(OUTPUT_DIR / video_name, fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print(f'Done! Animation saved to: {OUTPUT_DIR / video_name}')