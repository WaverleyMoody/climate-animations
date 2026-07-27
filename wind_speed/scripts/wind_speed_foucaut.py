"""
animate_foucaut_wind.py

10m wind speed climatology animation in the Foucaut projection, with mean
wind vectors overlaid.

Same Foucaut projection class as animate_foucaut_2m_temp.py / 
animate_foucaut_SLP.py. Longitude is normalized to [-180, 180) and sorted
for the WHOLE dataset up front (fixes the pcolormesh antimeridian
distortion seen on the temp version) -- this makes the original script's
quiver-only local longitude re-sort redundant, so that's been removed
rather than left in as dead logic.

Since Foucaut has a working inverse transform (unlike Nicolosi), Cartopy's
ax.quiver(transform=ccrs.PlateCarree()) handles wind vector rotation
automatically -- no manual Jacobian rotation needed here, unlike the
Nicolosi wind script.

NOTE: still not execution-tested end to end (ongoing tool access issue).
The Foucaut class and longitude fix carry over from the temp/SLP versions,
which the user confirmed working; wind-specific changes (variables,
colormap, quiver call) mirror animate_wind.py but haven't been run.
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
import imageio.v2 as imageio
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Arial'

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_FILE = Path('/Users/waverleymoody/Downloads/climate_data_by_variable/wind_climatology.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/wind_animation')
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

# ── Plot settings (identical to animate_wind.py) ──────────────────────────────
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

# ── Load pre-computed climatology ──────────
ds_clim = xr.open_dataset(DATA_FILE)

# Normalize longitude to [-180, 180) and sort ascending -- fixes the
# pcolormesh antimeridian distortion for the whole dataset up front,
# which also makes the quiver subsampling's local re-sort unnecessary.
lon_vals = ds_clim['longitude'].values
lon_norm = ((lon_vals + 180) % 360) - 180
ds_clim = ds_clim.assign_coords(longitude=lon_norm).sortby('longitude')

lats = ds_clim['latitude'].values
lons = ds_clim['longitude'].values
labels = ds_clim['frame_label'].values

frames_data = []
for i in range(ds_clim.sizes['frame']):
    mean_speed = ds_clim['wind_speed'].isel(frame=i).values
    frame_u = ds_clim['u10'].isel(frame=i).values
    frame_v = ds_clim['v10'].isel(frame=i).values
    label = str(labels[i])
    frames_data.append((label, mean_speed, frame_u, frame_v))
    print(f'Loaded: {label}')

# ── Generate frames + MP4 ──────────────────────────────────────────────────────
frames_dir = OUTPUT_DIR / PROJ_NAME / 'frames'
frames_dir.mkdir(parents=True, exist_ok=True)
frame_paths = []

for i, (label, mean_speed, frame_u, frame_v) in enumerate(frames_data):
    fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': PROJ_CRS})

    im = ax.pcolormesh(
        lons, lats, mean_speed,
        vmin=VMIN, vmax=VMAX,
        cmap=CMAP,
        transform=ccrs.PlateCarree(),
        shading='auto'
    )

    # lons/lats are already normalized+sorted at the dataset level above,
    # so subsampling here needs no separate re-sort (unlike the original
    # script, where lons was still 0-360 at this point).
    step = 30
    lons_sub = lons[::step]
    lats_sub = lats[::step]
    u_sub = frame_u[::step, ::step]
    v_sub = frame_v[::step, ::step]
    lons_grid, lats_grid = np.meshgrid(lons_sub, lats_sub)

    # Foucaut has a working inverse transform (unlike Nicolosi), so
    # Cartopy handles wind vector rotation automatically here -- same as
    # it already does for the Robinson output.
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

    fpath = frames_dir / f'frame_{i:03d}.png'
    plt.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close()
    frame_paths.append(str(fpath))
    print(f'Saved frame {i+1}/{len(frames_data)}')

video_name = 'wind_climatology_foucaut.mp4'
with imageio.get_writer(OUTPUT_DIR / video_name, fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print(f'Done! Animation saved to: {OUTPUT_DIR / video_name}')