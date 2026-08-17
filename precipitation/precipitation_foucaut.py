"""
SDSU Climate Informatics Lab
San Diego State University
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
Python Code Version 1.0.0

A reproduction of the University of Washington General Circulation
Animations Library by Professor John Michael Wallace.

Script: precipitation_foucaut.py
Description: Generates the precipitation climatology animation from ERA5 reanalysis data (1979-2000), rendered in the Foucaut projection. Weekly totals (cm water equivalent) computed server-side by the derived-era5-single-levels-daily-statistics CDS dataset, averaged into a 22-year climatological mean per 48-frame weekly cycle.
Note: For the Plate Carrée, Robinson, and Nicolosi projections, see the other scripts in the precipitation scripts folder.

--- Implementation notes ---
Uses PROJ's `fouc` operation (plain Foucaut, NOT `fouc_s`/Foucaut
Sinusoidal — a different, parameterized blend family). `fouc` has a
working inverse transform, so unlike Nicolosi this can use a normal
Cartopy Projection subclass and standard GeoAxes rendering rather than a
manual reprojection pipeline.

Foucaut has pointed poles: every longitude at lat=90/-90 collapses to a
single point, so there's no clean row to place longitude labels along
without overlap — longitude labels are omitted, and each pole gets one
label instead of five overlapping ones.
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
DATA_FILE = Path('/Volumes/CLIMATEDATA/precip/climatology/climatology_precip_48frame.nc')
OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/precipitation_animation')
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)


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
        # hardcoded guess (the earlier version's hardcoded Robinson-style
        # bounds is what cut off the top/bottom of the map), so this stays
        # correct even if the true extent doesn't match some other
        # projection's typical bounds.
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


# ── Projection to render ──────────────────────────────────────────────────────
PROJ_NAME = 'foucaut'
PROJ_CRS = Foucaut()

# ── Plot settings ─────────────────────────────────────────────────────────────
VMIN, VMAX = 0, 21
TICK_STEP = 3

# Smooth continuous gradient (not discrete flat blocks): white at 0 cm
# gradually blends into blue over just the first 0.3 cm (shrunk from 1 cm
# to cut down on how much of the colorbar reads as blank/white), then
# continues through the standard blue → green → yellow → orange → red →
# purple wet gradient up to 21 cm. Positions are fractions of VMAX.
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

# ERA5 stores longitude as 0-360; normalize to -180-180 and sort to prevent
# antimeridian reprojection distortion.
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
print(f'=== Rendering projection: {PROJ_NAME} ===')
frames_dir = OUTPUT_DIR / PROJ_NAME / 'frames'
frames_dir.mkdir(parents=True, exist_ok=True)
frame_paths = []

for i, (label, field) in enumerate(frames_data):
    fig, ax = plt.subplots(
        figsize=(12, 6),
        subplot_kw={'projection': PROJ_CRS}
    )

    im = ax.contourf(
        lons, lats, field,
        levels=np.linspace(VMIN, VMAX, 61),
        cmap=CMAP,
        extend='neither',
        transform=ccrs.PlateCarree()
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.8, edgecolor='black')
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')
    ax.set_global()

    gl = ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', alpha=0.5)
    gl.xlocator = matplotlib.ticker.FixedLocator([-180, -90, 0, 90, 180])
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -45, 0, 45, 90])

    # ── Lat/lon labels (projection-aware placement) ──────────────────────────
    # Foucaut has pointed poles: every longitude at lat=90/-90 collapses to
    # the same single point, so there's no clean row to place longitude
    # labels along without them overlapping. Following the confirmed-working
    # reference implementation, longitude labels are omitted; each pole gets
    # a single unambiguous label instead of five overlapping ones.
    x, y = PROJ_CRS.transform_point(0, 90, ccrs.PlateCarree())
    ax.annotate('90°N', xy=(x, y), xycoords='data',
                xytext=(0, 6), textcoords='offset points',
                ha='center', va='bottom', fontsize=9, annotation_clip=False)
    x, y = PROJ_CRS.transform_point(0, -90, ccrs.PlateCarree())
    ax.annotate('90°S', xy=(x, y), xycoords='data',
                xytext=(0, -6), textcoords='offset points',
                ha='center', va='top', fontsize=9, annotation_clip=False)

    # Latitude labels along the left edge — 90°N/90°S are excluded here
    # since they're already covered by the dedicated pole labels above.
    for lat, lat_label in [(-45, '45°S'), (0, '0°'), (45, '45°N')]:
        x, y = PROJ_CRS.transform_point(-180, lat, ccrs.PlateCarree())
        ax.annotate(lat_label, xy=(x, y), xycoords='data',
                    xytext=(-8, 0), textcoords='offset points',
                    ha='right', va='center', fontsize=9, annotation_clip=False)

    # ── Colorbar (rectangular, no pointed extend arrows) ─────────────────────
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal',
                         pad=0.12, fraction=0.05, aspect=20, extend='neither')
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
with imageio.get_writer(OUTPUT_DIR / 'precipitation_climatology_foucaut.mp4',
                         fps=4, codec='libx264') as writer:
    for fp in frame_paths:
        writer.append_data(imageio.imread(fp))

print('Done! Animation saved to:', OUTPUT_DIR / 'precipitation_climatology_foucaut.mp4')