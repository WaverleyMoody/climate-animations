"""
SDSU Climate Informatics
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
San Diego State University

Script: add_uv_to_climatology.py
Description: Reads the 48 per-slot u/v cache files produced by
download_uv700_climatology.py and merges them into the existing
climatology_wind700_48frame.nc as variables u700 and v700.

Run once, after download_uv700_climatology.py has completed all 48 slots.
"""

import xarray as xr
import numpy as np
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
CLIMO_FILE = Path('/Volumes/CLIMATEDATA/wind700/climatology/climatology_wind700_48frame.nc')
TEMP_DIR   = Path('/Volumes/CLIMATEDATA/wind700/uv_temp')

WEEK_STARTS = [1, 8, 15, 22]
MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# ── Verify all 48 slot files exist before doing anything ──────────────────────
missing = []
for month in range(1, 13):
    for week_start in WEEK_STARTS:
        slot_file = TEMP_DIR / f'uv700_climo_{month:02d}_{week_start:02d}.nc'
        if not slot_file.exists():
            missing.append(slot_file.name)

if missing:
    print(f'ERROR: {len(missing)} slot file(s) missing from {TEMP_DIR}:')
    for m in missing:
        print(f'  {m}')
    print('Re-run download_uv700_climatology.py to fill the gaps, then retry.')
    raise SystemExit(1)

print(f'All 48 slot files present in {TEMP_DIR}')

# ── Load all 48 slots in frame order ─────────────────────────────────────────
u_climo = []
v_climo = []

for month in range(1, 13):
    for week_start in WEEK_STARTS:
        slot_file = TEMP_DIR / f'uv700_climo_{month:02d}_{week_start:02d}.nc'
        label = f'{MONTH_NAMES[month-1]} {week_start}'
        ds = xr.open_dataset(slot_file)
        u_climo.append(ds['u700'].values)
        v_climo.append(ds['v700'].values)
        lats = ds['latitude'].values
        lons = ds['longitude'].values
        ds.close()
        print(f'Loaded: {label}')

u_array = np.stack(u_climo, axis=0)  # (48, lat, lon)
v_array = np.stack(v_climo, axis=0)

print(f'\nu_array shape: {u_array.shape}')
print(f'v_array shape: {v_array.shape}')

# ── Merge into the existing climatology file ──────────────────────────────────
print(f'\nOpening existing climatology file...')
ds_existing = xr.open_dataset(CLIMO_FILE)
print('Existing variables:', list(ds_existing.data_vars))
print('Existing dims:', dict(ds_existing.sizes))

frame_coord = np.arange(48)
lat_coord   = ds_existing['latitude']
lon_coord   = ds_existing['longitude']

ds_uv = xr.Dataset({
    'u700': xr.DataArray(
        u_array,
        dims=['frame', 'latitude', 'longitude'],
        coords={'frame': frame_coord, 'latitude': lat_coord, 'longitude': lon_coord},
        attrs={'long_name': 'Climatological mean u-component of wind at 700 hPa',
               'units': 'm s**-1'}
    ),
    'v700': xr.DataArray(
        v_array,
        dims=['frame', 'latitude', 'longitude'],
        coords={'frame': frame_coord, 'latitude': lat_coord, 'longitude': lon_coord},
        attrs={'long_name': 'Climatological mean v-component of wind at 700 hPa',
               'units': 'm s**-1'}
    ),
})

ds_merged = xr.merge([ds_existing, ds_uv])
ds_existing.close()

# Write to a temp file then replace atomically, so the original is never
# left in a partially-written state if the write is interrupted.
print('\nWriting merged file...')
tmp_out = CLIMO_FILE.with_suffix('.tmp.nc')
ds_merged.to_netcdf(tmp_out)
ds_merged.close()
tmp_out.replace(CLIMO_FILE)

# Verify the result
ds_check = xr.open_dataset(CLIMO_FILE)
print(f'\nDone! {CLIMO_FILE.name} now contains:')
print('  Variables:', list(ds_check.data_vars))
print('  Dims:', dict(ds_check.sizes))
ds_check.close()