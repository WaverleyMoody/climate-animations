"""
SDSU Climate Informatics
by Waverley Moody
Supervised by Distinguished Professor Samuel Shen
San Diego State University

Script: download_uv700_climatology.py
Description: Downloads ERA5 700mb u and v wind components for the 48
climatological week slots (1979-2000), computes the climatological mean
u and v per slot, and adds them to the existing climatology_wind700_48frame.nc
as variables u700 and v700.

This is a targeted re-download -- far smaller than the original pipeline
since it only fetches the specific days needed for each of the 48 week
slots rather than the full hourly archive. 48 CDS requests total.

Run once. Safe to re-run: already-completed slots are skipped.
"""

import cdsapi
import xarray as xr
import numpy as np
from pathlib import Path
import calendar

# ── Paths ─────────────────────────────────────────────────────────────────────
CLIMO_FILE = Path('/Volumes/CLIMATEDATA/wind700/climatology/climatology_wind700_48frame.nc')
TEMP_DIR   = Path('/Volumes/CLIMATEDATA/wind700/uv_temp')
TEMP_DIR.mkdir(parents=True, exist_ok=True)

YEARS      = list(range(1979, 2001))
WEEK_STARTS = [1, 8, 15, 22]
MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

client = cdsapi.Client()

def week_days(year, month, start_day):
    """Return list of day strings for the week starting on start_day."""
    last_day = calendar.monthrange(year, month)[1]
    if start_day == 22:
        end_day = last_day
    else:
        end_day = start_day + 6
    return [str(d) for d in range(start_day, min(end_day, last_day) + 1)]

# ── Download and compute climatological mean u/v per week slot ────────────────
u_climo = []
v_climo = []

frame_idx = 0
for month in range(1, 13):
    for week_start in WEEK_STARTS:
        label = f'{MONTH_NAMES[month-1]}_{week_start:02d}'
        out_nc = TEMP_DIR / f'uv700_climo_{month:02d}_{week_start:02d}.nc'

        if out_nc.exists():
            print(f'[{frame_idx+1}/48] {label}: loading cached result...')
            ds = xr.open_dataset(out_nc)
            u_climo.append(ds['u700'].values)
            v_climo.append(ds['v700'].values)
            ds.close()
            frame_idx += 1
            continue

        print(f'[{frame_idx+1}/48] {label}: downloading from CDS...')

        # Collect all days for this week slot across all years.
        # Use year 2000 (a leap year) to get day lists safely for Feb.
        days_for_slot = week_days(2000, month, week_start)

        grib_path = TEMP_DIR / f'_temp_{month:02d}_{week_start:02d}.grib'

        client.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'variable': [
                    'u_component_of_wind',
                    'v_component_of_wind',
                ],
                'pressure_level': '700',
                'year':  [str(y) for y in YEARS],
                'month': [f'{month:02d}'],
                'day':   days_for_slot,
                'time':  ['00:00', '06:00', '12:00', '18:00'],
                'data_format': 'grib',
                'download_format': 'unarchived',
            },
            str(grib_path)
        )

        # Load u and v from grib, average over all timesteps.
        import cfgrib
        ds_u = xr.open_dataset(str(grib_path), engine='cfgrib',
                               backend_kwargs={'filter_by_keys': {'shortName': 'u'}})
        ds_v = xr.open_dataset(str(grib_path), engine='cfgrib',
                               backend_kwargs={'filter_by_keys': {'shortName': 'v'}})

        u_mean = ds_u['u'].mean(dim='time').values  # (lat, lon)
        v_mean = ds_v['v'].mean(dim='time').values

        lats = ds_u['latitude'].values
        lons = ds_u['longitude'].values

        ds_u.close()
        ds_v.close()
        grib_path.unlink(missing_ok=True)

        # Cache this slot's result so we can resume if interrupted.
        ds_slot = xr.Dataset({
            'u700': xr.DataArray(u_mean, dims=['latitude', 'longitude'],
                                 coords={'latitude': lats, 'longitude': lons}),
            'v700': xr.DataArray(v_mean, dims=['latitude', 'longitude'],
                                 coords={'latitude': lats, 'longitude': lons}),
        })
        ds_slot.to_netcdf(out_nc)
        ds_slot.close()

        u_climo.append(u_mean)
        v_climo.append(v_mean)
        print(f'  Done: {label}')
        frame_idx += 1

# ── Merge u700 and v700 into the existing climatology file ────────────────────
print('\nMerging u700 and v700 into climatology file...')

u_array = np.stack(u_climo, axis=0)  # (48, lat, lon)
v_array = np.stack(v_climo, axis=0)

ds_existing = xr.open_dataset(CLIMO_FILE)
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

# Write to a temp file then replace, to avoid corrupting the original
# if the write is interrupted.
tmp_out = CLIMO_FILE.with_suffix('.tmp.nc')
ds_merged.to_netcdf(tmp_out)
ds_merged.close()
tmp_out.replace(CLIMO_FILE)

print(f'Done! climatology_wind700_48frame.nc now contains:')
ds_check = xr.open_dataset(CLIMO_FILE)
print('  Variables:', list(ds_check.data_vars))
print('  Dims:', dict(ds_check.sizes))
ds_check.close()