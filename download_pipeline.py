import cdsapi
import xarray as xr
import numpy as np
from pathlib import Path
import os

OUTPUT_DIR = Path('/Users/waverleymoody/Downloads/climate_daily_means')
OUTPUT_DIR.mkdir(exist_ok=True)

TEMP_DIR = Path('/Users/waverleymoody/Downloads/climate_temp')
TEMP_DIR.mkdir(exist_ok=True)

c = cdsapi.Client()

YEARS = [str(y) for y in range(1979, 2001)]
MONTHS = [f'{m:02d}' for m in range(1, 13)]
TIMES = ['00:00', '06:00', '12:00', '18:00']

WEEKS = {
    '01': ['01', '02', '03', '04', '05', '06', '07'],
    '08': ['08', '09', '10', '11', '12', '13', '14'],
    '15': ['15', '16', '17', '18', '19', '20', '21'],
    '22': ['22', '23', '24', '25', '26', '27', '28'],
}

VARIABLES = [
    '2m_temperature',
    'mean_sea_level_pressure',
    '10m_u_component_of_wind',
    '10m_v_component_of_wind',
    'total_column_water_vapour',
    'total_precipitation',
]

for year in YEARS:
    for month in MONTHS:
        for week_start, week_days in WEEKS.items():
            out_file = OUTPUT_DIR / f'daily_mean_{year}_{month}_{week_start}.nc'

            if out_file.exists():
                print(f'Already exists, skipping: {out_file.name}')
                continue

            print(f'Downloading: {year}-{month}-{week_start}...')

            temp_file = TEMP_DIR / f'temp_{year}_{month}_{week_start}.grib'

            try:
                c.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'variable': VARIABLES,
                        'year': year,
                        'month': month,
                        'day': week_days,
                        'time': TIMES,
                        'data_format': 'grib',
                        'download_format': 'unarchived',
                    },
                    str(temp_file)
                )

                print(f'Averaging: {year}-{month}-{week_start}...')

                ds_list = []
                filter_keys_list = [
                    {'shortName': '2t'},
                    {'shortName': 'msl'},
                    {'shortName': '10u'},
                    {'shortName': '10v'},
                    {'shortName': 'tcwv'},
                    {'shortName': 'tp'},
                ]

                for fk in filter_keys_list:
                    try:
                        ds = xr.open_dataset(str(temp_file), engine='cfgrib',
                                           backend_kwargs={'filter_by_keys': fk,
                                                         'indexpath': ''})
                        ds_list.append(ds)
                    except Exception:
                        pass

                ds_merged = xr.merge(ds_list, compat='override')
                time_dim = 'time' if 'time' in ds_merged.dims else 'valid_time'
                ds_mean = ds_merged.mean(dim=time_dim)
                ds_mean.to_netcdf(out_file)

                for ds in ds_list:
                    ds.close()

                os.remove(temp_file)
                print(f'Saved weekly mean: {out_file.name}')

            except Exception as e:
                print(f'Error on {year}-{month}-{week_start}: {e}')
                if temp_file.exists():
                    os.remove(temp_file)
                continue

print('All done! Weekly means saved to:', OUTPUT_DIR)