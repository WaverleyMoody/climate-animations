import cdsapi
import xarray as xr
from pathlib import Path
import calendar

# ── Config ───────────────────────────────────────────────────────────────────
OUTPUT_DIR = Path.home() / 'Downloads' / 'climate_daily_means_700mb_wind'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TEMP_GRIB = OUTPUT_DIR / '_temp_download.grib'

YEARS = range(1979, 2001)  # 1979-2000 inclusive
MONTHS = range(1, 13)

# Same 4 target weeks per month as the existing pipeline: week_start is the
# label used in the output filename; day_range is the actual days pulled
# and averaged for that week.
WEEKS = [
    (1, range(1, 8)),    # days 1-7
    (8, range(8, 15)),   # days 8-14
    (15, range(15, 22)), # days 15-21
    (22, range(22, 29)), # days 22-28
]

TIMES = ['00:00', '06:00', '12:00', '18:00']  # 6-hourly, matching the existing pipeline
PRESSURE_LEVEL = '700'

client = cdsapi.Client()


def download_week(year, month, week_start, days):
    """Download one week of 700mb u/v wind, average to a weekly mean,
    save as NetCDF, clean up the raw GRIB."""
    out_path = OUTPUT_DIR / f'daily_mean_700mb_{year}_{month:02d}_{week_start:02d}.nc'
    if out_path.exists():
        print(f'  Skipping {out_path.name} (already exists)')
        return

    # clamp day range to the actual number of days in this month (handles
    # Feb and the 22-28 week not overrunning into a nonexistent day 29-31
    # for months with fewer days)
    days_in_month = calendar.monthrange(year, month)[1]
    day_list = [f'{d:02d}' for d in days if d <= days_in_month]
    if not day_list:
        print(f'  Skipping {year}-{month:02d} week {week_start} (no valid days)')
        return

    request = {
        'product_type': 'reanalysis',
        'variable': ['u_component_of_wind', 'v_component_of_wind'],
        'pressure_level': PRESSURE_LEVEL,
        'year': str(year),
        'month': f'{month:02d}',
        'day': day_list,
        'time': TIMES,
        'data_format': 'grib',
        'download_format': 'unarchived',
    }

    print(f'  Downloading {year}-{month:02d} week {week_start}...')
    client.retrieve('reanalysis-era5-pressure-levels', request, str(TEMP_GRIB))

    # u and v both live on the same typeOfLevel (isobaricInhPa) at the same
    # single pressure level, so unlike the single-levels pipeline (which had
    # to filter_by_keys per variable to separate surface/meanSea/etc.), this
    # opens cleanly as one dataset with both variables already present.
    ds = xr.open_dataset(
        TEMP_GRIB, engine='cfgrib',
        backend_kwargs={'filter_by_keys': {'typeOfLevel': 'isobaricInhPa'}}
    )

    weekly_mean = ds.mean(dim='time')
    weekly_mean = weekly_mean.expand_dims(
        time=[xr.cftime_range(start=f'{year}-{month:02d}-{week_start:02d}', periods=1)[0]]
    )
    weekly_mean.to_netcdf(out_path)
    ds.close()

    TEMP_GRIB.unlink(missing_ok=True)
    # cfgrib also leaves a .idx sidecar file next to the source -- clean
    # that up too, same as the original pipeline should already be doing
    idx_file = Path(str(TEMP_GRIB) + '.923a8.idx')
    if idx_file.exists():
        idx_file.unlink()

    print(f'  Saved {out_path.name}')


def main():
    total = len(list(YEARS)) * len(list(MONTHS)) * len(WEEKS)
    done = 0
    for year in YEARS:
        for month in MONTHS:
            for week_start, days in WEEKS:
                download_week(year, month, week_start, days)
                done += 1
                print(f'Progress: {done}/{total}')

    print('All downloads complete.')


if __name__ == '__main__':
    main()