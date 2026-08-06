"""
download_pipeline_wind700.py

ERA5 700mb wind speed download + climatology pipeline for the
SDSU Animation Library.

Pipeline:
  1. Download hourly u/v wind components at 700hPa, one CDS request
     per year-month, for 1979-2000.
  2. Compute daily mean wind speed (mean of sqrt(u^2+v^2) at each
     hourly timestep -- NOT the magnitude of the mean vector).
  3. Compute weekly means per year, using the same week definition as
     the rest of the library: weeks start on day 1, 8, 15, 22 of each
     month (the last "week" of a month runs to month-end).
  4. Compute the climatological weekly mean across all years, producing
     a single 48-frame NetCDF file (frame 0 = Jan week 1 ... frame 47
     = Dec week 4), matching the precomputed-climatology format used
     by animate_2m_temp.py / animate_SLP.py / animate_wind.py /
     animate_precipitable_water.py.

SECURITY NOTE:
  This script does NOT contain your CDS API key. cdsapi.Client() reads
  credentials from ~/.cdsapirc (a file OUTSIDE this git repo). Create
  that file once with:

      url: https://cds.climate.copernicus.eu/api
      key: <your-key>

  Never commit an API key to climate-animations -- it's a public repo.

USAGE:
  python download_pipeline_wind700.py
  (safe to re-run/resume: already-downloaded/processed files are skipped)
"""

from __future__ import annotations

import calendar
import logging
from pathlib import Path

import cdsapi
import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR = Path("/Volumes/CLIMATEDATA/wind700")
RAW_DIR = BASE_DIR / "raw_hourly"
DAILY_DIR = BASE_DIR / "daily_mean"
WEEKLY_DIR = BASE_DIR / "weekly_mean"
CLIMO_DIR = BASE_DIR / "climatology"

for d in (RAW_DIR, DAILY_DIR, WEEKLY_DIR, CLIMO_DIR):
    d.mkdir(parents=True, exist_ok=True)

YEARS = list(range(1979, 2001))  # 1979-2000 inclusive
PRESSURE_LEVEL = "700"

# Free disk space as you go. Set False if you want to keep raw hourly
# files around for debugging/reuse.
DELETE_RAW_AFTER_DAILY = True

# Week start days matching the rest of the library (48 frames/year total)
WEEK_STARTS = [1, 8, 15, 22]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(BASE_DIR / "pipeline.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# Reads url/key from ~/.cdsapirc -- do not pass credentials here
client = cdsapi.Client()


def week_bounds(year: int, month: int, start_day: int) -> tuple[int, int]:
    """Return (start_day, end_day) for a week block within a month."""
    last_day = calendar.monthrange(year, month)[1]
    end_day = last_day if start_day == 22 else start_day + 6
    return start_day, end_day


# ---------------------------------------------------------------------------
# Step 1: Download raw hourly data (one request per year-month)
# ---------------------------------------------------------------------------
def download_month(year: int, month: int) -> Path:
    out_path = RAW_DIR / f"wind700_{year}_{month:02d}.nc"
    if out_path.exists():
        log.info(f"Raw file exists, skipping download: {out_path.name}")
        return out_path

    n_days = calendar.monthrange(year, month)[1]
    log.info(f"Requesting {year}-{month:02d} from CDS...")
    client.retrieve(
        "reanalysis-era5-pressure-levels",
        {
            "product_type": "reanalysis",
            "variable": ["u_component_of_wind", "v_component_of_wind"],
            "pressure_level": PRESSURE_LEVEL,
            "year": str(year),
            "month": f"{month:02d}",
            "day": [f"{d:02d}" for d in range(1, n_days + 1)],
            "time": [f"{h:02d}:00" for h in range(24)],
            "data_format": "netcdf",
            "download_format": "unarchived",
        },
        str(out_path),
    )
    log.info(f"Downloaded {out_path.name}")
    return out_path


# ---------------------------------------------------------------------------
# Step 2: Daily means
# ---------------------------------------------------------------------------
def daily_means_done(year: int, month: int) -> bool:
    """True if every daily-mean file for this month already exists."""
    return all(
        (DAILY_DIR / f"daily_mean_wind700_{year}_{month:02d}_{d:02d}.nc").exists()
        for d in range(1, calendar.monthrange(year, month)[1] + 1)
    )


def compute_daily_means(year: int, month: int, raw_path: Path) -> None:
    """Per-day mean wind speed from hourly u/v (mean of speed, not
    speed of the mean vector)."""
    if daily_means_done(year, month):
        log.info(f"Daily means already computed for {year}-{month:02d}, skipping")
        if DELETE_RAW_AFTER_DAILY and raw_path.exists():
            raw_path.unlink()
        return

    ds = xr.open_dataset(raw_path, chunks={"valid_time": 24})

    speed = np.sqrt(ds["u"] ** 2 + ds["v"] ** 2)
    speed.name = "wind_speed_700mb"

    daily = speed.resample(valid_time="1D").mean(skipna=True).compute()

    for day_val in daily.valid_time.values:
        day_da = daily.sel(valid_time=day_val)
        day_str = np.datetime_as_string(day_val, unit="D").replace("-", "_")
        out_path = DAILY_DIR / f"daily_mean_wind700_{day_str}.nc"
        day_da.to_netcdf(out_path)

    ds.close()
    log.info(f"Computed daily means for {year}-{month:02d}")

    if DELETE_RAW_AFTER_DAILY:
        raw_path.unlink()
        log.info(f"Deleted raw file {raw_path.name} to free disk space")


# ---------------------------------------------------------------------------
# Step 3: Weekly means (per year)
# ---------------------------------------------------------------------------
def compute_weekly_means(year: int) -> None:
    for month in range(1, 13):
        for start_day in WEEK_STARTS:
            start_day, end_day = week_bounds(year, month, start_day)
            out_path = WEEKLY_DIR / f"weekly_mean_wind700_{year}_{month:02d}_{start_day:02d}.nc"
            if out_path.exists():
                continue

            day_files = [
                DAILY_DIR / f"daily_mean_wind700_{year}_{month:02d}_{d:02d}.nc"
                for d in range(start_day, end_day + 1)
            ]
            day_files = [f for f in day_files if f.exists()]
            if not day_files:
                log.warning(f"No daily files found for week {year}-{month:02d}-{start_day:02d}")
                continue

            ds = xr.open_mfdataset(day_files, combine="nested", concat_dim="valid_time")
            week_mean = ds["wind_speed_700mb"].mean(dim="valid_time").compute()
            week_mean.to_netcdf(out_path)
            ds.close()

    log.info(f"Computed weekly means for {year}")


# ---------------------------------------------------------------------------
# Step 4: Climatological weekly mean across all years -> 48-frame file
# ---------------------------------------------------------------------------
def compute_climatology() -> Path:
    frames = []
    frame_idx = 0
    for month in range(1, 13):
        for start_day in WEEK_STARTS:
            files = [
                WEEKLY_DIR / f"weekly_mean_wind700_{y}_{month:02d}_{start_day:02d}.nc"
                for y in YEARS
            ]
            files = [f for f in files if f.exists()]
            if len(files) < len(YEARS) * 0.9:
                log.warning(
                    f"Only {len(files)}/{len(YEARS)} years available for "
                    f"{month:02d}-{start_day:02d}; climatology may be biased"
                )

            ds = xr.open_mfdataset(files, combine="nested", concat_dim="year")
            climo = ds["wind_speed_700mb"].mean(dim="year").compute()
            climo = climo.expand_dims(frame=[frame_idx])
            frames.append(climo)
            ds.close()
            frame_idx += 1

    ds_clim = xr.concat(frames, dim="frame")
    ds_clim.name = "wind_speed_700mb"
    out_path = CLIMO_DIR / "climatology_wind700_48frame.nc"
    ds_clim.to_netcdf(out_path)
    log.info(f"Saved climatology file: {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
FAILED_MONTHS_LOG = BASE_DIR / "failed_months.txt"


def process_month(year: int, month: int, max_retries: int = 1) -> bool:
    """Download + compute daily means for one month. If the raw file is
    corrupt (truncated download, HDF5 read errors, etc.), delete it and
    retry the download up to max_retries times. Returns True on success,
    False if it still fails after retries (logged for manual follow-up)."""
    if daily_means_done(year, month):
        log.info(f"{year}-{month:02d} already fully processed, skipping download")
        return True

    attempt = 0
    while True:
        raw_path = download_month(year, month)
        try:
            compute_daily_means(year, month, raw_path)
            return True
        except (OSError, RuntimeError) as e:
            log.error(f"Failed processing {year}-{month:02d} (attempt {attempt + 1}): {e}")
            if raw_path.exists():
                raw_path.unlink()  # remove corrupted file so next attempt re-downloads
            attempt += 1
            if attempt > max_retries:
                log.error(f"Giving up on {year}-{month:02d} after {attempt} attempts")
                with open(FAILED_MONTHS_LOG, "a") as f:
                    f.write(f"{year}-{month:02d}\n")
                return False


def main() -> None:
    for year in YEARS:
        for month in range(1, 13):
            process_month(year, month)
        compute_weekly_means(year)

    compute_climatology()
    log.info("Pipeline complete.")
    if FAILED_MONTHS_LOG.exists():
        log.warning(
            f"Some months failed after retries — see {FAILED_MONTHS_LOG} "
            "and rerun the script to pick them up (weekly means/climatology "
            "will be incomplete for those months until then)."
        )


if __name__ == "__main__":
    main()