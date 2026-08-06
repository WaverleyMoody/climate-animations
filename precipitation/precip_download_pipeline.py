"""
precip_download_pipeline.py

ERA5 total precipitation (tp) download + climatology pipeline, matching
the structure used for 2m_temp / SLP / 10m wind / precipitable water /
700mb wind:

    1. Download    -- one CDS request per year (22 total, 1979-2000),
                       pulling ERA5's server-side daily SUM of
                       total_precipitation via the
                       derived-era5-single-levels-daily-statistics
                       dataset. This replaces manual hourly download +
                       de-accumulation entirely -- CDS computes the
                       daily accumulated total for us.
    2. Weekly sum  -- for each year, SUM the daily totals into the same
                       four target weeks per month used by the rest of
                       the library (days 1-7, 8-14, 15-21, 22-month_end).
    3. Climatology -- average those weekly TOTALS across all 22 years,
                       producing one 48-frame climatological weekly-total
                       precipitation NetCDF (12 months x 4 weeks),
                       matching the precomputed-climatology format used
                       by animate_2m_temp.py / animate_SLP.py /
                       animate_wind.py / animate_precipitable_water.py /
                       animate_wind700.py.
    4. Units       -- ERA5 tp is archived in metres of water equivalent;
                       converted to cm (x100) here to match the Climate
                       Reanalyzer reference scale (0-50 cm).

CONFIRMED FIELD NAMES (checked against a live 1979 download):
    coords:    number, latitude, longitude, valid_time
    data_vars: tp
    grid:      721 x 1440 (standard ERA5 0.25 deg global)

SECURITY NOTE:
    This script does NOT contain your CDS API key. cdsapi.Client() reads
    credentials from ~/.cdsapirc (a file OUTSIDE this git repo). Create
    that file once with:

        url: https://cds.climate.copernicus.eu/api
        key: <your-key>

    Never commit an API key to climate-animations -- it's a public repo.

USAGE:
    python precip_download_pipeline.py
    (safe to re-run/resume: already-downloaded/processed years are skipped)

    Runs from wherever you place it -- BASE_DIR points at the external
    drive, not anywhere relative to this script's own location.
"""

from __future__ import annotations

import calendar
import logging
import time
from pathlib import Path

import cdsapi
import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR = Path("/Volumes/CLIMATEDATA/precip")
RAW_DIR = BASE_DIR / "raw_daily_sum"      # one file per year: all daily sums
WEEKLY_DIR = BASE_DIR / "weekly_sum"      # one file per year: 48 weekly totals
CLIMO_DIR = BASE_DIR / "climatology"

for d in (RAW_DIR, WEEKLY_DIR, CLIMO_DIR):
    d.mkdir(parents=True, exist_ok=True)

YEARS = list(range(1979, 2001))  # 1979-2000 inclusive
DATASET = "derived-era5-single-levels-daily-statistics"
VARIABLE = "total_precipitation"

MONTHS = [f"{m:02d}" for m in range(1, 13)]
ALL_DAYS = [f"{d:02d}" for d in range(1, 32)]

# Week start days matching the rest of the library (48 frames/year total)
WEEK_STARTS = [1, 8, 15, 22]

MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

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
    """Return (start_day, end_day) for a week block within a month,
    matching the rest of the library: 1-7, 8-14, 15-21, 22-month_end."""
    _, days_in_month = calendar.monthrange(year, month)
    idx = WEEK_STARTS.index(start_day)
    if idx == len(WEEK_STARTS) - 1:
        end_day = days_in_month
    else:
        end_day = WEEK_STARTS[idx + 1] - 1
    return start_day, end_day


def raw_file(year: int) -> Path:
    return RAW_DIR / f"precip_daily_sum_{year}.nc"


def weekly_file(year: int) -> Path:
    return WEEKLY_DIR / f"precip_weekly_sum_{year}.nc"


def weekly_done(year: int) -> bool:
    """Resume check: skip years whose weekly-sum output already exists
    and opens cleanly."""
    f = weekly_file(year)
    if not f.exists():
        return False
    try:
        with xr.open_dataset(f):
            pass
        return True
    except Exception:
        log.warning(f"{f.name} exists but is unreadable -- reprocessing {year}")
        f.unlink(missing_ok=True)
        return False


def _verify_readable(path: Path) -> None:
    """Force an actual data read (not just a metadata open) so corrupted
    HDF5 chunks are caught immediately, rather than ~40 minutes later
    during weekly-sum processing."""
    with xr.open_dataset(path) as ds_check:
        var = "tp" if "tp" in ds_check.data_vars else VARIABLE
        ds_check[var].load()


def download_year(year: int, max_retries: int = 3) -> Path:
    """Download server-side daily SUM of total_precipitation for one full
    year in a single CDS request. Retries on corrupt/truncated downloads."""
    out = raw_file(year)

    if out.exists():
        try:
            _verify_readable(out)
            log.info(f"{year}: raw daily-sum file already present and verified, skipping download")
            return out
        except Exception as e:
            log.warning(f"{year}: existing raw file failed verification ({e}), re-downloading")
            out.unlink(missing_ok=True)

    request = {
        "product_type": "reanalysis",
        "variable": VARIABLE,
        "year": str(year),
        "month": MONTHS,
        "day": ALL_DAYS,
        "daily_statistic": "daily_sum",
        "time_zone": "utc+00:00",
        "frequency": "1_hourly",
        "data_format": "netcdf",
    }

    for attempt in range(1, max_retries + 1):
        try:
            log.info(f"{year}: submitting CDS request (attempt {attempt})")
            client.retrieve(DATASET, request).download(str(out))
            _verify_readable(out)  # forces a real read -- catches corruption now, not later
            log.info(f"{year}: download OK and verified -> {out.name}")
            return out
        except Exception as e:
            log.warning(f"{year}: attempt {attempt} failed ({e})")
            out.unlink(missing_ok=True)
            if attempt < max_retries:
                time.sleep(30)
            else:
                raise


def process_year(year: int) -> None:
    if weekly_done(year):
        log.info(f"{year}: weekly sums already computed, skipping")
        return

    raw_path = download_year(year)
    ds = xr.open_dataset(raw_path, chunks={"time": 50})

    # Confirmed naming from a live download: valid_time / tp. Fallbacks
    # kept in case CDS changes the schema on this newer dataset.
    time_name = "time" if "time" in ds.coords else "valid_time"
    var_name = "tp" if "tp" in ds.data_vars else VARIABLE
    tp = ds[var_name]

    frame_labels = []
    weekly_fields = []

    for month in range(1, 13):
        for start_day in WEEK_STARTS:
            s, e = week_bounds(year, month, start_day)
            mask = (
                (tp[time_name].dt.month == month)
                & (tp[time_name].dt.day >= s)
                & (tp[time_name].dt.day <= e)
            )
            week_slice = tp.sel({time_name: mask})
            # Weekly TOTAL = sum of daily sums across the week's days,
            # converted metres -> cm
            weekly_total_cm = (week_slice.sum(dim=time_name) * 100.0).compute()
            label = f"{MONTH_NAMES[month - 1]} {start_day}"
            frame_labels.append(label)
            weekly_fields.append(weekly_total_cm.values)
            log.info(f"{year} {label}: weekly total computed")

    out_ds = xr.Dataset(
        {"tp_weekly_total_cm": (("frame", "latitude", "longitude"), np.stack(weekly_fields))},
        coords={
            "frame": np.arange(48),
            "frame_label": ("frame", frame_labels),
            "latitude": ds["latitude"].values,
            "longitude": ds["longitude"].values,
        },
        attrs={
            "description": f"Weekly TOTAL precipitation (cm water equiv.), year {year}, "
                            "from ERA5 total_precipitation daily sums",
            "units": "cm",
        },
    )
    out_ds.to_netcdf(weekly_file(year))
    log.info(f"{year}: weekly-sum file written -> {weekly_file(year).name}")
    ds.close()


def build_climatology() -> None:
    climo_out = CLIMO_DIR / "climatology_precip_48frame.nc"
    if climo_out.exists():
        log.info("Climatology file already exists, skipping build step")
        return

    missing = [y for y in YEARS if not weekly_file(y).exists()]
    if missing:
        log.warning(f"Cannot build climatology yet -- missing years: {missing}")
        return

    log.info("Building 48-frame climatology (mean of weekly totals across 1979-2000)...")
    all_years = xr.open_mfdataset(
        [str(weekly_file(y)) for y in YEARS],
        combine="nested",
        concat_dim="year",
    )

    climo = all_years["tp_weekly_total_cm"].mean(dim="year").compute()
    frame_labels = all_years["frame_label"].isel(year=0).values

    climo_ds = xr.Dataset(
        {"tp_weekly_total_cm": (("frame", "latitude", "longitude"), climo.values)},
        coords={
            "frame": np.arange(48),
            "frame_label": ("frame", frame_labels),
            "latitude": all_years["latitude"].values,
            "longitude": all_years["longitude"].values,
        },
        attrs={
            "description": "48-frame climatological weekly TOTAL precipitation "
                            "(cm water equiv.), 1979-2000 mean",
            "units": "cm",
        },
    )
    climo_ds.to_netcdf(climo_out)
    log.info(f"Climatology written -> {climo_out}")


def main() -> None:
    for year in YEARS:
        try:
            process_year(year)
        except Exception as e:
            log.error(f"{year}: failed after retries -- {e}")
            continue
    build_climatology()


if __name__ == "__main__":
    main()