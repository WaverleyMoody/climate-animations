#!/usr/bin/env python3
"""
split_wind_variable.py

Wind-only variant of split_variables_by_type.py.

You've already produced 2m_temp_data.nc, slp_data.nc, and precip_water_data.nc
from earlier runs — this script skips re-deriving those and only builds
wind_data.nc (u10 + v10 + derived wind_speed) from the raw daily_mean_*.nc files.

It also drops all non-wind variables at file-open time (via preprocess),
so dask never has to carry t2m/msl/tcwv/tp through memory at all.

Input:  <UPDATE_ME>/daily_mean_YYYY_MM_DD.nc   (1,056 files, now on external disk)
Output: ~/Downloads/climate_data_by_variable/wind_data.nc

Run:
    python split_wind_variable.py
"""

import os
import re
import glob
import xarray as xr
import pandas as pd
import numpy as np

# UPDATE THIS to wherever the disk mounted the folder, e.g.:
# "/Volumes/MyDisk/climate_daily_means/"
INPUT_DIR = "/Volumes/CLIMATEDATA/climate_daily_means/"
OUTPUT_DIR = os.path.expanduser("~/Downloads/climate_data_by_variable/")

WIND_COMPONENTS = ("u10", "v10")
WIND_OUTPUT_NAME = "wind_speed"
OUTPUT_FILENAME = "wind_data.nc"

FNAME_RE = re.compile(r"daily_mean_(\d{4})_(\d{2})_(\d{2})\.nc$")


def parse_date_from_filename(path):
    m = FNAME_RE.search(os.path.basename(path))
    if not m:
        return None
    y, mo, d = m.groups()
    return pd.Timestamp(int(y), int(mo), int(d))


def keep_only_wind(ds):
    """Preprocess hook for open_mfdataset: drop everything except u10/v10
    as soon as each file is opened, so other variables never get loaded."""
    drop_vars = [v for v in ds.data_vars if v not in WIND_COMPONENTS]
    return ds.drop_vars(drop_vars, errors="ignore")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    files = sorted(glob.glob(os.path.join(INPUT_DIR, "daily_mean_*.nc")))
    if not files:
        raise SystemExit(f"No files found in {INPUT_DIR} — check the path.")
    print(f"Found {len(files)} input files.")

    with xr.open_dataset(files[0]) as sample:
        missing = [v for v in WIND_COMPONENTS if v not in sample.data_vars]
        if missing:
            raise SystemExit(
                f"Expected wind components {WIND_COMPONENTS} but missing {missing}. "
                f"Variables in sample file: {list(sample.data_vars)}"
            )

    dates = [parse_date_from_filename(f) for f in files]
    if any(d is None for d in dates):
        bad = [f for f, d in zip(files, dates) if d is None][0]
        raise SystemExit(f"Could not parse a date from filename: {bad} — check FNAME_RE.")

    print("Opening all files (lazy, dask-backed, wind-only)...")
    ds_all = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="time",
        coords="minimal",
        compat="override",
        chunks={"time": 20},
        preprocess=keep_only_wind,
    )
    ds_all = ds_all.assign_coords(time=("time", dates))

    out_path = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME)
    print(f"Writing {out_path} ...")

    u, v = ds_all[WIND_COMPONENTS[0]], ds_all[WIND_COMPONENTS[1]]
    speed = np.hypot(u, v).rename(WIND_OUTPUT_NAME)
    data = xr.merge([u, v, speed])

    data.to_netcdf(out_path, engine="netcdf4")
    print(f"  -> done ({os.path.getsize(out_path) / 1e6:.1f} MB)")
    print("\nWind variable file written to:", out_path)


if __name__ == "__main__":
    main()