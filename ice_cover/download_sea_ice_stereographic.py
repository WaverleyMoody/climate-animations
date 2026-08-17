# SDSU Climate Informatics Lab / San Diego State University / by Waverley Moody / Supervised by Distinguished Professor Samuel Shen / Python Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation Animations Library by Professor John Michael Wallace.
#
# Script: download_sea_ice_stereographic.py
# Description: Downloads ERA5 sea_ice_cover for the 24 biweekly target dates of 2004
#              (1st and 15th of each month), for use as the sea ice overlay in the
#              northern/southern hemisphere stereographic seasonal-cycle animation.
# Note: This script is projection-agnostic. Reprojection to the polar stereographic
#       grids (north and south) happens downstream in the per-hemisphere render scripts.

import cdsapi
import xarray as xr
from pathlib import Path

OUTPUT_DIR = Path("/Volumes/CLIMATEDATA/era5_sea_ice/2004")
OUTPUT_FILE = OUTPUT_DIR / "era5_sea_ice_cover_2004_biweekly.nc"

MONTHS = [f"{m:02d}" for m in range(1, 13)]
DAYS = ["01", "15"]
TIME = ["00:00"]  # sea ice cover has negligible diurnal cycle; one timestep/day is sufficient


def download_sea_ice():
    """Single CDS batch request for all 24 target dates in 2004."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_FILE.exists():
        print(f"File already exists at {OUTPUT_FILE}, skipping download.")
        return

    c = cdsapi.Client()

    print("Requesting ERA5 sea_ice_cover for 2004 (1st and 15th of each month)...")
    c.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": "reanalysis",
            "variable": "sea_ice_cover",
            "year": "2004",
            "month": MONTHS,
            "day": DAYS,
            "time": TIME,
            "format": "netcdf",
        },
        str(OUTPUT_FILE),
    )
    print(f"Download complete: {OUTPUT_FILE}")


def verify_file():
    """Deep-load check to catch corrupted HDF5 files that pass a shallow open."""
    print("Verifying downloaded file integrity...")
    ds = xr.open_dataset(OUTPUT_FILE)
    ds.load()  # forces a full read; open_dataset() alone can succeed on a corrupted file
    time_dim = "valid_time" if "valid_time" in ds.dims else "time"
    n_times = ds.dims[time_dim]
    expected = len(MONTHS) * len(DAYS)
    status = "OK" if n_times == expected else "MISMATCH — check request/response"
    print(f"Verified {status} — {n_times} time steps present (expected {expected}).")
    ds.close()


if __name__ == "__main__":
    download_sea_ice()
    verify_file()