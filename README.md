# Climate Animations

Reproductions of the seasonally varying climatology animations from the [University of Washington General Circulation Animations Library](https://animations.atmos.uw.edu/), originally created by Mike Wallace and David Battisti.

## Data Source
All animations are based on **ERA5 reanalysis data** from the [Copernicus Climate Data Store (CDS)](https://cds.climate.copernicus.eu/), covering the period **1979–2000**. Each frame represents a **weekly mean** (weeks starting on the 1st, 8th, 15th, and 22nd of each month).

## How to Reproduce

### 1. Set up your environment
Install Python 3.13+ and run:

    python -m venv venv
    source venv/bin/activate
    pip install xarray matplotlib cartopy numpy imageio scipy dask cfgrib cdsapi netCDF4 pandas

### 2. Set up your CDS API key
Create a file at ~/.cdsapirc with your Copernicus credentials:

    url: https://cds.climate.copernicus.eu/api
    key: YOUR-API-KEY-HERE

### 3. Download the data
Run the pipeline script to download and average the ERA5 data:

    python downlaod_pipeline.py

### 4. Run the animation scripts
Each animation has its own script in the corresponding folder.

---

## Animations

### 1. Temperature at 2 m
Shows the seasonally varying climatology of air temperature 2 metres above the surface. Range: -60C to 50C.

### 2. Mean Sea Level Pressure
Shows the seasonally varying climatology of atmospheric pressure adjusted to sea level. Range: 970-1040 mb.

### 3. Surface Wind Speed at 10 m
Shows the seasonally varying climatology of wind speed 10 metres above the surface, with static wind direction arrows showing the overall mean flow. Range: 0-16 m/s.

### 4. Precipitable Water
Shows the seasonally varying climatology of total column water vapour. Range: 0-80 kg/m2.

---

## References
Wallace, J.M. et al. (2023). The Atmospheric General Circulation. Cambridge University Press.
