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

**About**
This animation shows the seasonally varying climatology of air temperature 2 metres above the surface, averaged over the period 1979–2000. Each frame represents a weekly mean, cycling through the 1st, 8th, 15th, and 22nd of each month. The animation reveals the dramatic seasonal contrast between the Northern and Southern Hemispheres, the persistent warmth of the tropics, and the extreme cold of the polar regions and high-altitude areas like the Tibetan Plateau and Antarctica.

**Animation**

[![2m Temperature Animation](https://img.youtube.com/vi/NoBMwxeQTJQ/0.jpg)](https://youtu.be/NoBMwxeQTJQ)

**Data**
- Source: ERA5 Reanalysis, European Centre for Medium-Range Weather Forecasts (ECMWF)
- Variable: 2 metre temperature (t2m)
- Period: 1979–2000
- Temporal resolution: Weekly means (days 1, 8, 15, 22 of each month)
- Spatial resolution: 0.25° x 0.25°
- Units: °C (converted from Kelvin)
- Download: [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels)

**Python Script**
Located at: `2m_temperature/scripts/animate_2m_temp.py`

Key libraries: xarray, matplotlib, cartopy, imageio, scipy

**R Script**
Located at: `2m_temperature/scripts/animate_2m_temp.R`

Coming soon.

### 2. Mean Sea Level Pressure

**About**
This animation shows the seasonally varying climatology of atmospheric pressure adjusted to sea level, averaged over the period 1979–2000. It reveals the dominant pressure systems that drive global weather patterns, including the subtropical high pressure belts, the Icelandic and Aleutian lows in winter, and the Asian monsoon low in summer.

**Animation**

[![SLP Animation](https://img.youtube.com/vi/m6FXkEbiyRc/0.jpg)](https://youtu.be/m6FXkEbiyRc)

**Data**
- Source: ERA5 Reanalysis, ECMWF
- Variable: Mean sea level pressure (msl)
- Period: 1979–2000
- Temporal resolution: Weekly means (days 1, 8, 15, 22 of each month)
- Spatial resolution: 0.25° x 0.25°
- Units: mb (converted from Pa)
- Download: [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels)

**Python Script**
Located at: `sea_level_pressure/scripts/animate_SLP.py`

**R Script**
Located at: `sea_level_pressure/scripts/animate_SLP.R`

Coming soon.

### 3. Surface Wind Speed at 10 m

**About**
This animation shows the seasonally varying climatology of wind speed 10 metres above the surface, averaged over the period 1979–2000. Static arrows show the overall mean wind direction across the full period. The animation reveals the trade winds, the westerlies, the jet streams, and the seasonal reversal of winds associated with the monsoon systems.

**Animation**

[![Wind Speed Animation](https://img.youtube.com/vi/hfYU_kA1568/0.jpg)](https://youtu.be/hfYU_kA1568)

**Data**
- Source: ERA5 Reanalysis, ECMWF
- Variables: 10m u-component of wind (u10), 10m v-component of wind (v10)
- Period: 1979–2000
- Temporal resolution: Weekly means (days 1, 8, 15, 22 of each month)
- Spatial resolution: 0.25° x 0.25°
- Units: m/s
- Download: [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels)

**Python Script**
Located at: `wind_speed/scripts/animate_wind.py`

**R Script**
Located at: `wind_speed/scripts/animate_wind.R`

Coming soon.

### 4. Precipitable Water

**About**
This animation shows the seasonally varying climatology of total column water vapour, averaged over the period 1979–2000. It reveals the concentration of atmospheric moisture in the tropics, the dry subtropical deserts, and the seasonal migration of the moist tropical belt associated with the Intertropical Convergence Zone (ITCZ).

**Animation**

[![Precipitable Water Animation](https://img.youtube.com/vi/YjMiDUWG2gc/0.jpg)](https://youtu.be/YjMiDUWG2gc)

**Data**
- Source: ERA5 Reanalysis, ECMWF
- Variable: Total column water vapour (tcwv)
- Period: 1979–2000
- Temporal resolution: Weekly means (days 1, 8, 15, 22 of each month)
- Spatial resolution: 0.25° x 0.25°
- Units: kg/m²
- Download: [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels)

**Python Script**
Located at: `precipitable_water/scripts/animate_precipitable_water.py`

**R Script**
Located at: `precipitable_water/scripts/animate_precipitable_water.R`

Coming soon.

---

## References
Wallace, J.M. et al. (2023). The Atmospheric General Circulation. Cambridge University Press.
