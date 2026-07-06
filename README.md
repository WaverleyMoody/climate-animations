# SDSU Animation Library
#### Version 1.0 (2026)
Our library is a reproduction of the climatology animations from the [University of Washington General Circulation Animations Library](https://animations.atmos.uw.edu/), originally created by Mike Wallace and David Battisti.

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

<table>
  <tr>
    <td width="50%"><video width="100%" src="https://github.com/user-attachments/assets/210ac59f-f6a6-47b1-8614-f767630d222a"></video></td>
  </tr>
</table>
<table width="100%">
  <tr>
    <td><a href="https://youtu.be/NoBMwxeQTJQ">YouTube</a></td>
    <td><a href="2m_temperature/scripts/animate_2m_temp.py">Python Code</a></td>
    <td><a href="2m_temperature/scripts/animate_2m_temp.R">R Code</a></td>
    <td colspan="2"><details><summary><b><i>About:</i></b> ERA5 Weekly Mean 2m Temperature 1979-2000. Shows the climatology of air temperature 2 meters above the surface. Each date</summary> marks the start of the week being averaged.</details></td>
  <tr>
</table>

## 2. Pressure at Mean Sea Level

<table>
  <tr>
    <td width="50%"><video width="100%" src="https://github.com/user-attachments/assets/b84ec4eb-7b53-4f7a-8e69-d69d2b2f5ca7"></video></td>
  </tr>
</table>
<table width="100%">
  <tr>
    <td><a href="https://youtu.be/m6FXkEbiyRc">YouTube</a></td>
    <td><a href="sea_level_pressure/scripts/animate_SLP.py">Python Code</a></td>
    <td><a href="sea_level_pressure/scripts/animate_SLP.R">R Code</a></td>
   <td colspan="2"><details><summary><b><i>About:</i></b> ERA5 Weekly Mean Sea Level Pressure 1979-2000. Shows the climatology of atmospheric pressure adjusted to sea level. Each date </summary>marks the start of the week being averaged.</details></td>
  </tr>
</table>

## 3. Wind Speed at 10m

<table>
  <tr>
    <td width="50%"><video width="100%" src="https://github.com/user-attachments/assets/882e939a-c7c4-4a3d-8cfb-c971e7744a0c"></video></td>
  </tr>
</table>
<table width="100%">
  <tr>
    <td><a href="https://youtu.be/hfYU_kA1568">YouTube</a></td>
    <td><a href="wind_speed/scripts/animate_wind.py">Python Code</a></td>
    <td><a href="wind_speed/scripts/animate_wind.R">R Code</a></td>
    <td colspan="2"><details><summary><b><i>About:</i></b> ERA5 Weekly Mean Surface Wind Speed at 10m 1979-2000. Shows the climatology of wind speed with static arrows showing mean wind</summary>direction. Each date marks the start of the week being averaged.</details></td>
  </tr>
</table>

## 4. Precipitable Water
<table>
  <tr>
    <td width="50%"><video width="100%" src="https://github.com/user-attachments/assets/05ac134b-8d80-464c-8e4e-45f737b57afe"></video></td>
  </tr>
</table>
<table width="100%">
  <tr>
    <td><a href="https://youtu.be/YjMiDUWG2gc">YouTube</a></td>
    <td><a href="precipitable_water/scripts/animate_precipitable_water.py">Python Code</a></td>
    <td><a href="precipitable_water/scripts/animate_precipitable_water.R">R Code</a></td>
    <td colspan="2"><details><summary><b><i>About:</i></b> ERA5 Weekly Mean Precipitable Water 1979-2000. Shows the climatology of total column water vapour. Each date marks the start of</summary>the week being averaged.</details></td>
  </tr>
</table>