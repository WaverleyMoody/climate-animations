# SDSU Climate Dynamics Animation Library

**Version 1.0.0 (2026)**
SDSU Climate Informatics Lab · San Diego State University
By Waverley Moody · Supervised by Distinguished Professor Samuel Shen

An open-source reproduction of the [University of Washington General Circulation Animations Library](https://animations.atmos.uw.edu/) by Professor John Michael Wallace. Each animation is rebuilt from publicly available reanalysis and satellite data, with both **Python** and **R** code so anyone can regenerate it.

### ▶ [View the library](https://waverleymoody.github.io/climate-animations/)

The website is the main way to explore the library: every animation can be played in the browser (with speed and loop controls), switched between map projections, and comes with links to its data, Python code, R code, YouTube video, and a downloadable MP4.

---

## Animations

### Seasonally Varying Climatologies

| Animation | Data | Period | Projections | Code |
|---|---|---|---|---|
| Temperature at 2 m | ERA5 | 1979–2000 weekly climatology | Plate Carrée, Robinson, Foucaut, Nicolosi | [`2m_temperature/`](2m_temperature/scripts) |
| Sea Level Pressure | ERA5 | 1979–2000 weekly climatology | Plate Carrée, Robinson, Foucaut, Nicolosi | [`sea_level_pressure/`](sea_level_pressure/scripts) |
| Surface Wind at 10 m | ERA5 | 1979–2000 weekly climatology | Plate Carrée, Robinson, Foucaut, Nicolosi | [`wind_speed/`](wind_speed/scripts) |
| Wind Speed at 700 mb | ERA5 | 1979–2000 weekly climatology | Plate Carrée, Robinson, Foucaut, Nicolosi | [`wind_700mb/`](wind_700mb) |
| Precipitable Water | ERA5 | 1979–2000 weekly climatology | Plate Carrée, Robinson, Foucaut, Nicolosi | [`precipitable_water/`](precipitable_water/scripts) |
| Precipitation | ERA5 | 1979–2000 weekly climatology | Plate Carrée, Robinson, Foucaut, Nicolosi | [`precipitation/`](precipitation) |
| Lightning | WWLLN Global Lightning Climatology (WGLC) | 2010–2025 daily climatology | Plate Carrée, Robinson, Foucaut, Nicolosi | [`lightning/`](lightning) |
| Sea Ice — Polar Stereographic Seasonal Cycle | ERA5 sea ice cover over NASA Blue Marble (2004) | Seasonal cycle | North & South Polar Stereographic | [`ice_cover/`](ice_cover) |

### Evolving Fields of Surface Temperature

| Animation | Data | Period | Projection | Code |
|---|---|---|---|---|
| Land Surface Temperature (Day) | NASA MODIS MOD11C3 v6.1 | Feb 2000 – Dec 2013, monthly | Mollweide | [`land_surface_temperature/`](land_surface_temperature/scripts) |
| Sea Surface Temperature | NASA JPL MUR25 | 2016–2020, daily | Orthographic (spinning globe) | [`sea_surface_temp/`](sea_surface_temp/scripts) |

### The Biosphere

| Animation | Data | Period | Views | Code |
|---|---|---|---|---|
| Global Biosphere | Copernicus GlobColour ocean chlorophyll + MODIS MOD13C2 land NDVI | 2000–2010, monthly | Australia, North Pacific | [`biosphere/`](biosphere) |
| Gross Primary Productivity | NASA MODIS MOD17A2HGF v061 | 2000–2009, 8-day | Plate Carrée, Robinson, Foucaut, Nicolosi | [`global_primary_productivity/`](global_primary_productivity) |

### The Cryosphere

| Animation | Data | Period | View | Code |
|---|---|---|---|---|
| Sea Ice Motion | IABP drifting buoys over NOAA/NSIDC G02202 V6 sea ice concentration | 1979–2022, daily | Arctic | [`arctic_buoy/`](arctic_buoy) |
| Snow Cover & Sea Ice | NOAA/NSIDC IMS 4 km snow & ice + G02202 V6 Antarctic sea ice over NASA Blue Marble | 2019–2021, 5-day steps | North & South polar globes | [`snow_and_sea_ice/`](snow_and_sea_ice) |

### Fires

| Animation | Data | Period | Projections | Code |
|---|---|---|---|---|
| Fires Measured by Carbon Emission | GFED5.1 daily emissions | Mar 2003 – Jan 2022, weekly | Plate Carrée, Robinson, Foucaut, Nicolosi | [`fires/`](fires) |

### Aerosol Transport

| Animation | Data | Period | Projection | Code |
|---|---|---|---|---|
| Saharan Dust, June 2020 | NASA MERRA-2 dust optical depth (DUEXTTAU) | June 2020, hourly | Plate Carrée | [`saharan_dust/`](saharan_dust) |
| Smoke from Arctic Boreal Forest Wildfires | NASA MERRA-2 total aerosol optical depth (TOTEXTTAU) | Jun–Aug 2019, hourly | Orthographic | [`arctic_smoke/`](arctic_smoke) |
| Global Aerosols | NASA MERRA-2 aerosol optical thickness by species | Aug 2018 – Jan 2020, 3-hourly | Plate Carrée | [`global_aerosols/`](global_aerosols) |

### Transport of Trace Gases: Carbon Species

| Animation | Data | Period | Projection | Code |
|---|---|---|---|---|
| CO from Arctic Wildfires | NASA MERRA-2 CO column (COCL) | Jun–Aug 2019, hourly | Orthographic | [`co_arctic/`](co_arctic) |
| CO, December 2006 | NASA MERRA-2 surface CO (COSC) | Dec 2006, hourly | Plate Carrée (Pacific-centered) | [`co_surface/`](co_surface) |
| CO and CO₂, 2006 | CAMS EGG4 column-mean CO₂ + MERRA-2 CO column | 2006, 3-hourly | Plate Carrée | [`co_co2_2006/`](co_co2_2006) |

Full data-source descriptions, units, and processing notes for each animation are on the [website](https://waverleymoody.github.io/climate-animations/).

---

## Repository Layout

```
climate-animations/
├── docs/                    # GitHub Pages site (index.html, style.css, script.js, data pages)
│   └── videos/              # Rendered MP4s shown on the site
├── <animation folder>/      # One folder per animation: Python (.py) and R (.R) scripts,
│                            #   one script per projection where applicable
└── download_pipeline.py     # ERA5 download + weekly averaging for the climatology animations
```

Every script opens with a header naming the lab, author, supervisor, version, the animation it produces, its data source and date range, and its projection.

## Data

Processed input data is published as [GitHub Releases](https://github.com/WaverleyMoody/climate-animations/releases) rather than stored in the repository:

| Release | Contents |
|---|---|
| `climatology_v1.0_data`, `climatology_v1.1_data` | ERA5 climatologies (2m temperature, SLP, winds, precipitable water, precipitation, sea ice cover), lightning, Blue Marble imagery, GPP, land surface temperature, and sea surface temperature |
| `biosphere_v1.0_data` | Merged ocean chlorophyll + land NDVI, one NetCDF per year (2000–2010) |
| `cryosphere_v1.0_data` | Sea ice concentration, filled gaps, and projected buoy tracks; snow & sea ice inputs |
| `fire_climatology_v1.0` | GFED5.1 weekly carbon emissions, 2003–2022 |
| `aerosol_transport_v1.0` | MERRA-2 dust, smoke, and global aerosol files |
| `carbon_species_v1.0` | MERRA-2 CO files and the CO/CO₂ 2006 bundle |

Each animation's **Data Download** button on the website links directly to the files it uses. Original sources are the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/) (ERA5, CAMS), [NASA Earthdata](https://www.earthdata.nasa.gov/) (MODIS, MERRA-2, MUR), [Copernicus Marine Service](https://marine.copernicus.eu/) (GlobColour), [NSIDC](https://nsidc.org/) (sea ice, IMS), [IABP](https://iabp.apl.uw.edu/), [GFED](https://www.globalfiredata.org/), and [WGLC/WWLLN](https://wwlln.net/).

## How to Reproduce

### 1. Set up the environment

**Python** (3.11+, a conda environment is recommended):

```bash
conda create -n climate -c conda-forge python=3.11 \
    xarray dask netcdf4 h5netcdf numpy pandas scipy matplotlib cartopy \
    pyproj shapely pyhdf imageio imageio-ffmpeg tqdm
conda activate climate
```

`pyhdf` is required for the MODIS HDF4 files. For downloading raw data yourself, also install `cdsapi`, `earthaccess`, and `copernicusmarine`.

**R:**

```r
install.packages(c("terra", "ncdf4", "sf", "rnaturalearth", "av"))
```

### 2. Get the data

The fastest route is to download the processed files from the [Releases](https://github.com/WaverleyMoody/climate-animations/releases) page (or the Data Download links on the website). To rebuild from raw sources instead:

- **ERA5 climatologies** — add your CDS key to `~/.cdsapirc`, then run `python download_pipeline.py`.
- **NASA products** (MODIS, MERRA-2, MUR) — log in with a free [Earthdata](https://urs.earthdata.nasa.gov/) account via `earthaccess`.
- **Ocean chlorophyll** — log in to Copernicus Marine via `copernicusmarine`.

### 3. Render an animation

Open the script for the animation and projection you want, point its input path at the downloaded data, and run it:

```bash
python fires/animate_fires_robinson.py
```

or, in R:

```bash
Rscript fires/animate_fires_robinson.R
```

Some scripts include a test mode that renders a few frames for review before a full render; see the notes in each script's header. Full renders can take a long time — on macOS, wrap them in `caffeinate -i` to keep the machine awake.

---

## Credits

- **Original library:** Professor John Michael Wallace, University of Washington — [General Circulation Animations Library](https://animations.atmos.uw.edu/)
- **Reproduction:** Waverley Moody, SDSU Climate Informatics Lab
- **Supervision:** Distinguished Professor Samuel Shen, San Diego State University

Several animations follow NASA Scientific Visualization Studio originals referenced by the UW library. Data credits for each animation are listed on the website.
