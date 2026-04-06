# Spring temperature remains the dominant driver of leaf-out in a warming world

**Zohner, C.M., Wu, Z., Mo, L., Crowther, T.W., Fu, Y.H., Renner, S.S., Vitasse, Y. & Rebindaine, D.**

*Science* (submitted)

---

## Overview

This repository contains all code used to analyse spring leaf-out phenology across Northern Hemisphere deciduous broadleaf forests. The study integrates ground-based observations from the PEP725 network with satellite-derived phenology from MODIS, and uses Bayesian hierarchical models and growing-degree day simulations to quantify the drivers of spring leaf-out and the mechanisms behind apparent declines in temperature sensitivity.

---

## Repository structure

Scripts are numbered to reflect the order of execution within each analytical branch (`1_` data preparation → `2_` driver assembly → `3_` statistical models → `4_` figures). HPC scripts (`.r`) are designed to run on a computing cluster; all other scripts (`.Rmd`, `.qmd`) can be run locally.

### PEP725 analysis (`PEP725/`)

The PEP725 branch processes ground-based spring leaf-out observations for 10 European woody species from 1951–2023. Data preparation scripts clean the raw PEP725 records, extract daily climate variables from E-OBS (temperature, radiation, humidity), and compute daylength for each site. Driver assembly scripts then extract species- and site-specific preseason climate windows, compute seasonal temperature summaries (TQ2–TQ4) and previous-year covariates, and run four variants of growing-degree day (GDD) models — forcing-only, low-chilling, high-chilling, and chilling + photoperiod — for each site × species time series on the HPC. A leave-one-out cross-validation framework evaluates temporal trends in GDD model performance. Analysis scripts fit Bayesian hierarchical models (brms/CmdStan) with site random intercepts and year random slopes, including main models, robustness checks (Student-t, unscaled response), and sensitivity analyses under alternative model specifications. Figure scripts generate all moving-window trend plots, preseason length diagnostics, simulation assumption diagrams, and model performance figures.

### MODIS analysis (`MODIS/`)

The MODIS branch processes satellite-derived phenology for deciduous broadleaf forests across North America, Europe, and Asia (2001–2023). Data preparation scripts extract land cover classifications (MCD12Q1), quality-filter phenological metrics (MCD12Q2; SOS₁₅, SOS₅₀, EOS₅₀), compute daylength per pixel, and process daily GLDAS-2.1 climate fields (temperature, radiation, specific humidity). A Python notebook handles bulk GLDAS NetCDF extraction. Driver assembly converts specific to relative humidity, computes preseason climate windows and seasonal temperature summaries, and assembles the complete site × year driver table. Analysis and figure scripts mirror the PEP725 branch: main Bayesian models, robustness checks, and sensitivity analyses under alternative specifications.

---

## Data sources

All input data are publicly available and are **not** included in this repository due to file size.

| Dataset | Description | Access |
|---------|-------------|--------|
| PEP725 | Pan-European Phenology Database; spring leaf-out and autumn senescence for 10 woody species, 1951–2023 | [www.pep725.eu](http://www.pep725.eu) |
| MCD12Q2 v006 | MODIS Land Cover Dynamics (phenology), 500 m, 2001–2023 | [NASA EOSDIS LP DAAC](https://doi.org/10.5067/MODIS/MCD12Q2.006) |
| MCD12Q1 v006 | MODIS Land Cover Type, 500 m, 2001–2023 | [NASA EOSDIS LP DAAC](https://lpdaac.usgs.gov) |
| GLDAS-2.1 | Global Land Data Assimilation System daily climate, 0.25°, 2001–2023 | [NASA GES DISC](https://disc.gsfc.nasa.gov) |
| E-OBS v30.0e | Gridded European daily climate, 0.1°, 1950–2023 | [Copernicus/ECA&D](https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php) |
| WWF Biomes | Terrestrial Ecoregions of the World | [WWF](https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world) |

---

## Software requirements

All R code was developed and tested with **R 4.4**. Key packages:

| Package | Purpose |
|---------|---------|
| `brms`, `cmdstanr` | Bayesian hierarchical models via Stan |
| `lme4` | Mixed-effects models for moving-window analyses |
| `chillR` | Hourly temperature reconstruction and chilling calculations |
| `data.table`, `tidyverse` | Data wrangling |
| `terra`, `raster`, `sf` | Spatial data handling |
| `geosphere` | Photoperiod calculation |
| `ggplot2`, `patchwork` | Figures |

---

## Running the code

Scripts should be run in numerical order within each branch (1 → 2 → 3 → 4). Set your local data paths in the `set directories` chunk at the top of each script. HPC scripts assume a SLURM environment and use parallel processing via the `parallel` and `pbmcapply` packages.

---

## Citation

If you use this code, please cite:

> Zohner, C.M., Wu, Z., Mo, L., Crowther, T.W., Fu, Y.H., Renner, S.S., Vitasse, Y. & Rebindaine, D. (2025). Spring temperature remains the dominant driver of leaf-out in a warming world. *Science* (in review). DOI: [to be added upon publication]

---

## Funding

This work was supported by a European Research Council (ERC) Consolidator Grant under the European Union's Horizon Europe research and innovation programme (grant agreement No. 101229851, CHILL-TIME) awarded to C.M.Z.

---

## Contact

Constantin Zohner — constantin.zohner@branch.eco  
BRANCH Institute, Zug, Switzerland  
Institute for Future Initiatives, University of Tokyo, Japan
