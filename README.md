# Spring leaf-out keeps pace with warming across the Northern Hemisphere

**Zohner, C.M., Wu, Z., Mo, L., Crowther, T.W., Fu, Y.H., Renner, S.S., Vitasse, Y. & Rebindaine, D.**

*Nature* (submitted)

---

## Overview

This repository contains all code used to analyse spring leaf-out phenology across Northern Hemisphere deciduous broadleaf forests. The study integrates ground-based observations from the PEP725 network with satellite-derived phenology from MODIS, and uses Bayesian hierarchical models and growing-degree day simulations to quantify the drivers of spring leaf-out and the mechanisms behind apparent declines in temperature sensitivity.

---

## License

This code is released under the [MIT License](LICENSE). You are free to use, modify, and distribute it with attribution.

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

## System requirements

All R code was developed and tested with **R 4.4** on macOS and Linux. The HPC scripts additionally require a SLURM-based computing cluster. No non-standard hardware is required for local scripts.

Key packages:

| Package | Purpose |
|---------|---------|
| `brms`, `cmdstanr` | Bayesian hierarchical models via Stan |
| `lme4` | Mixed-effects models for moving-window analyses |
| `chillR` | Hourly temperature reconstruction and chilling calculations |
| `data.table`, `tidyverse` | Data wrangling |
| `terra`, `raster`, `sf` | Spatial data handling |
| `geosphere` | Photoperiod calculation |
| `ggplot2`, `patchwork` | Figures |
| `parallel`, `pbmcapply` | Parallel processing (HPC scripts only) |

---

## Installation

**1. Install R packages**

```r
install.packages(c(
  "brms", "lme4", "chillR",
  "data.table", "tidyverse",
  "terra", "raster", "sf", "geosphere",
  "ggplot2", "patchwork",
  "parallel", "pbmcapply"
))
```

**2. Install cmdstanr and CmdStan**

`brms` requires CmdStan as a backend. Install `cmdstanr` from the Stan r-universe, then install CmdStan itself:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()  # downloads and compiles CmdStan (~5–10 min)
```

**Typical install time:** approximately 15–20 minutes on a standard desktop, dominated by CmdStan compilation.

---

## Demo

To verify the installation, run the local (non-HPC) scripts on the provided example data subset:

```r
# From the PEP725/ directory, run scripts in order:
# 1_data_prep.Rmd  →  2_driver_assembly.Rmd  →  3_models.Rmd  →  4_figures.Rmd
```

**Expected output:** Fitted model objects (`.rds`), coefficient summary tables (`.csv`), and figure files (`.pdf`/`.png`) written to the `output/` subdirectory within each branch.

**Expected run time:** Local scripts (data preparation, driver assembly, figures) run in minutes to a few hours on a standard desktop. The Bayesian model fitting scripts and GDD simulations are designed for HPC and require days to weeks of compute time depending on cluster configuration and available cores. Full reproduction of all results therefore requires access to a computing cluster.

---

## Instructions for use

1. Download the required input datasets from the sources listed above and place them in a local data directory.
2. Open each script and set your local paths in the `set directories` chunk at the top of the file.
3. Run scripts in numerical order within each branch (`1_` → `2_` → `3_` → `4_`).
4. HPC scripts (`.r`) should be submitted via SLURM. Local scripts (`.Rmd`, `.qmd`) can be run interactively in RStudio or rendered via `rmarkdown::render()` / `quarto::quarto_render()`.

---

## Citation

If you use this code, please cite:

> Zohner, C.M., Wu, Z., Mo, L., Crowther, T.W., Fu, Y.H., Renner, S.S., Vitasse, Y. & Rebindaine, D. (2025). Spring leaf-out keeps pace with warming across the Northern Hemisphere. *Nature* (in review). DOI: [to be added upon publication]

---

## Funding

This work was supported by a European Research Council (ERC) Consolidator Grant under the European Union's Horizon Europe research and innovation programme (grant agreement No. 101229851, CHILL-TIME) awarded to C.M.Z.

---

## Contact

Constantin Zohner — constantin.zohner@branch.eco  
BRANCH Institute, Zug, Switzerland  
Institute for Future Initiatives, University of Tokyo, Japan
