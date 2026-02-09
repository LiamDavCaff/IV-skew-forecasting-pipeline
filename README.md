# IV Skew Forecasting Pipelines (Monthly, Daily, Cross-Market)

This repository contains the full empirical pipeline used to evaluate whether the implied volatility (IV) skew measures predict future equity index excess returns under a consistent out-of-sample (OOS) forecasting framework.

Where applicable, the OOS framework uses:

-   Expanding/rolling window regressions with no look-ahead bias

-   Forecast performance is gauged against like-for-like expanding/rolling historical mean benchmarks

-   Inference is assessed using Clark-West (CW) test

-   Multiple testing bias corrected using BH-FDR on CW p-values

It includes:

\- **Monthly** forecasting pipeline (panel build → EDA → in-sample → OOS → heatmaps/cumulative $R^2$ → rolling/expanding betas)

\- **Daily** OOS forecasting heatmaps

\- **Cross-market** (multi-index) monthly OOS heatmap + significance tables

\- An **IV surface builder** example (3D surface for a selected date)

> **Important:** Input data sets are created from Bloomberg outputs and are **not included** in the repo. The repo is set up so that **no data, figures, or tables are pushed to GitHub**.

------------------------------------------------------------------------

## Repository structure

-   `scripts/`

    -   `monthly_forecasting_pipeline.R` — builds monthly panel, conducts EDA, in-sample results, runs OOS forecasts, produces heatmaps + path plots + rolling/expanding betas\
    -   `daily_forecasting_pipeline.R` — runs daily OOS forecasting heatmaps\
    -   `cross_market_pipeline.R` — multi-index monthly OOS forecasting heatmap + significance tables\
    -   `spx_surface_builder.R` — IV surface visualisation example (3D surface for a selected date)

-   `data/`\
    Place required `.rds` input files here (see **Data requirements** below).

-   `outputs/`\
    Figures/tables are written here. Outputs are ignored by Git (folders kept via `.gitkeep`).

-   `bloomberg_scripts/`\
    Helper scripts used to build/pull datasets.\
    **Must run these first** to ensure the relevant `.rds` files are created in `data/` before running the pipelines in `scripts/`.

------------------------------------------------------------------------

## Requirements

### Software

-   R (recommended: R 4.2+)
-   RStudio (recommended)
-   Bloomberg terminal/API subscription

### Packages

The scripts load common packages (e.g., `dplyr`, `ggplot2`, `tidyr`, `lubridate`, `sandwich`, `lmtest`, `scico`, etc.).

If you’re missing packages, install them with:

``` r
 install.packages(c( "dplyr","tidyr","ggplot2","lubridate","zoo",`

`"slider","scales","stringr","cowplot","scico","ggtext","lmtest",`

`"sandwich","broom","np","knitr","kableExtra","tibble","purrr",`

`"readr","roll","plot3D", "viridisLite", "here"))
```

### Data requirements

All required inputs must exist as `.rds` files in `data/` to successfully run below scripts.

#### 1) Core monthly/daily data required

S&P 500 data required for `scripts/monthly_forecasting_pipeline.R` and `scripts/daily_forecasting_pipeline.R`

-   `data/bbg_spx_data.rds`

#### 2) Cross market data required

Required for `scripts/cross_market_pipeline.R`. Script expects one `.rds` file per index, with the below naming convention (see **Quick Start** below)

-   `data/bbg_spx_data.rds`

-   `data/bbg_ftse_data.rds`

-   `data/bbg_dax_data.rds`

-   `data/bbg_asx_data.rds`

-   `data/bbg_estox_data.rds`

-   `data/bbg_kospi_data.rds`

-   `data/bbg_nasdaq_data.rds`

-   `data/bbg_nifty_data.rds`

-   `data/bbg_nikkei_data.rds`

-   `data/bbg_smi_data.rds`

#### 3) IV Surface data required

Used by `scripts/spx_surface_builder.R` to visualise the 3D IV surface (IV, moneyness and tenor)

-   `data/spx_iv_surface_long.rds`

## Quick Start

#### 1) Generate data sets (required)

Run the scripts in `bloomberg_scripts/` first to create all the necessary `.rds` files in `data/`

Note, users without subscribing to Bloomberg will not be able to run required scripts.

#### 2) Run pipelines from the repo root

In R/RStudio (from the project root directory)

``` r
source("scripts/monthly_forecasting_pipeline.R")
source("scripts/daily_forecasting_pipeline.R")
source("scripts/cross_market_pipeline.R")
source("scripts/spx_surface_builder.R")  # optional
```

#### 3) Outputs

Once scripts are run, figures and tables will be written to the following folders:

-   `outputs/monthly/figures/` and `outputs/monthly/tables`

-   `outputs/daily/figures/`

-   `outputs/cross_market/figures/` and `outputs/cross_market/appendix/`

Folders exist in the repo via `.gitkeep` , but **generated outputs are ignored by Git.**

## Notes on publication/licensing of data

This repo intentionally excludes Bloomberg data and any derived tables and figures from the proprietary data.

-   **Code only** is published here.

-   **No raw data** is committed.

-   **No generated figures/tables** are committed.

If you are reproducing the results, you must generate your own input `.rds` files from Bloomberg (or another reputable source)
