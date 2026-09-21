# Carbon at Risk (CaR) — Replication Package

Replication code for Lee et al. 2026.

The default pipeline reproduces every figure in the paper and the SI in roughly 15 minutes, using the pre-computed intermediate data (MODIS spatial correlations and DACCS Monte Carlo results) included in this repository. Code to download and process the raw MODIS burned-area data from NASA Earthdata is also included but is not run by default; see Step 0 below.

## Requirements

- **R** >= 4.1.
- **NASA Earthdata account**: only required if re-running the MODIS spatial correlation pipeline from scratch (Step 0). Not needed for the default run.
- **Python** >= 3.9 with `pandas` and `requests`: only required for `code/si/si_defor_driver_audit.py`, which is optional and is not called by `run_all.R`.

Install the R packages:

```r
install.packages(c("tidyverse", "patchwork", "truncnorm", "sf", "terra",
                   "evd", "ggrepel", "ggrastr", "scales", "httr", "jsonlite"))
```

`grid`, `parallel`, `MASS` and `stats` ship with R and do not need installing.

### Environment used

The shipped outputs in `outputs/` were produced with the following versions. Any
reasonably recent combination should reproduce them, but these are the exact ones.

| Package | Version | Package | Version | Package | Version |
|---|---|---|---|---|---|
| `tidyverse` | 2.0.0 | `sf` | 1.0.21 | `httr` | 1.4.7 |
| `ggplot2` | 4.0.3 | `terra` | 1.8.70 | `jsonlite` | 2.0.0 |
| `dplyr` | 1.1.4 | `evd` | 2.3.7.1 | `MASS` | 7.3.65 |
| `tidyr` | 1.3.1 | `ggrepel` | 0.9.6 | `grid` | 4.5.0 |
| `purrr` | 1.1.0 | `ggrastr` | 1.0.2 | `parallel` | 4.5.0 |
| `patchwork` | 1.3.2 | `scales` | 1.4.0 | `truncnorm` | 1.0.9 |

- **R** 4.5.0 (2025-04-11), platform `aarch64-apple-darwin20`
- **`sf` system libraries**: GEOS 3.13.0, GDAL 3.8.5, PROJ 9.5.1
- **Python** 3.10.11, `pandas` 2.3.3, `requests` 2.32.5 (for the optional driver audit only)

`ggplot2` 4.0.3 was built under R 4.5.2 while the run used R 4.5.0. This produces a
harmless startup warning and does not affect results.

## Quick Start

From the **repository root**:

```bash
Rscript run_all.R
```

This runs the full pipeline end-to-end and writes all figures to `outputs/`.

### Runtime

Measured on an Apple M4 MacBook Pro with default settings (`SKIP_MODIS <- TRUE`, `OVERWRITE_DACCS_FLAG <- FALSE`), which use the pre-computed MODIS correlation estimates and cached DACCS simulations included in this repository:

| Step | Script | Time |
|------|--------|------|
| 1a Figure 1 (CaR definition) | `figure1.R` | 2 sec |
| 1b Figure 2 (forest CaR, 3 regions + 21-point correlation sweep at K = 1, 10, 100) | `figure2.R` | **470 sec** |
| 1c Figure 3 (DACCS, loads cached MC draws) | `figure3.R` | 2 sec |
| 1d Figure 4 (portfolio grid search + 6 SI panels) | `figure4.R` | 95 sec |
| 2a SI VaR illustration | `si_1_var.R` | < 1 sec |
| 2b SI distribution assumption (GPD) | `si_distribution_assumption.R` | 1 sec |
| 2c SI K-rho convergence (54 cells, K up to 500) | `si_k_rho_convergence.R` | **382 sec** |
| 2d SI correlation impact | `si_correlation_impact.R` | 26 sec |
| 2e SI gamma sensitivity | `si_gamma_sensitivity.R` | 2 sec |
| 2f SI regrowth sensitivity | `si_regrowth_sensitivity.R` | < 1 sec |
| 2g SI fire history | `si_fire_history.R` | < 1 sec |
| 2h SI CaR phases | `si_car_phases.R` | < 1 sec |
| 2i SI conversion record | `si_deforestation_history.R` | < 1 sec |
| 2j SI conversion risk | `si_deforestation.R` | 5 sec |
| 2k SI DACCS correlation sweep | `si_rho_daccs_sweep.R` | 4 sec |
| 2l SI CaR decomposition when Q < mu | `si_car_negative_gap.R` | < 1 sec |
| **Total** | | **~16 minutes** |

Two steps account for almost all of the runtime: **`figure2.R` (about 8 minutes)** and
**`si_k_rho_convergence.R` (about 6 minutes)**. Both run large Monte Carlo sweeps at
`N_SIMULATIONS = 5000` over many projects (`figure2.R` sweeps 21 correlation values at
K = 1, 10 and 100; `si_k_rho_convergence.R` covers 9 values of K up to 500 by 6
correlations). They print progress as they go, but each individual cell can take tens of
seconds, so long gaps between lines are expected and do not mean the run has hung.

Two timed runs on this machine gave 13 and 16 minutes, so treat these figures as
indicative rather than exact; they move with background load.

Setting `SKIP_MODIS <- FALSE` adds the MODIS download and processing pipeline (Step 0), which requires a NASA Earthdata account and takes approximately 1-2 hours on first run depending on network speed. This is not needed for replication, as the pre-computed correlation outputs are included in the repository.

Setting `OVERWRITE_DACCS_FLAG <- TRUE` reruns the DACCS Monte Carlo simulations from scratch rather than loading cached results (adds roughly 2 minutes).

## Pipeline Steps

| Step | Description | Script(s) | Data inputs | Output |
|------|-------------|-----------|-------------|--------|
| 0a | Download MODIS burned-area tiles | `code/main/figure2/spatial_correlation/01_download_modis.R` | NASA Earthdata (remote) | Raw HDF files |
| 0b | Process burned area to grid cells | `code/main/figure2/spatial_correlation/02_process_burned_area.R` | Step 0a output | Processed rasters |
| 0c | Estimate pairwise spatial correlations | `code/main/figure2/spatial_correlation/03_estimate_correlations.R` | Step 0b output | `outputs/intermediate/correlation_results/` |
| 1a | CaR definition and buffer interpretation | `code/main/figure1.R` | None (schematic) | `outputs/main/figure1.pdf`, `subfigs/figure1{a,b}.pdf` |
| 1b | Forest fire CaR, diversification, spatial correlation | `code/main/figure2/figure2.R` | EFFIS fire + forest cover, Zang regrowth rates, MODIS correlations (Step 0c) | `outputs/main/figure2.pdf`, `subfigs/figure2_{a..e}.pdf`, `outputs/si/si_rho_distance.pdf`, `si_density_k.pdf`, `outputs/intermediate/figure2_numbers.csv` |
| 1c | DACCS geological storage CaR | `code/main/figure3/figure3.R` | SSC parameters (hard-coded from Alcalde et al.) | `outputs/main/figure3.pdf`, `outputs/si/si_daccs_{1000,10000}yr.pdf`, `outputs/intermediate/daccs_mc_{raw.rds,results.csv}` |
| 1d | Portfolio design and effective cost | `code/main/figure4/figure4.R` | Calibrated from Figs 2-3 (survival probabilities, costs) | `outputs/main/figure4.pdf`, six `outputs/si/si_portfolio_*.pdf` |
| 2a | VaR illustration | `code/si/si_1_var.R` | None (schematic) | `outputs/si/si_1_var.pdf` |
| 2b | Fire distribution assumption, empirical vs spliced GPD | `code/si/si_distribution_assumption.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_distribution_assumption.pdf`, `si_car_empirical_vs_gpd.{pdf,csv}` |
| 2c | Portfolio CaR convergence in K and rho | `code/si/si_k_rho_convergence.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_k_rho_convergence.{pdf,csv}` |
| 2d | Effect of inter-project correlation on the diversification benefit | `code/si/si_correlation_impact.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_correlation_impact.pdf` |
| 2e | Sensitivity to the climate trend in the burn rate (gamma) | `code/si/si_gamma_sensitivity.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_gamma_sensitivity.{pdf,csv}` |
| 2f | Sensitivity to the regrowth rate | `code/si/si_regrowth_sensitivity.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_regrowth_sensitivity.{pdf,csv}` |
| 2g | Annual burn fraction record for the three focal regions | `code/si/si_fire_history.R` | EFFIS fire | `outputs/si/si_fire_history.pdf` |
| 2h | The three phases of the CaR curve against the moving equilibrium | `code/si/si_car_phases.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_car_phases.pdf` |
| 2i | Conversion (non-fire) hazard record for the three regions | `code/si/si_deforestation_history.R` | `data/conversion_rates.csv`, EFFIS fire | `outputs/si/si_deforestation_history.pdf` |
| 2j | Forest CaR extended with a conversion hazard | `code/si/si_deforestation.R` | `data/conversion_rates.csv`, EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_deforestation.pdf`, `si_deforestation_results.csv` |
| 2k | Sensitivity of the portfolio result to the within-DACCS correlation | `code/si/si_rho_daccs_sweep.R` | Portfolio calibration (as Fig 4) | `outputs/si/si_rho_daccs_sweep.pdf` |
| 2l | CaR decomposition when contracted volume falls below expected delivery | `code/si/si_car_negative_gap.R` | None (schematic) | `outputs/si/si_car_negative_gap.pdf` |

`code/si/si_deforestation_common.R` is a shared helper sourced by steps 2g, 2h, 2i and 2j. It is not a pipeline step and is not run on its own.

`code/si/si_defor_driver_audit.py` is an optional audit of the GFW driver attribution behind `data/conversion_rates.csv`. It is not called by `run_all.R`; run it with `python3 code/si/si_defor_driver_audit.py` to regenerate `outputs/si/si_defor_driver_audit.csv`.

Step 0 is skipped by default (`SKIP_MODIS <- TRUE`) because intermediate correlation outputs are included in the repository. Set `SKIP_MODIS <- FALSE` in `run_all.R` to rerun from scratch.

## Configuration Flags

Set at the top of `run_all.R`:

| Flag | Default | Effect |
|------|---------|--------|
| `SKIP_MODIS` | `TRUE` | Skip MODIS download/processing; use pre-computed correlation outputs |
| `OVERWRITE_DACCS_FLAG` | `FALSE` | Skip DACCS Monte Carlo if `outputs/intermediate/daccs_mc_raw.rds` exists |

## Directory Structure

```
├── run_all.R                          # Master pipeline script
├── README.md
├── code/
│   ├── 0_funcs/
│   │   ├── fire_funcs.R               # Fire simulation, copula sampling, EFFIS data functions
│   │   ├── regrowth_funcs.R           # Zang et al. (2024) regrowth rate calibration
│   │   ├── portfolio_funcs.R          # Bernoulli portfolio model, cost optimisation
│   │   └── prepare_gpkg_subset.R      # Extract 3-region subset from the full global GeoPackage (transparency only; the full file is not shipped)
│   ├── main/
│   │   ├── figure1.R                  # Figure 1: CaR definition schematic
│   │   ├── figure2/
│   │   │   ├── figure2.R             # Figure 2: Forest CaR (3-panel, plus 2 SI figures)
│   │   │   ├── spatial_correlation/   # MODIS correlation pipeline (Steps 0a-0c)
│   │   │   │   ├── config.R
│   │   │   │   ├── helpers.R
│   │   │   │   ├── 01_download_modis.R
│   │   │   │   ├── 02_process_burned_area.R
│   │   │   │   └── 03_estimate_correlations.R
│   │   ├── figure3/
│   │   │   ├── figure3.R             # Figure 3: DACCS geological storage CaR
│   │   │   ├── ssc_common.R          # Storage Security Calculator core functions
│   │   │   ├── ssc_offshore.R        # Offshore scenario parameters
│   │   │   └── ssc_onshore.R         # Onshore scenario parameters
│   │   └── figure4/
│   │       └── figure4.R             # Figure 4: Portfolio design (plus 6 SI figures)
│   └── si/                            # Supplementary Information figures (steps 2a-2l)
│       ├── si_1_var.R                  # 2a VaR illustration
│       ├── si_distribution_assumption.R  # 2b empirical vs spliced GPD tail
│       ├── si_k_rho_convergence.R     # 2c portfolio CaR against K and rho
│       ├── si_correlation_impact.R    # 2d correlation vs diversification benefit
│       ├── si_gamma_sensitivity.R     # 2e climate trend sensitivity
│       ├── si_regrowth_sensitivity.R  # 2f regrowth rate sensitivity
│       ├── si_fire_history.R          # 2g burn fraction record
│       ├── si_car_phases.R            # 2h CaR curve phases
│       ├── si_deforestation_history.R # 2i conversion hazard record
│       ├── si_deforestation.R         # 2j CaR with conversion hazard
│       ├── si_rho_daccs_sweep.R       # 2k within-DACCS correlation sweep
│       ├── si_car_negative_gap.R      # 2l decomposition when Q < mu
│       ├── si_deforestation_common.R  # shared helper for 2g-2j (not a step)
│       └── si_defor_driver_audit.py   # optional GFW driver audit (not in run_all.R)
├── data/
│   ├── admin_regrowth_with_gpp.gpkg   # Region boundaries and geo-IDs (3 regions only)
│   ├── conversion_rates.csv           # Non-fire conversion hazard panel (shipped input)
│   ├── effis_cache/                   # Cached EFFIS API responses (3 regions, 2002-2023)
│   │   ├── effis_fire_USA_5_1.csv
│   │   ├── effis_fire_BRA_12_1.csv
│   │   ├── effis_fire_IDN_23_1.csv
│   │   ├── effis_forest_USA_5_1.rds
│   │   ├── effis_forest_BRA_12_1.rds
│   │   └── effis_forest_IDN_23_1.rds
│   └── gfw_cache/                     # Cached GFW tree-cover-loss-by-driver responses
│       ├── USA_change.json
│       ├── BRA_change.json
│       └── IDN_change.json
└── outputs/
    ├── main/                          # Main paper figures (PDFs)
    │   ├── figure1.pdf
    │   ├── figure2.pdf
    │   ├── figure3.pdf
    │   ├── figure4.pdf
    │   └── subfigs/                   # Individual panels
    ├── si/                            # SI figures, tables, and CSVs
    └── intermediate/                  # Cached simulation outputs
        ├── correlation_results/       # MODIS spatial correlation estimates
        ├── daccs_mc_raw.rds           # DACCS Monte Carlo raw output
        └── daccs_mc_results.csv       # DACCS CaR summary table
```

## Data Sources

| Dataset | Source | Files | Used by | Purpose |
|---------|--------|-------|---------|---------|
| EFFIS fire data | [EFFIS API](https://effis.jrc.ec.europa.eu/) | `data/effis_cache/effis_fire_*.csv` | Figures 2, SI fire history, all forest sensitivity analyses | Annual burned area (hectares) for California, Mato Grosso, Papua, 2002-2023. Downloaded at runtime via `fetch_fire_data()` in `fire_funcs.R`, which calls the EFFIS API and caches responses locally. Pre-cached files are included so no API call is needed on first run. |
| EFFIS forest cover | [EFFIS API](https://effis.jrc.ec.europa.eu/) | `data/effis_cache/effis_forest_*.rds` | Figures 2, all forest analyses | Total forest area (land-cover class 1) per region, used as the denominator to compute annual burn-area fractions. Downloaded via `fetch_forest_indicators()` in `fire_funcs.R` with the same caching logic. Pre-cached. |
| Admin boundaries | [GADM](https://gadm.org/) | `data/admin_regrowth_with_gpp.gpkg` | Figures 2, all forest analyses | GeoPackage with GADM Level 1 admin boundaries and geo-IDs for the three study regions (California, Mato Grosso, Papua; ~800 KB). Used to look up EFFIS geo-IDs and region names. See `code/0_funcs/prepare_gpkg_subset.R` for the extraction script. |
| Regrowth rates | [Zang et al. (2024)](https://doi.org/10.1038/s41597-024-03896-8) | Computed in `code/0_funcs/regrowth_funcs.R` | Figures 2, all forest analyses | Post-fire regrowth rates calibrated from satellite-derived height-recovery equations. California rate adjusted to 2.0%/yr based on local estimates from [Cook-Patton et al. (2020)](https://doi.org/10.1038/s41586-020-2686-x). See Methods in the paper. |
| MODIS MCD64A1 | [NASA Earthdata](https://earthdata.nasa.gov/) | Downloaded in Step 0 | Figure 2 (panel c: spatial correlation) | Monthly 500m burned-area product, 2002-2023. Processed to 1-degree grid cells to estimate pairwise Spearman correlations within California. Pre-computed outputs included in `outputs/intermediate/correlation_results/`; raw download only needed if `SKIP_MODIS = FALSE`. |
| DACCS/SSC parameters | [Alcalde et al. (2018)](https://doi.org/10.1038/s41467-018-04423-1) | Hard-coded in `code/main/figure3/ssc_*.R` | Figure 3 | Geological storage leakage parameters for offshore (high-integrity) and onshore (low-integrity) scenarios, based on the Storage Security Calculator. |
| Tree-cover loss by driver | [Global Forest Watch](https://www.globalforestwatch.org/) | `data/gfw_cache/{USA,BRA,IDN}_change.json` | SI conversion-risk extension | Annual tree-cover loss in hectares by dominant driver, per admin unit, 2001-2025. GFW dataset `gadm__tcl__adm1_change` v20260424 at a 30% canopy threshold, with the WRI/Google DeepMind 1 km dominant-driver layer of [Sims et al. (2025)](https://doi.org/10.1088/1748-9326/add606). Pre-cached; `code/si/si_defor_driver_audit.py` reports the share of loss falling to each driver. |
| Conversion rates | Derived from the GFW cache above | `data/conversion_rates.csv` | `code/si/si_deforestation.R` | **Shipped input, not built at run time.** See "Deriving `conversion_rates.csv`" below. |

### Deriving `conversion_rates.csv`

This panel is a shipped input: no script in the package rebuilds it, because the
denominator was fetched separately from the GFW API and is not cached. It is
documented here so the values can be checked without re-running that fetch.

Each row is `iso, adm1, year, delta` for the three focal units (USA/5 California,
BRA/12 Mato Grosso, IDN/23 Papua) over 2001-2025, 75 rows. The rate is

    delta[j,y] = (non-fire anthropogenic tree-cover loss in unit j, year y)
                 / (unit j's forest area in 2000)

The **numerator** is fully reproducible from `data/gfw_cache/{USA,BRA,IDN}_change.json`:
sum `loss_ha` over the five non-fire anthropogenic drivers, namely Hard commodities,
Logging, Permanent agriculture, Settlements & Infrastructure, and Shifting cultivation.
Loss attributed to Wildfire, Other natural disturbances, or Unknown is excluded, so
the conversion hazard and the EFFIS fire hazard are disjoint and can be composed
without double-counting.

The **denominator** is the unit's 2000 tree-cover extent at the same 30% canopy
threshold, taken from GFW's extent endpoint at the time the panel was built:

| Unit | 2000 forest area used |
|------|----------------------|
| USA/5 (California)   |  9.843 Mha |
| BRA/12 (Mato Grosso) | 56.979 Mha |
| IDN/23 (Papua)       | 29.568 Mha |

To verify the shipped panel, divide the numerator computed from the cache by the
corresponding figure above; this reproduces `delta` in all 75 rows, to the six
decimal places at which `delta` is stored.

Note the two vintages differ: the loss panel runs to 2025, while the driver layer
(v1.2) covers 2001-2024, so loss in the final year carries the last available
driver classification.

## Key Parameters

| Parameter | Value | Set in |
|-----------|-------|--------|
| Forest MC simulations | 5,000 | `code/main/figure2/figure2.R` (`N_SIMULATIONS`) |
| DACCS MC simulations | 5,000 | `code/main/figure3/figure3.R` (`N_SIMULATIONS`) |
| Climate trend (gamma) | 0.5%/yr | `code/0_funcs/fire_funcs.R` |
| Regrowth rate, California | 2.0%/yr (Zang et al. 2024, adjusted) | `code/0_funcs/regrowth_funcs.R` |
| Regrowth rate, Mato Grosso | 3.0%/yr (Zang et al. 2024) | `code/0_funcs/regrowth_funcs.R` |
| Regrowth rate, Papua | 1.5%/yr (Zang et al. 2024) | `code/0_funcs/regrowth_funcs.R` |
| CaR confidence level | 95% (with 90%, 98% reported) | Throughout |
| Time horizons | 1-200 years (forest), 200-10,000 years (DACCS) | Per script |
| Spatial correlation (California) | rho ~ 0.07 (median Spearman) | Estimated in Step 0c |
| Portfolio project size | 500 tCO2 | `code/main/figure4/figure4.R` |

## Outputs

All generated figures are written to `outputs/`. Main-text figures go to `outputs/main/`, supplementary figures to `outputs/si/`, and intermediate data to `outputs/intermediate/`.

The pipeline also writes CSV summary tables alongside some SI figures (e.g., `si_gamma_sensitivity.csv`, `si_k_rho_convergence.csv`, `si_car_empirical_vs_gpd.csv`).

To regenerate all intermediate outputs from scratch, delete `outputs/intermediate/` and set `SKIP_MODIS <- FALSE` and `OVERWRITE_DACCS_FLAG <- TRUE` in `run_all.R`.

## License

This code is released under the MIT License. See `LICENSE` for details.
