# Carbon at Risk (CaR) — Replication Package

Replication code for Lee et al. 2026.

The default pipeline reproduces all figures in approximately 2 minutes using pre-computed intermediate data (MODIS spatial correlations and DACCS Monte Carlo results) included in this repository. Code to download and process the raw MODIS burned-area data from NASA Earthdata is also included but is not run by default; see Step 0 below.

## Requirements

- **R** (>= 4.1; tested with R 4.5.0)
- **R packages**: `tidyverse` (2.0.0), `patchwork` (1.3.0), `truncnorm` (1.0-9), `sf` (1.0-21), `terra` (1.8-70), `evd` (2.3-7.1), `ggrepel` (0.9.6), `scales` (1.4.0), `httr`, `jsonlite`, `grid`, `MASS`, `parallel`
- **NASA Earthdata account**: Only required if re-running the MODIS spatial correlation pipeline from scratch (Step 0). Not needed for the default run.

Install all R packages:

```r
install.packages(c("tidyverse", "patchwork", "truncnorm", "sf", "httr",
                    "jsonlite", "terra", "MASS", "evd", "ggrepel",
                    "scales", "parallel"))
```

## Quick Start

From the **repository root** (not from inside ``):

```bash
Rscript run_all.R
```

This runs the full pipeline end-to-end and writes all figures to `outputs/`.

### Runtime

Measured on an Apple M4 MacBook Pro. With default settings (`SKIP_MODIS <- TRUE`, `OVERWRITE_DACCS_FLAG <- FALSE`), which use the pre-computed MODIS correlation estimates and cached DACCS simulations included in this repository:

| Step | Script | Time |
|------|--------|------|
| Figure 1 (CaR definition) | `figure1.R` | 2 sec |
| Figure 2 (Forest CaR, 1,000 MC sims x 3 regions + correlation sweep) | `figure2.R` | 53 sec |
| Figure 3 (DACCS, loads cached 10,000 MC sims) | `figure3.R` | 1 sec |
| Figure 4 (Portfolio design, grid search) | `figure4.R` | 16 sec |
| SI — VaR illustration | `si_1_var.R` | < 1 sec |
| SI — Distribution assumption (GPD) | `si_distribution_assumption.R` | 2 sec |
| SI — K-rho convergence (45 cells) | `si_k_rho_convergence.R` | 40 sec |
| SI — Correlation impact | `si_correlation_impact.R` | 5 sec |
| SI — Gamma sensitivity | `si_gamma_sensitivity.R` | 1 sec |
| SI — Regrowth sensitivity | `si_regrowth_sensitivity.R` | 1 sec |
| SI — Fire history | `si_fire_history.R` | 1 sec |
| **Total** | | **~2 minutes** |

Setting `SKIP_MODIS <- FALSE` adds the MODIS download and processing pipeline (Step 0), which requires a NASA Earthdata account and takes approximately 1-2 hours on first run depending on network speed. This is not needed for replication, as the pre-computed correlation outputs are included in the repository.

Setting `OVERWRITE_DACCS_FLAG <- TRUE` reruns the DACCS Monte Carlo simulations from scratch rather than loading cached results (adds ~2 minutes).

## Pipeline Steps

| Step | Description | Script(s) | Data inputs | Output |
|------|-------------|-----------|-------------|--------|
| 0a | Download MODIS burned-area tiles | `code/main/figure2/spatial_correlation/01_download_modis.R` | NASA Earthdata (remote) | Raw HDF files |
| 0b | Process burned area to grid cells | `code/main/figure2/spatial_correlation/02_process_burned_area.R` | Step 0a output | Processed rasters |
| 0c | Estimate pairwise spatial correlations | `code/main/figure2/spatial_correlation/03_estimate_correlations.R` | Step 0b output | `outputs/intermediate/correlation_results/` |
| 1a | Figure 1 — CaR definition & buffer interpretation | `code/main/figure1.R` | None (schematic) | `outputs/main/figure1.pdf` |
| 1b | Figure 2 — Forest fire CaR, diversification & spatial correlation | `code/main/figure2/figure2.R` | EFFIS fire + forest cover, Zang regrowth rates, MODIS correlations (Step 0c) | `outputs/main/figure2.pdf` |
| 1c | Figure 3 — DACCS/BECCS geological storage CaR | `code/main/figure3/figure3.R` | SSC parameters (hard-coded from Alcalde et al.) | `outputs/main/figure3.pdf` |
| 1d | Figure 4 — Portfolio design & effective cost | `code/main/figure4/figure4.R` | Calibrated from Figs 2-3 outputs (survival probabilities, costs) | `outputs/main/figure4.pdf` |
| 2a | SI — VaR illustration (Fig S1) | `code/si/si_1_var.R` | None (schematic) | `outputs/si/si_1_var.pdf` |
| 2b | SI — Distribution assumption sensitivity (Fig S3) | `code/si/si_distribution_assumption.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_distribution_assumption.pdf` |
| 2c | SI — K-rho convergence (Fig S6) | `code/si/si_k_rho_convergence.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_k_rho_convergence.pdf` |
| 2d | SI — Correlation impact on CaR (Fig S7) | `code/si/si_correlation_impact.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_correlation_impact.pdf` |
| 2e | SI — Gamma (climate trend) sensitivity (Fig S5) | `code/si/si_gamma_sensitivity.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_gamma_sensitivity.pdf` |
| 2f | SI — Regrowth rate sensitivity (Fig S8) | `code/si/si_regrowth_sensitivity.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_regrowth_sensitivity.pdf` |
| 2g | SI — Fire history bar charts (Fig S4) | `code/si/si_fire_history.R` | EFFIS fire | `outputs/si/si_fire_history.pdf` |

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
│   │   └── prepare_gpkg_subset.R      # Extract 3-region subset from full GeoPackage (transparency only)
│   ├── main/
│   │   ├── figure1.R                  # Figure 1: CaR definition schematic
│   │   ├── figure2/
│   │   │   ├── figure2.R             # Figure 2: Forest CaR (composite 5-panel)
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
│   │       └── figure4.R             # Figure 4: Portfolio design (6-panel)
│   └── si/                            # Supplementary Information figures
│       ├── si_1_var.R
│       ├── si_correlation_impact.R
│       ├── si_distribution_assumption.R
│       ├── si_fire_history.R
│       ├── si_gamma_sensitivity.R
│       ├── si_k_rho_convergence.R
│       └── si_regrowth_sensitivity.R
├── data/
│   ├── admin_regrowth_with_gpp.gpkg   # Region boundaries and geo-IDs (3 regions only)
│   └── effis_cache/                   # Cached EFFIS API responses (3 regions, 2002-2023)
│       ├── effis_fire_USA_5_1.csv
│       ├── effis_fire_BRA_12_1.csv
│       ├── effis_fire_IDN_23_1.csv
│       ├── effis_forest_USA_5_1.rds
│       ├── effis_forest_BRA_12_1.rds
│       └── effis_forest_IDN_23_1.rds
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

## Key Parameters

| Parameter | Value | Set in |
|-----------|-------|--------|
| Forest MC simulations | 1,000 | `code/0_funcs/fire_funcs.R` |
| DACCS MC simulations | 10,000 | `code/main/figure3/figure3.R` |
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
