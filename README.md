# Carbon at Risk (CaR) — Replication Package

Replication code for Lee et al. 2026.

The default pipeline reproduces all figures in approximately 20 minutes using pre-computed intermediate data (MODIS spatial correlations and DACCS Monte Carlo results) included in this repository. Code to download and process the raw MODIS burned-area data from NASA Earthdata is also included but is not run by default; see Step 0 below.

## Requirements

- **R** (>= 4.1; tested with R 4.5.0, on Mac OS (M2))
- **R packages**: `tidyverse` (2.0.0), `patchwork` (1.3.0), `truncnorm` (1.0-9), `sf` (1.0-21), `terra` (1.8-70), `evd` (2.3-7.1), `ggrepel` (0.9.6), `ggrastr` (1.0.2), `scales` (1.4.0), `httr`, `jsonlite`, `grid`, `MASS`, `parallel`
- **NASA Earthdata account**: Only required if re-running the MODIS spatial correlation pipeline from scratch (Step 0). Not needed for the default run.
- **Python 3** with `pandas` and `requests`: Only required for the optional driver audit (`code/si/si_defor_driver_audit.py`), which is not part of `run_all.R`.

Install all R packages:

```r
install.packages(c("tidyverse", "patchwork", "truncnorm", "sf", "httr",
                    "jsonlite", "terra", "MASS", "evd", "ggrepel",
                    "ggrastr", "scales", "parallel"))
```

## Quick Start

From the **repository root** (the directory containing `run_all.R`):

```bash
Rscript run_all.R
```

This runs the full pipeline end-to-end and writes all figures to `outputs/`.

### Runtime

Measured on an Apple M4 MacBook Pro, on a clean checkout with `outputs/` deleted apart from the
pre-computed MODIS correlation estimates. Two steps account for most of the wall clock; the rest
are seconds each.

| Step | Script | Time |
|------|--------|------|
| Figure 2 (forest CaR, 5,000 draws x 3 regions, plus a 21-point correlation sweep at K = 1, 10, 100) | `figure2.R` | ~5.5 min |
| SI — portfolio CaR against K and rho (54 cells, up to K = 500) | `si_k_rho_convergence.R` | ~4 min |
| Figure 3 (DACCS, regenerating 5,000 SSC draws for two scenarios) | `figure3.R` | ~3 min |
| Figure 4 (portfolio design, grid search plus SI panels) | `figure4.R` | ~20 sec |
| SI — correlation impact | `si_correlation_impact.R` | ~30 sec |
| Everything else (Figure 1, and nine SI scripts) | — | ~1 min combined |
| **Total, DACCS regenerated** | | **13 min (measured, 781 sec)** |
| **Total, DACCS loaded from cache (the default)** | | **~10 min** |

The default is `OVERWRITE_DACCS_FLAG <- FALSE`, which loads the shipped DACCS results rather than
regenerating them; the cached file reproduces the paper exactly, so the shorter run is sufficient
for replication. Set it `TRUE` to regenerate from scratch, as in the measured run above.

Setting `SKIP_MODIS <- FALSE` adds the MODIS download and processing pipeline (Step 0), which
requires a NASA Earthdata account and takes approximately 1-2 hours on first run depending on
network speed. This is not needed for replication, as the pre-computed correlation outputs are
included in the repository.

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
| 2a | SI — VaR illustration | `code/si/si_1_var.R` | None (schematic) | `outputs/si/si_1_var.pdf` |
| 2b | SI — Distribution assumption sensitivity | `code/si/si_distribution_assumption.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_distribution_assumption.pdf` |
| 2c | SI — K-rho convergence | `code/si/si_k_rho_convergence.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_k_rho_convergence.pdf` |
| 2d | SI — Correlation impact on CaR | `code/si/si_correlation_impact.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_correlation_impact.pdf` |
| 2e | SI — Gamma (climate trend) sensitivity | `code/si/si_gamma_sensitivity.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_gamma_sensitivity.pdf` |
| 2f | SI — Regrowth rate sensitivity | `code/si/si_regrowth_sensitivity.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_regrowth_sensitivity.pdf` |
| 2g | SI — Fire history bar charts | `code/si/si_fire_history.R` | EFFIS fire | `outputs/si/si_fire_history.pdf` |
| 2h | SI — Three phases of the CaR curve | `code/si/si_car_phases.R` | EFFIS fire + forest cover, Zang regrowth rates | `outputs/si/si_car_phases.pdf` |
| 2i | SI — Conversion record (descriptive) | `code/si/si_deforestation_history.R` | `data/conversion_rates.csv`, EFFIS fire | `outputs/si/si_deforestation_history.pdf` |
| 2j | SI — Conversion risk (model and results) | `code/si/si_deforestation.R` | `data/conversion_rates.csv`, EFFIS fire + forest cover | `outputs/si/si_deforestation.pdf`, `si_deforestation_results.csv` |
| 2k | SI — Within-DACCS correlation sweep | `code/si/si_rho_daccs_sweep.R` | None (analytic) | `outputs/si/si_rho_daccs_sweep.pdf` |
| 2l | SI — CaR decomposition when Q < mu | `code/si/si_car_negative_gap.R` | None (schematic) | `outputs/si/si_car_negative_gap.pdf` |

Supplementary figure numbers are deliberately not given here: the SI figure order is not pinned by this repository, so a number quoted in the README would go stale. Match on the output filename instead.

Step 0 is skipped by default (`SKIP_MODIS <- TRUE`) because intermediate correlation outputs are included in the repository. Set `SKIP_MODIS <- FALSE` in `run_all.R` to rerun from scratch.

## Configuration Flags

Set at the top of `run_all.R`:

| Flag | Default | Effect |
|------|---------|--------|
| `SKIP_MODIS` | `TRUE` | Skip MODIS download/processing; use pre-computed correlation outputs |
| `OVERWRITE_DACCS_FLAG` | `FALSE` | Skip the DACCS Monte Carlo if `outputs/intermediate/daccs_mc_raw.rds` exists. Set `TRUE` to regenerate; required if `N_SIMULATIONS` in `figure3.R` is ever changed, because the cache is keyed on file existence alone |

## Reproducibility

Results are deterministic. Every simulation call site is seeded individually with
`CAR_SEED = 101` (defined in `code/0_funcs/fire_funcs.R`) immediately before the call, rather
than once at the top of each script. This matters more than it might appear: a single seed at
the top of a script makes each result depend on how much randomness earlier calls happened to
consume, so the same quantity computed in two places, or the same region processed second
rather than first, comes out different. Seeding per call site removes that dependence.

The observable consequence is that quantities computed in more than one place agree exactly.
For example, the California 100-year 95% CaR is 648.8 kg per tonne contracted in Figure 2's
single-project panel, in its delivery-density panel, in its correlation sweep, in
`si_k_rho_convergence.R` at K=1, and in `si_distribution_assumption.R`'s empirical arm. The
fire-only conversion-risk baseline in `si_deforestation.R` reproduces the Figure 2 values to
zero difference across all three regions and all three horizons, despite reaching them through
a different code path.

Forest and DACCS simulations both use N = 5,000 draws. At that size the seed-to-seed standard
deviation of the California 100-year CaR is about 1.8 kg per tonne, so the reported figures are
stable to their final digit.

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
│       ├── si_1_var.R                  # VaR illustration
│       ├── si_car_phases.R             # Three phases of the CaR curve
│       ├── si_car_negative_gap.R       # CaR decomposition when Q < mu
│       ├── si_correlation_impact.R     # Diversification benefit against rho
│       ├── si_deforestation_common.R   # Shared loaders for the conversion analysis
│       ├── si_deforestation_history.R  # Descriptive conversion-rate figure
│       ├── si_deforestation.R          # Conversion-risk model and results
│       ├── si_defor_driver_audit.py    # Driver shares behind the conversion panel (Python; not in run_all.R)
│       ├── si_distribution_assumption.R # Spliced GPD tail sensitivity
│       ├── si_fire_history.R           # Annual burn fractions
│       ├── si_gamma_sensitivity.R      # Sensitivity to the climate burn-rate trend
│       ├── si_k_rho_convergence.R      # Portfolio CaR against K and rho
│       ├── si_regrowth_sensitivity.R   # Sensitivity to the regrowth rate
│       └── si_rho_daccs_sweep.R        # Sensitivity to within-DACCS correlation
├── data/
│   ├── admin_regrowth_with_gpp.gpkg   # Region boundaries and geo-IDs (3 regions only)
│   ├── conversion_rates.csv           # Annual non-fire anthropogenic loss rates (3 regions, 2001-2025)
│   ├── gfw_cache/                     # Cached GFW driver-attribution API responses
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
        ├── daccs_mc_results.csv       # DACCS CaR and CVaR summary table
        └── figure2_numbers.csv        # Every forest number quoted in the paper
```

## Data Sources

| Dataset | Source | Files | Used by | Purpose |
|---------|--------|-------|---------|---------|
| EFFIS fire data | [EFFIS API](https://effis.jrc.ec.europa.eu/) | `data/effis_cache/effis_fire_*.csv` | Figures 2, SI fire history, all forest sensitivity analyses | Annual burned area (hectares) for California, Mato Grosso, Papua, 2002-2023. Downloaded at runtime via `fetch_fire_data()` in `fire_funcs.R`, which calls the EFFIS API and caches responses locally. Pre-cached files are included so no API call is needed on first run. |
| EFFIS forest cover | [EFFIS API](https://effis.jrc.ec.europa.eu/) | `data/effis_cache/effis_forest_*.rds` | Figures 2, all forest analyses | Total forest area (land-cover class 1) per region, used as the denominator to compute annual burn-area fractions. Downloaded via `fetch_forest_indicators()` in `fire_funcs.R` with the same caching logic. Pre-cached. |
| Admin boundaries | [GADM](https://gadm.org/) | `data/admin_regrowth_with_gpp.gpkg` | Figures 2, all forest analyses | GeoPackage with GADM Level 1 admin boundaries and geo-IDs for the three study regions (California, Mato Grosso, Papua; ~800 KB). Used to look up EFFIS geo-IDs and region names. See `code/0_funcs/prepare_gpkg_subset.R` for the extraction script. |
| Regrowth rates | [Zang et al. (2024)](https://doi.org/10.1038/s41597-024-03896-8) | Computed in `code/0_funcs/regrowth_funcs.R` | Figures 2, all forest analyses | Post-fire regrowth rates calibrated from satellite-derived height-recovery equations. California rate adjusted to 2.0%/yr based on local estimates from [Cook-Patton et al. (2020)](https://doi.org/10.1038/s41586-020-2686-x). See Methods in the paper. |
| MODIS MCD64A1 | [NASA Earthdata](https://earthdata.nasa.gov/) | Downloaded in Step 0 | Figure 2 (panel c: spatial correlation) | Monthly 500m burned-area product, 2002-2023. Processed to 1-degree grid cells to estimate pairwise Spearman correlations within California. Pre-computed outputs included in `outputs/intermediate/correlation_results/`; raw download only needed if `SKIP_MODIS = FALSE`. |
| Tree-cover loss (conversion) | [Global Forest Watch](https://www.globalforestwatch.org/) | `data/conversion_rates.csv`, `data/gfw_cache/*.json` | SI conversion-risk analysis | Annual non-fire anthropogenic tree-cover loss as a fraction of 2000 forest area, for the three regions, 2001-2025. Derived from GFW `gadm__tcl__adm1_change` v20260424 at a 30% canopy threshold, restricted to five non-fire anthropogenic drivers using the dominant-driver attribution of Sims et al. (2025), so wildfire is excluded at source and the fire and conversion hazards are disjoint. `code/si/si_defor_driver_audit.py` reports the driver shares behind it. |
| DACCS/SSC parameters | [Alcalde et al. (2018)](https://doi.org/10.1038/s41467-018-04423-1) | Hard-coded in `code/main/figure3/ssc_*.R` | Figure 3 | Geological storage leakage parameters for offshore (high-integrity) and onshore (low-integrity) scenarios, based on the Storage Security Calculator. |

## Key Parameters

| Parameter | Value | Set in |
|-----------|-------|--------|
| Forest MC simulations | 5,000 | `code/0_funcs/fire_funcs.R` |
| DACCS MC simulations | 5,000 | `code/main/figure3/figure3.R` |
| Random seed | 101 (`CAR_SEED`) | `code/0_funcs/fire_funcs.R` |
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

The pipeline also writes CSV summary tables alongside some SI figures (`si_gamma_sensitivity.csv`, `si_k_rho_convergence.csv`, `si_car_empirical_vs_gpd.csv`, `si_regrowth_sensitivity.csv`, `si_deforestation_results.csv`).

Two of these are the machine-readable record of the numbers quoted in the paper, and are the right place to check a value rather than reading it off a figure:

- `outputs/intermediate/figure2_numbers.csv` — every forest CaR, CVaR and diversification figure, by region, horizon and correlation level.
- `outputs/intermediate/daccs_mc_results.csv` — DACCS CaR at the 90th, 95th and 98th percentiles plus CVaR, for both regulatory scenarios at 200, 1,000 and 10,000 years.

To regenerate all intermediate outputs from scratch, delete `outputs/intermediate/` and set `SKIP_MODIS <- FALSE` and `OVERWRITE_DACCS_FLAG <- TRUE` in `run_all.R`.

## License

This code is released under the MIT License. See `LICENSE` for details.
