# Master pipeline for Nature Climate Change submission
# Run from repo root: Rscript run_all.R
#
# Pipeline:
#   Step 0: Spatial correlation (MODIS download + processing + estimation)
#   Step 1: Main figures (figure1, figure2, figure3/DACCS, figure4/portfolio)
#   Step 2: SI figures (VaR, distribution sensitivity, correlation, gamma,
#           regrowth, fire history, K-rho convergence)
#
# Note: Step 0 requires NASA Earthdata credentials and takes hours on first run.
#       Set SKIP_MODIS = TRUE below to skip if outputs already exist.

SKIP_MODIS <- TRUE             # Set to FALSE to rerun MODIS pipeline from scratch
OVERWRITE_DACCS_FLAG <- FALSE  # Set to TRUE to rerun DACCS MC simulations

cat("═══ NCC Pipeline ═══\n\n")

# ── Step 0: Spatial correlation pipeline ──────────────────────────────────

corr_output <- "outputs/intermediate/correlation_results/cell_1.00deg/within_region.rds"

if (!SKIP_MODIS || !file.exists(corr_output)) {
  cat("── Step 0a: Downloading MODIS data ──\n")
  source("code/main/figure2/spatial_correlation/01_download_modis.R")

  cat("\n── Step 0b: Processing burned area ──\n")
  source("code/main/figure2/spatial_correlation/02_process_burned_area.R")

  cat("\n── Step 0c: Estimating correlations ──\n")
  source("code/main/figure2/spatial_correlation/03_estimate_correlations.R")
} else {
  cat("Step 0: Skipping MODIS pipeline (outputs exist)\n\n")
}

# ── Step 1: Main figures ─────────────────────────────────────────────────

cat("── Step 1a: Figure 1 (CaR definition) ──\n")
source("code/main/figure1.R")

cat("\n── Step 1b: Figure 2 (Forest CaR) ──\n")
source("code/main/figure2/figure2.R")

cat("\n── Step 1c: Figure 3 (DACCS CaR) ──\n")
source("code/main/figure3/figure3.R")

cat("\n── Step 1d: Figure 4 (Portfolio design) ──\n")
source("code/main/figure4/figure4.R")

# ── Step 2: SI figures ───────────────────────────────────────────────────

cat("\n── Step 2a: SI VaR illustration ──\n")
source("code/si/si_1_var.R")

cat("\n── Step 2b: SI Distribution assumption ──\n")
source("code/si/si_distribution_assumption.R")

cat("\n── Step 2c: SI K-rho convergence ──\n")
source("code/si/si_k_rho_convergence.R")

cat("\n── Step 2d: SI Correlation impact ──\n")
source("code/si/si_correlation_impact.R")

cat("\n── Step 2e: SI Gamma sensitivity ──\n")
source("code/si/si_gamma_sensitivity.R")

cat("\n── Step 2f: SI Regrowth sensitivity ──\n")
source("code/si/si_regrowth_sensitivity.R")

cat("\n── Step 2g: SI Fire history ──\n")
source("code/si/si_fire_history.R")

cat("\n═══ Pipeline complete ═══\n")
cat("Main figures: outputs/main/\n")
cat("SI figures:   outputs/si/\n")
