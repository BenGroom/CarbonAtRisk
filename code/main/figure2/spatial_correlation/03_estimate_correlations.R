# Estimate spatial correlations in fire risk
#
# Prerequisites: Run 02_process_burned_area.R first
#
# Run from repo root:
#   Rscript code/main/figure2/spatial_correlation/03_estimate_correlations.R
#   Rscript code/main/figure2/spatial_correlation/03_estimate_correlations.R --cell_size 0.5
#
# Output:
#   - outputs/intermediate/correlation_results/cell_{size}/within_region.rds
#   - outputs/intermediate/correlation_results/cell_{size}/between_region.rds
#   - outputs/intermediate/correlation_results/cell_{size}/correlation_estimates_{size}.csv

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(httr)
library(jsonlite)

source("code/main/figure2/spatial_correlation/config.R")
source("code/main/figure2/spatial_correlation/helpers.R")
source("code/0_funcs/fire_funcs.R")

ensure_dirs()

output_dir <- get_cell_path("correlation_results")
diag_dir <- get_cell_path("diagnostics")

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: Within-Region Correlations (MODIS cell-level)
# ══════════════════════════════════════════════════════════════════════════════

cat("═══ PART 1: Within-Region Correlations ═══\n\n")

within_results <- list()
cell_burn_dir <- get_cell_path("cell_burn_fractions")

for (reg in REGIONS) {
  region <- reg$name
  rds_path <- file.path(cell_burn_dir, sprintf("%s.rds", region))

  if (!file.exists(rds_path)) {
    cat(sprintf("  %s: RDS not found, skipping\n", region))
    next
  }

  cat(sprintf("\n─── %s ───\n", region))

  cell_burns <- readRDS(rds_path)
  centroids <- attr(cell_burns, "centroids")
  n_cells <- attr(cell_burns, "n_cells")

  burn_wide <- cell_burns %>%
    select(cell_id, year, burn_fraction) %>%
    pivot_wider(names_from = cell_id, values_from = burn_fraction,
                names_prefix = "cell_")

  burn_mat <- as.matrix(burn_wide[, -1])

  valid <- apply(burn_mat, 2, function(x) {
    !all(is.na(x)) && sd(x, na.rm = TRUE) > 0
  })
  burn_mat <- burn_mat[, valid, drop = FALSE]
  valid_idx <- which(valid)

  cat(sprintf("  %d cells with valid data (of %d total)\n", ncol(burn_mat), n_cells))

  if (ncol(burn_mat) < 3) {
    cat("  Too few valid cells, skipping\n")
    next
  }

  cor_mat <- cor(burn_mat, use = "pairwise.complete.obs", method = "spearman")
  upper <- upper.tri(cor_mat)
  pairs_cor <- cor_mat[upper]

  cat(sprintf("  Pairwise Spearman correlations:\n"))
  cat(sprintf("    Median: %.4f\n", median(pairs_cor, na.rm = TRUE)))
  cat(sprintf("    Mean:   %.4f\n", mean(pairs_cor, na.rm = TRUE)))
  cat(sprintf("    IQR:    [%.4f, %.4f]\n",
              quantile(pairs_cor, 0.25, na.rm = TRUE),
              quantile(pairs_cor, 0.75, na.rm = TRUE)))
  cat(sprintf("    N pairs: %d\n", sum(!is.na(pairs_cor))))

  cent_valid <- centroids[valid_idx, , drop = FALSE]
  cent_sf <- st_as_sf(data.frame(lon = cent_valid[, 1], lat = cent_valid[, 2]),
                      coords = c("lon", "lat"), crs = 4326)
  dist_m <- st_distance(cent_sf)
  dist_km <- as.matrix(units::drop_units(dist_m)) / 1000

  pair_idx <- which(upper, arr.ind = TRUE)
  pairs_df <- data.frame(
    cell_i = pair_idx[, 1],
    cell_j = pair_idx[, 2],
    cor_spearman = pairs_cor,
    distance_km = dist_km[upper]
  )

  pairs_clean <- pairs_df %>% filter(!is.na(cor_spearman), is.finite(distance_km))

  decay_fit <- tryCatch({
    nls(cor_spearman ~ a * exp(-distance_km / lambda),
        data = pairs_clean,
        start = list(a = 0.3, lambda = 100),
        control = nls.control(maxiter = 200))
  }, error = function(e) {
    cat(sprintf("    Distance-decay fit failed: %s\n", e$message))
    NULL
  })

  if (!is.null(decay_fit)) {
    a_hat <- coef(decay_fit)["a"]
    lambda_hat <- coef(decay_fit)["lambda"]
    cat(sprintf("  Distance-decay: rho(d) = %.3f * exp(-d / %.0f km)\n",
                a_hat, lambda_hat))
  } else {
    a_hat <- NA
    lambda_hat <- NA
  }

  # Diagnostic plots
  p_decay <- ggplot(pairs_clean, aes(x = distance_km, y = cor_spearman)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = sprintf("Within-Region Correlation vs Distance - %s", region),
         x = "Distance between cell centroids (km)",
         y = "Spearman correlation") +
    theme_minimal()

  if (!is.null(decay_fit)) {
    pred_d <- seq(0, max(pairs_clean$distance_km), length.out = 200)
    pred_rho <- a_hat * exp(-pred_d / lambda_hat)
    p_decay <- p_decay +
      geom_line(data = data.frame(distance_km = pred_d, cor_spearman = pred_rho),
                color = "darkgreen", linewidth = 1) +
      annotate("text", x = max(pred_d) * 0.6, y = a_hat * 0.9,
               label = sprintf("rho = %.3f * exp(-d / %.0f km)", a_hat, lambda_hat),
               color = "darkgreen", size = 4)
  }

  ggsave(file.path(diag_dir, sprintf("correlation_vs_distance_within_%s_%s.pdf",
                                      region, PARAMS$cell_size_label)),
         p_decay, width = 8, height = 6)

  med_cor <- median(pairs_cor, na.rm = TRUE)
  p_hist <- ggplot(pairs_df, aes(x = cor_spearman)) +
    geom_histogram(bins = 40, fill = "steelblue", color = "white") +
    geom_vline(xintercept = med_cor, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = sprintf("Within-Region Pairwise Correlations - %s", region),
         subtitle = sprintf("Median Spearman = %.4f", med_cor),
         x = "Spearman correlation", y = "Count") +
    theme_minimal()

  ggsave(file.path(diag_dir, sprintf("correlation_histogram_within_%s_%s.pdf",
                                      region, PARAMS$cell_size_label)),
         p_hist, width = 7, height = 5)

  within_results[[region]] <- list(
    cor_matrix = cor_mat,
    pairs = pairs_df,
    decay_fit = decay_fit,
    rho_within = median(pairs_cor, na.rm = TRUE),
    rho_within_se = sd(pairs_cor, na.rm = TRUE) / sqrt(sum(!is.na(pairs_cor))),
    lambda_km = lambda_hat,
    n_cells = ncol(burn_mat),
    n_pairs = sum(!is.na(pairs_cor))
  )
}

saveRDS(within_results, file.path(output_dir, "within_region.rds"))

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: Between-Region Correlations (EFFIS ADM1-level)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n═══ PART 2: Between-Region Correlations ═══\n\n")

gdf <- st_read(get_path("gpkg"), quiet = TRUE)

TARGET_COUNTRIES <- c("United States", "Brazil", "Indonesia")

# Cache EFFIS responses in data/effis_cache/
effis_cache <- "data/effis_cache"

fetch_fire_safe <- function(geo_id) {
  tryCatch(fetch_fire_data(geo_id, cache_dir = effis_cache), error = function(e) NULL)
}

fetch_forest_safe <- function(geo_id) {
  tryCatch(fetch_forest_indicators(geo_id, cache_dir = effis_cache), error = function(e) NULL)
}

centroid_distance_km <- function(sf_df) {
  old_s2 <- sf_use_s2()
  on.exit(sf_use_s2(old_s2))
  sf_use_s2(FALSE)

  centroids <- st_centroid(st_make_valid(sf_df$geom))
  dist_m <- st_distance(centroids)
  units::drop_units(dist_m) / 1000
}

between_results <- list()

for (country in TARGET_COUNTRIES) {
  cat(sprintf("\n─── %s ───\n", country))

  regions <- gdf %>% filter(NAME_0 == country)
  cat(sprintf("  %d ADM1 units found\n", nrow(regions)))

  fire_list <- list()

  for (i in seq_len(nrow(regions))) {
    row <- regions[i, ]
    geo_id <- paste0(row$country_code, ".", row$ID_1, "_1")
    label <- row$NAME_1

    cat(sprintf("  [%d/%d] Fetching %s...", i, nrow(regions), label))

    fires_df <- fetch_fire_safe(geo_id)
    forest_json <- fetch_forest_safe(geo_id)

    if (is.null(fires_df) || is.null(forest_json) || nrow(fires_df) == 0) {
      cat(" SKIP\n")
      next
    }

    burn_fracs <- tryCatch({
      calculate_burn_fractions(fires_df, forest_json$lc1, project_area = 100000)
    }, error = function(e) NULL)

    if (is.null(burn_fracs)) {
      cat(" SKIP\n")
      next
    }

    fire_list[[label]] <- data.frame(
      year = fires_df$year,
      burn_fraction = burn_fracs,
      ba_area_ha = fires_df$ba_area_ha
    )

    cat(sprintf(" OK (%d years)\n", nrow(fires_df)))
    Sys.sleep(0.5)
  }

  if (length(fire_list) < 3) {
    cat("  Too few regions with data, skipping country\n")
    next
  }

  burn_df <- bind_rows(fire_list, .id = "region") %>%
    select(region, year, burn_fraction) %>%
    pivot_wider(names_from = region, values_from = burn_fraction)

  region_names <- setdiff(names(burn_df), "year")

  nonzero_counts <- colSums(burn_df[, region_names] > 0, na.rm = TRUE)
  keep <- names(nonzero_counts[nonzero_counts >= PARAMS$min_fire_years])
  cat(sprintf("  %d regions with >= %d fire years\n",
              length(keep), PARAMS$min_fire_years))

  if (length(keep) < 3) {
    cat("  Too few qualifying regions, skipping\n")
    next
  }

  burn_mat <- as.matrix(burn_df[, keep])

  cor_pearson <- cor(burn_mat, use = "pairwise.complete.obs")
  cor_spearman <- cor(burn_mat, use = "pairwise.complete.obs", method = "spearman")

  upper <- upper.tri(cor_pearson)
  pairs_pearson <- cor_pearson[upper]
  pairs_spearman <- cor_spearman[upper]

  cat(sprintf("  Pairwise correlations:\n"))
  cat(sprintf("    Pearson median:  %.3f\n", median(pairs_pearson, na.rm = TRUE)))
  cat(sprintf("    Spearman median: %.3f\n", median(pairs_spearman, na.rm = TRUE)))

  regions_kept <- regions %>% filter(NAME_1 %in% keep)
  regions_kept <- regions_kept[match(keep, regions_kept$NAME_1), ]

  dist_km <- centroid_distance_km(regions_kept)

  pair_idx <- which(upper, arr.ind = TRUE)
  pairs_df <- data.frame(
    region_i = keep[pair_idx[, 1]],
    region_j = keep[pair_idx[, 2]],
    cor_pearson = pairs_pearson,
    cor_spearman = pairs_spearman,
    distance_km = dist_km[upper]
  )

  decay_fit <- tryCatch({
    nls(cor_pearson ~ a * exp(-distance_km / lambda),
        data = pairs_df,
        start = list(a = 0.5, lambda = 500),
        control = nls.control(maxiter = 200))
  }, error = function(e) {
    cat(sprintf("    Distance-decay fit failed: %s\n", e$message))
    NULL
  })

  if (!is.null(decay_fit)) {
    a_hat <- coef(decay_fit)["a"]
    lambda_hat <- coef(decay_fit)["lambda"]
    cat(sprintf("  Distance-decay: rho(d) = %.3f * exp(-d / %.0f km)\n",
                a_hat, lambda_hat))
  } else {
    a_hat <- NA
    lambda_hat <- NA
  }

  p_decay <- ggplot(pairs_df, aes(x = distance_km, y = cor_pearson)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = TRUE, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = sprintf("Between-Region Correlation vs Distance - %s", country),
         x = "Distance between centroids (km)",
         y = "Pearson correlation") +
    theme_minimal()

  if (!is.null(decay_fit)) {
    pred_d <- seq(0, max(pairs_df$distance_km), length.out = 200)
    pred_rho <- a_hat * exp(-pred_d / lambda_hat)
    p_decay <- p_decay +
      geom_line(data = data.frame(distance_km = pred_d, cor_pearson = pred_rho),
                color = "darkgreen", linewidth = 1)
  }

  ggsave(file.path(diag_dir, sprintf("correlation_vs_distance_between_%s_%s.pdf",
                                     gsub(" ", "_", country), PARAMS$cell_size_label)),
         p_decay, width = 8, height = 6)

  between_results[[country]] <- list(
    burn_matrix = burn_mat,
    cor_pearson = cor_pearson,
    cor_spearman = cor_spearman,
    pairs = pairs_df,
    decay_fit = decay_fit,
    rho_between = median(pairs_pearson, na.rm = TRUE),
    lambda_km = lambda_hat,
    n_regions = length(keep),
    n_pairs = nrow(pairs_df)
  )
}

saveRDS(between_results, file.path(output_dir, "between_region.rds"))

# ══════════════════════════════════════════════════════════════════════════════
# Create Summary CSV
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n═══ SUMMARY: Correlation Estimates ═══\n\n")

summary_rows <- list()

for (region in names(within_results)) {
  r <- within_results[[region]]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    region = region,
    cell_size_deg = PARAMS$cell_size_deg,
    rho_within = round(r$rho_within, 4),
    rho_within_se = round(r$rho_within_se, 4),
    lambda_km = if (is.null(r$lambda_km) || is.na(r$lambda_km)) NA else round(r$lambda_km, 0),
    n_cells = r$n_cells,
    n_pairs = r$n_pairs,
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)

csv_path <- file.path(output_dir, sprintf("correlation_estimates_%s.csv", PARAMS$cell_size_label))
write.csv(summary_df, csv_path, row.names = FALSE)

cat("Within-region correlation estimates:\n")
print(summary_df, row.names = FALSE)

cat(sprintf("\nSaved to %s\n", csv_path))

cat("\n═══ Correlation estimation complete ═══\n")
