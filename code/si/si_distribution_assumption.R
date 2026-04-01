# SI: Distribution Assumption Comparison — Empirical Resampling vs GPD
# Compares bootstrap approach with fitted fat-tailed (GPD) model for California
#
# Run from repo root:
#   Rscript code/si/si_distribution_assumption.R
#
# Output:
#   - outputs/si/si_distribution_assumption.pdf (2-panel comparison)
#   - outputs/si/si_car_empirical_vs_gpd.pdf (CaR curves comparison)
#   - outputs/si/si_car_empirical_vs_gpd.csv (numeric results)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(evd)
library(sf)
library(httr)
library(jsonlite)

source("code/0_funcs/fire_funcs.R")
source("code/0_funcs/regrowth_funcs.R")

set.seed(42)

# Configuration ---------------------------------------------------------------
GPKG_PATH <- "data/admin_regrowth_with_gpp.gpkg"
EFFIS_CACHE <- "data/effis_cache"
CLIMATE_RATE <- 0.005
N_SIMULATIONS <- 5000
N_DRAWS <- 100000
REGROWTH_RATES <- get_regrowth_rates()

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load California burn fractions -----------------------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)
ca_regions <- list(list(country = "United States", subcountry = "California"))
results <- process_selected_geometries(
  gdf_gadm, ca_regions,
  regrowth_rates = REGROWTH_RATES,
  rescale_firesize = FALSE,
  cache_dir = EFFIS_CACHE
)

ca_burn <- results$empirical_burn_fractions %>%
  filter(grepl("California", geo_label)) %>%
  pull(burn_fraction)

cat(sprintf("California: %d observations, range [%.5f, %.5f], mean=%.5f\n",
            length(ca_burn), min(ca_burn), max(ca_burn), mean(ca_burn)))

# Fit GPD to exceedances above threshold ---------------------------------------
threshold <- median(ca_burn)
gpd_fit <- fpot(ca_burn, threshold = threshold, model = "gpd")
gpd_scale <- gpd_fit$estimate["scale"]
gpd_shape <- gpd_fit$estimate["shape"]

cat(sprintf("GPD fit: threshold=%.5f, scale=%.5f, shape=%.4f\n",
            threshold, gpd_scale, gpd_shape))

# Profile likelihood confidence intervals on GPD parameters
gpd_se <- gpd_fit$std.err
if (!is.null(gpd_se)) {
  cat(sprintf("GPD shape SE: %.4f, 95%% CI: [%.4f, %.4f]\n",
              gpd_se["shape"],
              gpd_shape - 1.96 * gpd_se["shape"],
              gpd_shape + 1.96 * gpd_se["shape"]))
  cat(sprintf("GPD scale SE: %.5f, 95%% CI: [%.5f, %.5f]\n",
              gpd_se["scale"],
              gpd_scale - 1.96 * gpd_se["scale"],
              gpd_scale + 1.96 * gpd_se["scale"]))
}

# Goodness-of-fit: Kolmogorov-Smirnov test on exceedances
exceedances <- ca_burn[ca_burn > threshold] - threshold
ks_result <- ks.test(exceedances, "pgpd", loc = 0, scale = gpd_scale, shape = gpd_shape)
cat(sprintf("KS test on GPD fit: D=%.4f, p-value=%.4f\n", ks_result$statistic, ks_result$p.value))

# Spliced sampling function
sample_spliced <- function(n, burn_data, threshold, scale, shape) {
  below <- burn_data[burn_data <= threshold]
  prob_below <- length(below) / length(burn_data)

  draws <- numeric(n)
  is_below <- runif(n) < prob_below

  draws[is_below] <- sample(below, sum(is_below), replace = TRUE)

  n_above <- sum(!is_below)
  gpd_draws <- rgpd(n_above, loc = 0, scale = scale, shape = shape)
  draws[!is_below] <- threshold + gpd_draws

  pmin(pmax(draws, 0), 1)
}

# Panel 1: GPD fit overlaid on histogram ---------------------------------------
x_grid <- seq(min(ca_burn) * 0.5, max(ca_burn) * 2, length.out = 500)
prob_exceed <- mean(ca_burn > threshold)
gpd_density <- ifelse(
  x_grid > threshold,
  prob_exceed * dgpd(x_grid - threshold, loc = 0, scale = gpd_scale, shape = gpd_shape),
  NA
)

gpd_overlay_df <- data.frame(x = x_grid, density = gpd_density) %>% filter(!is.na(density))

p2 <- ggplot(data.frame(burn = ca_burn), aes(x = burn)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 fill = "steelblue", color = "white", alpha = 0.5) +
  geom_line(data = gpd_overlay_df, aes(x = x, y = density),
            color = "red3", linewidth = 1) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "gray40") +
  annotate("text", x = threshold, y = Inf, label = "Threshold",
           hjust = -0.1, vjust = 2, size = 5, color = "gray40") +
  annotate("text", x = max(ca_burn) * 1.3, y = max(gpd_overlay_df$density) * 0.7,
           label = sprintf("xi = %.3f\nsigma = %.4f", gpd_shape, gpd_scale),
           size = 5, color = "red3", hjust = 0) +
  labs(
    title = "GPD tail fit (red) overlaid on historic data",
    x = "Annual burn fraction",
    y = "Density"
  ) +
  theme_classic(base_size = 16)

# Panel 2: Large-sample draws from both methods --------------------------------
empirical_draws <- sample(ca_burn, N_DRAWS, replace = TRUE)
gpd_draws <- sample_spliced(N_DRAWS, ca_burn, threshold, gpd_scale, gpd_shape)

draws_df <- bind_rows(
  data.frame(burn = empirical_draws, method = "Empirical resampling"),
  data.frame(burn = gpd_draws, method = "Spliced GPD")
)

p3 <- ggplot(draws_df, aes(x = burn, fill = method, color = method)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = c("steelblue", "red3")) +
  scale_color_manual(values = c("steelblue", "red3")) +
  labs(
    title = sprintf("Density of %s draws from each method", format(N_DRAWS, big.mark = ",")),
    x = "Annual burn fraction",
    y = "Density",
    fill = NULL, color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.9))

# Assemble and save ------------------------------------------------------------
fig <- (p2 | p3) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(file.path(out_dir, "si_distribution_assumption.pdf"), fig,
       width = 14, height = 5)

cat("\nSaved to outputs/si/si_distribution_assumption.pdf\n")

# ============================================================================
# CaR Comparison: Empirical Resampling vs Spliced GPD
# ============================================================================

ca_regrowth <- REGROWTH_RATES[["California"]]
cat(sprintf("\nCalifornia regrowth rate: %.3f\n", ca_regrowth))

N_SYNTHETIC <- 10000
gpd_synthetic <- sample_spliced(N_SYNTHETIC, ca_burn, threshold, gpd_scale, gpd_shape)

car_horizons <- seq(1, 200, by = 1)

cat("Running CaR simulation: Empirical resampling...\n")
car_empirical <- run_car_simulation(
  ca_burn, ca_regrowth,
  time_horizons = car_horizons,
  n_simulations = N_SIMULATIONS,
  climate_rate = CLIMATE_RATE
)
car_empirical$method <- "Empirical resampling"

cat("Running CaR simulation: Spliced GPD...\n")
car_gpd <- run_car_simulation(
  gpd_synthetic, ca_regrowth,
  time_horizons = car_horizons,
  n_simulations = N_SIMULATIONS,
  climate_rate = CLIMATE_RATE
)
car_gpd$method <- "Spliced GPD"

car_comparison <- bind_rows(car_empirical, car_gpd) %>%
  mutate(
    car_95_kg = car_95 * 1000,
    mean_car_kg = mean_car * 1000
  )

cat("\n95% CaR at selected horizons (kg/tonne):\n")
for (t in c(1, 10, 50, 100, 200)) {
  emp_val <- car_comparison %>% filter(method == "Empirical resampling", time_horizon == t) %>% pull(car_95_kg)
  gpd_val <- car_comparison %>% filter(method == "Spliced GPD", time_horizon == t) %>% pull(car_95_kg)
  cat(sprintf("  T=%3d: Empirical=%.1f, GPD=%.1f, ratio=%.2f\n",
              t, emp_val, gpd_val, gpd_val / emp_val))
}

p_car <- ggplot(car_comparison, aes(x = time_horizon, y = car_95_kg,
                                     color = method, linetype = method)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Empirical resampling" = "steelblue",
                                "Spliced GPD" = "red3")) +
  scale_linetype_manual(values = c("Empirical resampling" = "solid",
                                   "Spliced GPD" = "dashed")) +
  labs(
    x = "Time Horizon (years)",
    y = expression("95% CaR (kg/tCO"[2]*"e)"),
    color = NULL, linetype = NULL
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.65, 0.3),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

ggsave(file.path(out_dir, "si_car_empirical_vs_gpd.pdf"), p_car,
       width = 8, height = 6)

write.csv(car_comparison, file.path(out_dir, "si_car_empirical_vs_gpd.csv"), row.names = FALSE)

cat("\nSaved CaR comparison to outputs/si/\n")
