# SI: Sensitivity of CaR to Climate Trend Parameter (gamma)
# Compares CaR curves for gamma in {0, 0.5%, 1%} across three regions
#
# Run from repo root:
#   Rscript code/si/si_gamma_sensitivity.R
#
# Output:
#   - outputs/si/si_gamma_sensitivity.pdf
#   - outputs/si/si_gamma_sensitivity.csv

library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(httr)
library(jsonlite)

source("code/0_funcs/fire_funcs.R")
source("code/0_funcs/regrowth_funcs.R")

set.seed(42)

# Configuration ---------------------------------------------------------------
GPKG_PATH <- "data/admin_regrowth_with_gpp.gpkg"
EFFIS_CACHE <- "data/effis_cache"
N_SIMULATIONS <- 1000
ESTATE_AREA <- 1000
REGROWTH_RATES <- get_regrowth_rates()

GAMMA_VALUES <- c(0, 0.005, 0.01)
GAMMA_LABELS <- c("0%", "0.5%", "1.0%")

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Define regions ---------------------------------------------------------------
regions <- list(
  list(country = "United States", subcountry = "California"),
  list(country = "Brazil", subcountry = "Mato Grosso"),
  list(country = "Indonesia", subcountry = "Papua")
)

# Load data --------------------------------------------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)

car_horizons <- seq(1, 200, by = 1)

# Run simulations for each region x gamma combination --------------------------
all_results <- data.frame()

for (reg in regions) {
  row <- gdf_gadm %>%
    filter(NAME_0 == reg$country, NAME_1 == reg$subcountry) %>%
    slice(1)

  geo_id <- paste0(row$country_code, ".", row$ID_1, "_1")
  fires_df <- fetch_fire_data(geo_id, cache_dir = EFFIS_CACHE)
  forest_json <- fetch_forest_indicators(geo_id, cache_dir = EFFIS_CACHE)
  regrowth_rate <- REGROWTH_RATES[[reg$subcountry]]

  burn_fracs <- calculate_burn_fractions(
    fires_df, forest_json$lc1, project_area = ESTATE_AREA,
    rescale_firesize = FALSE
  )

  region_label <- paste0(reg$country, " - ", reg$subcountry)
  cat(sprintf("\n%s (regrowth=%.1f%%)\n", region_label, regrowth_rate * 100))

  for (gi in seq_along(GAMMA_VALUES)) {
    gamma <- GAMMA_VALUES[gi]
    cat(sprintf("  gamma = %s ... ", GAMMA_LABELS[gi]))

    car_result <- run_car_simulation(
      burn_fracs, regrowth_rate,
      time_horizons = car_horizons,
      n_simulations = N_SIMULATIONS,
      climate_rate = gamma
    )

    car_result$region <- region_label
    car_result$gamma <- gamma
    car_result$gamma_label <- GAMMA_LABELS[gi]

    all_results <- bind_rows(all_results, car_result)
    cat(sprintf("CaR_95 at 200yr = %.1f%%\n", car_result$car_95[car_result$time_horizon == 200] * 100))
  }
}

# Save results -----------------------------------------------------------------
write.csv(all_results, file.path(out_dir, "si_gamma_sensitivity.csv"), row.names = FALSE)

# Plot -------------------------------------------------------------------------
plot_df <- all_results %>%
  mutate(car_95_pct = car_95 * 100)

p <- ggplot(plot_df, aes(x = time_horizon, y = car_95_pct,
                          color = gamma_label, linetype = gamma_label)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ region, scales = "free_y") +
  scale_color_manual(
    values = c("0%" = "steelblue", "0.5%" = "#E67E22", "1.0%" = "red3")
  ) +
  scale_linetype_manual(
    values = c("0%" = "dotted", "0.5%" = "solid", "1.0%" = "dashed")
  ) +
  labs(
    x = "Time Horizon (years)",
    y = expression("95% CaR (% of initial carbon)"),
    color = expression(gamma ~ "(annual trend)"),
    linetype = expression(gamma ~ "(annual trend)")
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave(file.path(out_dir, "si_gamma_sensitivity.pdf"), p, width = 14, height = 5)

cat("\nSaved to outputs/si/si_gamma_sensitivity.pdf\n")

# Summary table ----------------------------------------------------------------
summary_table <- all_results %>%
  filter(time_horizon %in% c(50, 100, 200)) %>%
  select(region, gamma_label, time_horizon, car_95) %>%
  mutate(car_95_pct = sprintf("%.1f%%", car_95 * 100)) %>%
  select(-car_95) %>%
  pivot_wider(names_from = time_horizon, values_from = car_95_pct,
              names_prefix = "T=")

cat("\nSummary:\n")
print(summary_table)
