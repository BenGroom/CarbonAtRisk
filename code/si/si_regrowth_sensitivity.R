# SI: Sensitivity of CaR to Regrowth Rate Assumptions
# Compares CaR curves when regrowth is scaled by {0.5, 1.0, 1.5} for California
#
# Run from repo root:
#   Rscript code/si/si_regrowth_sensitivity.R
#
# Output:
#   - outputs/si/si_regrowth_sensitivity.pdf
#   - outputs/si/si_regrowth_sensitivity.csv

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
CLIMATE_RATE <- 0.005
N_SIMULATIONS <- 1000
ESTATE_AREA <- 1000
REGROWTH_RATES <- get_regrowth_rates()

REGROWTH_MULTIPLIERS <- c(0.5, 1.0, 1.5)
REGROWTH_LABELS <- c("0.5x baseline", "1.0x baseline", "1.5x baseline")

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load California data ---------------------------------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)
row <- gdf_gadm %>%
  filter(NAME_0 == "United States", NAME_1 == "California") %>%
  slice(1)

geo_id <- paste0(row$country_code, ".", row$ID_1, "_1")
fires_df <- fetch_fire_data(geo_id, cache_dir = EFFIS_CACHE)
forest_json <- fetch_forest_indicators(geo_id, cache_dir = EFFIS_CACHE)
base_regrowth <- REGROWTH_RATES[["California"]]

burn_fracs <- calculate_burn_fractions(
  fires_df, forest_json$lc1, project_area = ESTATE_AREA,
  rescale_firesize = FALSE
)

cat(sprintf("California baseline regrowth rate: %.1f%%/yr\n", base_regrowth * 100))

# Run simulations for each regrowth multiplier ---------------------------------
car_horizons <- seq(1, 200, by = 1)
all_results <- data.frame()

for (mi in seq_along(REGROWTH_MULTIPLIERS)) {
  mult <- REGROWTH_MULTIPLIERS[mi]
  regrowth <- min(base_regrowth * mult, 0.1)
  cat(sprintf("  Regrowth = %s (%.1f%%/yr) ... ", REGROWTH_LABELS[mi], regrowth * 100))

  car_result <- run_car_simulation(
    burn_fracs, regrowth,
    time_horizons = car_horizons,
    n_simulations = N_SIMULATIONS,
    climate_rate = CLIMATE_RATE
  )

  car_result$multiplier <- mult
  car_result$label <- REGROWTH_LABELS[mi]
  car_result$regrowth_pct <- regrowth * 100

  all_results <- bind_rows(all_results, car_result)
  cat(sprintf("CaR_95 at 200yr = %.1f%%\n", car_result$car_95[car_result$time_horizon == 200] * 100))
}

# Save results -----------------------------------------------------------------
write.csv(all_results, file.path(out_dir, "si_regrowth_sensitivity.csv"), row.names = FALSE)

# Plot -------------------------------------------------------------------------
plot_df <- all_results %>%
  mutate(car_95_pct = car_95 * 100)

p <- ggplot(plot_df, aes(x = time_horizon, y = car_95_pct,
                          color = label, linetype = label)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("0.5x baseline" = "red3", "1.0x baseline" = "#E67E22", "1.5x baseline" = "steelblue")
  ) +
  scale_linetype_manual(
    values = c("0.5x baseline" = "dashed", "1.0x baseline" = "solid", "1.5x baseline" = "dotted")
  ) +
  labs(
    title = "California: Sensitivity to regrowth rate",
    x = "Time Horizon (years)",
    y = expression("95% CaR (% of initial carbon)"),
    color = "Regrowth rate",
    linetype = "Regrowth rate"
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.3),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

ggsave(file.path(out_dir, "si_regrowth_sensitivity.pdf"), p, width = 8, height = 6)

cat("\nSaved to outputs/si/si_regrowth_sensitivity.pdf\n")

# Summary table ----------------------------------------------------------------
summary_table <- all_results %>%
  filter(time_horizon %in% c(50, 100, 200)) %>%
  select(label, regrowth_pct, time_horizon, car_95) %>%
  mutate(car_95_pct = sprintf("%.1f%%", car_95 * 100)) %>%
  select(-car_95) %>%
  pivot_wider(names_from = time_horizon, values_from = car_95_pct,
              names_prefix = "T=")

cat("\nSummary:\n")
print(summary_table)
