# SI: How Many Projects to Tame Tail Risk?
# CaR convergence as K increases, for different correlation levels
#
# Run from repo root:
#   Rscript code/si/si_k_rho_convergence.R
#
# Output:
#   - outputs/si/si_k_rho_convergence.pdf
#   - outputs/si/si_k_rho_convergence.csv

library(sf)
library(httr)
library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
select <- dplyr::select

source("code/0_funcs/fire_funcs.R")
source("code/0_funcs/regrowth_funcs.R")

# Configuration ---------------------------------------------------------------
GPKG_PATH <- "data/admin_regrowth_with_gpp.gpkg"
EFFIS_CACHE <- "data/effis_cache"
ESTATE_AREA <- 1000
CLIMATE_RATE <- 0.005
N_SIMULATIONS <- 1000
T_HORIZON <- 100
REGROWTH_RATES <- get_regrowth_rates()

K_VALUES <- c(1, 2, 5, 10, 20, 50, 100, 200, 500)
RHO_VALUES <- c(0, 0.05, 0.1, 0.25, 0.5)

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(42)

# Load data and extract California parameters ---------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)
row <- gdf_gadm %>% filter(NAME_0 == "United States", NAME_1 == "California") %>% slice(1)

geo_id <- paste0(row$country_code, ".", row$ID_1, "_1")
fires_df <- fetch_fire_data(geo_id, cache_dir = EFFIS_CACHE)
forest_json <- fetch_forest_indicators(geo_id, cache_dir = EFFIS_CACHE)
regrowth_rate <- REGROWTH_RATES[["California"]]

empirical_burn_fractions <- calculate_burn_fractions(
  fires_df, forest_json$lc1, project_area = ESTATE_AREA,
  rescale_firesize = FALSE
)

cat(sprintf("California regrowth rate: %.1f%%\n", regrowth_rate * 100))
cat(sprintf("Burn fraction range: %.4f - %.4f\n",
            min(empirical_burn_fractions), max(empirical_burn_fractions)))

# Sweep over K x rho ----------------------------------------------------------
cat("\nRunning K x rho sweep...\n")
sweep_results <- expand.grid(K = K_VALUES, rho = RHO_VALUES) %>%
  as_tibble() %>%
  mutate(car_95 = NA_real_, expected = NA_real_)

for (i in seq_len(nrow(sweep_results))) {
  k <- sweep_results$K[i]
  rho <- sweep_results$rho[i]
  cat(sprintf("  K=%3d, rho=%.2f ... ", k, rho))

  result <- simulate_diversified_car(
    empirical_burn_fractions, regrowth_rate, CLIMATE_RATE,
    years = c(T_HORIZON),
    n_simulations = N_SIMULATIONS,
    n_projects = k,
    estate_area = ESTATE_AREA,
    correlation = rho
  )

  sweep_results$car_95[i] <- result$car_95
  sweep_results$expected[i] <- result$mean_loss
  cat(sprintf("CaR=%.1f, E=%.1f\n", result$car_95, result$mean_loss))
}

# Save results
write.csv(sweep_results, file.path(out_dir, "si_k_rho_convergence.csv"), row.names = FALSE)

# Plot -------------------------------------------------------------------------
CA_COLOR <- "#E67E22"
expected_val <- sweep_results$expected[1]

sweep_results <- sweep_results %>%
  mutate(rho_label = factor(sprintf("rho == %.2f", rho),
                            levels = sprintf("rho == %.2f", RHO_VALUES)))

p <- ggplot(sweep_results, aes(x = K, y = car_95, color = rho_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = expected_val, linetype = "dashed", color = "gray40") +
  annotate("text", x = max(K_VALUES), y = expected_val,
           label = "Expected loss", hjust = 1, vjust = -0.5,
           size = 6, color = "gray40") +
  scale_x_log10(breaks = K_VALUES, labels = K_VALUES) +
  scale_color_manual(
    values = c(CA_COLOR, "#C0392B", "#2C3E50", "#7F8C8D", "#95A5A6"),
    labels = function(x) parse(text = x)
  ) +
  labs(
    title = NULL,
    x = "Number of Projects (K)",
    y = "95th Percentile CaR (ha)",
    color = NULL
  ) +
  theme_classic(base_size = 22) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(out_dir, "si_k_rho_convergence.pdf"), p, width = 8, height = 5)
cat("\nPlot saved to outputs/si/si_k_rho_convergence.pdf\n")
