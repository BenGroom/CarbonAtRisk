# SI: Impact of Correlation on Diversification Benefit
# Shows CaR trajectories for K=1 vs K=100 at rho = 0, 0.5, 1
#
# Run from repo root:
#   Rscript code/si/si_correlation_impact.R
#
# Output:
#   - outputs/si/si_correlation_impact.pdf

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
N_SIMULATIONS <- 1000
REGROWTH_RATES <- get_regrowth_rates()

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(42)

# Load data -------------------------------------------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)

# Run diversification at three correlation levels
cat("Running diversification: rho = 0 ...\n")
div_0 <- run_diversification_analysis(gdf_gadm, regrowth_rates = REGROWTH_RATES,
                                       correlation = 0,
                                       n_simulations = N_SIMULATIONS,
                                       cache_dir = EFFIS_CACHE)
cat("Running diversification: rho = 0.5 ...\n")
div_05 <- run_diversification_analysis(gdf_gadm, regrowth_rates = REGROWTH_RATES,
                                        correlation = 0.5,
                                        n_simulations = N_SIMULATIONS,
                                        cache_dir = EFFIS_CACHE)
cat("Running diversification: rho = 1 ...\n")
div_1 <- run_diversification_analysis(gdf_gadm, regrowth_rates = REGROWTH_RATES,
                                       correlation = 1,
                                       n_simulations = N_SIMULATIONS,
                                       cache_dir = EFFIS_CACHE)

# Build comparison data frame -------------------------------------------------
# The ncc version of run_diversification_analysis uses
# time_horizons = c(1, seq(10, max_years, by=10)) internally,
# so results vectors are shorter than div$years (which is 0:max_years).
time_horizons <- c(1, seq(10, 200, by = 10))

correlation_comparison_df <- bind_rows(
  data.frame(
    year = time_horizons,
    mean_loss = div_0$results[["1"]]$mean_loss,
    car_95_K1 = div_0$results[["1"]]$car_95,
    car_95_K100 = div_0$results[["100"]]$car_95,
    correlation = "0 (Independent)"
  ),
  data.frame(
    year = time_horizons,
    mean_loss = div_05$results[["1"]]$mean_loss,
    car_95_K1 = div_05$results[["1"]]$car_95,
    car_95_K100 = div_05$results[["100"]]$car_95,
    correlation = "0.5 (Moderate)"
  ),
  data.frame(
    year = time_horizons,
    mean_loss = div_1$results[["1"]]$mean_loss,
    car_95_K1 = div_1$results[["1"]]$car_95,
    car_95_K100 = div_1$results[["100"]]$car_95,
    correlation = "1 (Perfect)"
  )
)

# Reshape for plotting
correlation_plot_df <- correlation_comparison_df %>%
  pivot_longer(
    cols = c(mean_loss, car_95_K1, car_95_K100),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = factor(
      series,
      levels = c("mean_loss", "car_95_K1", "car_95_K100"),
      labels = c("Expected Value", "1 Project (95% CaR)", "100 Projects (95% CaR)")
    ),
    correlation = factor(correlation,
                         levels = c("0 (Independent)", "0.5 (Moderate)", "1 (Perfect)"))
  )

# Plot ------------------------------------------------------------------------
p <- ggplot(
  correlation_plot_df,
  aes(x = year, y = value, color = series, linetype = series)
) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Expected Value" = "black",
    "1 Project (95% CaR)" = "firebrick",
    "100 Projects (95% CaR)" = "steelblue"
  )) +
  scale_linetype_manual(values = c(
    "Expected Value" = "solid",
    "1 Project (95% CaR)" = "dashed",
    "100 Projects (95% CaR)" = "dashed"
  )) +
  labs(
    title = "Impact of Correlation on Diversification Benefit",
    subtitle = div_0$region,
    x = "Time (years)",
    y = "Carbon at Risk (ha)",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  facet_wrap(~correlation)

ggsave(file.path(out_dir, "si_correlation_impact.pdf"), p,
       width = 10, height = 6)

cat("\nSaved to outputs/si/si_correlation_impact.pdf\n")
