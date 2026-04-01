# Figure 3: Carbon at Risk for DACCS
# Two panels per horizon: (a) well-regulated offshore, (b) poorly-regulated onshore
# Uses the Storage Security Calculator from Alcalde et al. (2018)
#
# Run from repo root:
#   Rscript code/main/figure3/figure3.R
#
# Output:
#   - outputs/main/figure3.pdf          (200-year horizon, main text)
#   - outputs/si/si_daccs_1000yr.pdf    (1000-year horizon, SI)
#   - outputs/si/si_daccs_10000yr.pdf   (10000-year horizon, SI)
#   - outputs/intermediate/daccs_mc_results.csv  (all CaR statistics)
#
# Caching: if daccs_mc_results.csv exists, skips MC simulations unless
#   OVERWRITE_DACCS is set to TRUE.

library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(2100)

# Configuration ---------------------------------------------------------------
N_SIMULATIONS <- 10000
OVERWRITE_DACCS <- exists("OVERWRITE_DACCS_FLAG") && OVERWRITE_DACCS_FLAG

main_dir <- "outputs/main"
si_dir   <- "outputs/si"
cache_dir <- "outputs/intermediate"
for (d in c(main_dir, si_dir, cache_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cache_file <- file.path(cache_dir, "daccs_mc_results.csv")

# SSCMC output list indices: 1=1yr, 2=3yr, 3=10yr, 4=30yr, 5=100yr,
#   6=200yr, 7=500yr, 8=1000yr, ..., 17=10000yr
HORIZONS <- data.frame(
  year = c(200, 1000, 10000),
  index = c(6, 8, 17),
  stringsAsFactors = FALSE
)

# Colors
GREEN_MEDIUM <- "#33a02c"
RED_MEDIUM   <- "#fb6a4a"

# Run or load MC simulations ---------------------------------------------------
if (!file.exists(cache_file) || OVERWRITE_DACCS) {
  cat("Running SSC Monte Carlo simulations (N =", N_SIMULATIONS, ")...\n")

  cat("  Offshore scenario...\n")
  source("code/main/figure3/ssc_offshore.R")
  temp_offshore <- SSCMC(N_SIMULATIONS)

  cat("  Onshore scenario...\n")
  source("code/main/figure3/ssc_onshore.R")
  temp_onshore <- SSCMC(N_SIMULATIONS)

  # Extract % loss at each horizon and compute CaR statistics
  results <- list()
  for (i in seq_len(nrow(HORIZONS))) {
    yr <- HORIZONS$year[i]
    idx <- HORIZONS$index[i]

    loss_off <- temp_offshore[[idx]][, 4] / temp_offshore[[idx]][, 2] * 100
    loss_on  <- temp_onshore[[idx]][, 4] / temp_onshore[[idx]][, 2] * 100

    for (scenario in c("offshore", "onshore")) {
      loss <- if (scenario == "offshore") loss_off else loss_on
      stored <- 100 - loss
      results[[length(results) + 1]] <- data.frame(
        scenario = scenario,
        horizon_yr = yr,
        car_90 = round(100 - as.numeric(quantile(stored, 0.10)), 4),
        car_95 = round(100 - as.numeric(quantile(stored, 0.05)), 4),
        car_98 = round(100 - as.numeric(quantile(stored, 0.02)), 4),
        mean_loss = round(mean(loss), 4),
        n_simulations = length(loss)
      )
    }
  }

  car_df <- bind_rows(results)
  write.csv(car_df, cache_file, row.names = FALSE)
  cat("Saved CaR statistics to", cache_file, "\n")

  # Also save the raw loss vectors for plotting
  raw_cache <- file.path(cache_dir, "daccs_mc_raw.rds")
  raw_data <- list()
  for (i in seq_len(nrow(HORIZONS))) {
    yr <- HORIZONS$year[i]
    idx <- HORIZONS$index[i]
    raw_data[[paste0("offshore_", yr)]] <-
      temp_offshore[[idx]][, 4] / temp_offshore[[idx]][, 2] * 100
    raw_data[[paste0("onshore_", yr)]] <-
      temp_onshore[[idx]][, 4] / temp_onshore[[idx]][, 2] * 100
  }
  saveRDS(raw_data, raw_cache)

} else {
  cat("Loading cached DACCS results from", cache_file, "\n")
  car_df <- read.csv(cache_file)
  raw_data <- readRDS(file.path(cache_dir, "daccs_mc_raw.rds"))
}

# Print summary ----------------------------------------------------------------
cat("\nDACCS CaR Summary:\n")
print(car_df)

# Plotting function -----------------------------------------------------------
plot_daccs_car <- function(loss_pct, subtitle) {
  stored_pct <- 100 - loss_pct

  threshold_95 <- as.numeric(quantile(stored_pct, 0.05))
  car_95 <- 100 - threshold_95

  n_sim <- length(stored_pct)
  hist_data <- hist(stored_pct, breaks = 30, plot = FALSE)
  bin_width <- hist_data$breaks[2] - hist_data$breaks[1]

  dens <- density(stored_pct)
  density_df <- data.frame(x = dens$x, y = dens$y * bin_width * n_sim) %>%
    filter(x <= 100)

  label_95 <- paste0("95% CaR = ",
                     ifelse(car_95 >= 1, paste0(round(car_95), "%"),
                            paste0(round(car_95, 2), "%")))
  shade_95 <- density_df %>% filter(x <= threshold_95) %>%
    mutate(band = label_95)
  shade_95$band <- factor(shade_95$band)

  half_y <- max(density_df$y) * 0.5
  arrow_y <- max(density_df$y) * 0.22

  ggplot(data.frame(stored = stored_pct), aes(x = stored)) +
    geom_histogram(aes(y = after_stat(count)), bins = 43,
                   fill = "gray85", color = "white", alpha = 0.6) +
    geom_area(data = shade_95, aes(x = x, y = y, fill = band),
              alpha = 0.6, position = "identity") +
    geom_line(data = density_df, aes(x = x, y = y), color = "black", linewidth = 0.7) +
    annotate("segment", x = threshold_95, xend = threshold_95, y = 0, yend = half_y,
             color = GREEN_MEDIUM, linetype = "dashed", linewidth = 1) +
    annotate("text", x = threshold_95, y = half_y * 1.05, label = "95%",
             color = GREEN_MEDIUM, hjust = 1.1, size = 5.5) +
    annotate("segment", x = threshold_95, xend = 100, y = arrow_y, yend = arrow_y,
             color = RED_MEDIUM, linewidth = 0.8,
             arrow = arrow(length = unit(0.15, "cm"), ends = "both", type = "closed")) +
    geom_vline(xintercept = 100, color = "gray50", linetype = "dotted", linewidth = 0.5) +
    scale_fill_manual(name = "Carbon at Risk (%)",
                      values = setNames(GREEN_MEDIUM, label_95)) +
    labs(title = subtitle,
         x = expression("% CO"[2]*" Stored"),
         y = "Frequency (Number of Simulations)") +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.key.size = unit(0.6, "cm"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 17),
      legend.background = element_rect(fill = "white", color = NA)
    )
}

# Helper to build and save a two-panel figure for a given horizon
save_daccs_figure <- function(horizon_yr, output_path) {
  loss_off <- raw_data[[paste0("offshore_", horizon_yr)]]
  loss_on  <- raw_data[[paste0("onshore_", horizon_yr)]]

  pa <- plot_daccs_car(loss_off, "Offshore")
  pb <- plot_daccs_car(loss_on, "Onshore (unregulated)")

  fig <- (pa | pb) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 18))

  ggsave(output_path, fig, width = 14, height = 6)
  cat("Saved", output_path, "\n")
}

# Generate figures -------------------------------------------------------------

# Main text: 200-year
save_daccs_figure(200, file.path(main_dir, "figure3.pdf"))

# SI: 1000-year and 10000-year
save_daccs_figure(1000, file.path(si_dir, "si_daccs_1000yr.pdf"))
save_daccs_figure(10000, file.path(si_dir, "si_daccs_10000yr.pdf"))

cat("\nDone.\n")
