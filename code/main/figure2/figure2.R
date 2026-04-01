# Figure 2: Fire Risk, Diversification, and Spatial Correlation
# 6-panel figure for Nature Climate Change submission
#
# Prerequisites:
#   - Spatial correlation pipeline (01 -> 02 -> 03) must be run first
#   - EFFIS API access (or cached responses in data/effis_cache/)
#
# Run from repo root:
#   Rscript code/main/figure2/figure2.R
#
# Output:
#   - outputs/main/figure2.pdf (composite figure)
#   - outputs/main/figure2_{a..e}.pdf (individual panels)

library(sf)
library(httr)
library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
select <- dplyr::select

source("code/0_funcs/fire_funcs.R")
source("code/0_funcs/regrowth_funcs.R")

# Configuration ---------------------------------------------------------------
GPKG_PATH <- "data/admin_regrowth_with_gpp.gpkg"
EFFIS_CACHE <- "data/effis_cache"
PROJECT_AREA <- 100000
CLIMATE_RATE <- 0.005
RESCALE_FIRESIZE <- FALSE
N_SIMULATIONS <- 1000
TIME_HORIZONS <- c(1, seq(10, 200, 10))
REGROWTH_RATES <- get_regrowth_rates()

SELECTED_REGIONS <- list(
  list(country = "United States", subcountry = "California"),
  list(country = "Brazil", subcountry = "Mato Grosso"),
  list(country = "Indonesia", subcountry = "Papua")
)

# Color scheme
CA_COLOR <- "#E67E22"
CA_COLOR_LIGHT <- "#F5B041"

# Output directory
out_dir <- "outputs/main"
subfig_dir <- file.path(out_dir, "subfigs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(subfig_dir)) dir.create(subfig_dir, recursive = TRUE)

# Figure width for Nature Climate Change (183mm single-column)
FIG_WIDTH_MM <- 183
FIG_WIDTH_IN <- FIG_WIDTH_MM / 25.4  # ~7.2 inches

set.seed(101)

# Load data -------------------------------------------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)
results <- process_selected_geometries(
  gdf_gadm, SELECTED_REGIONS,
  regrowth_rates = REGROWTH_RATES,
  rescale_firesize = RESCALE_FIRESIZE,
  project_area = PROJECT_AREA,
  n_simulations = N_SIMULATIONS,
  time_horizons = TIME_HORIZONS,
  cache_dir = EFFIS_CACHE
)

# Load correlation analysis results
corr_rds <- "outputs/intermediate/correlation_results/cell_1.00deg/within_region.rds"
within_region <- readRDS(corr_rds)

ca_decay_fit <- within_region$California$decay_fit
ca_rho0 <- coef(ca_decay_fit)["a"]
ca_lambda_km <- coef(ca_decay_fit)["lambda"]

# Panel (a): Burn boxplots with regrowth rate ---------------------------------
years_to_2100 <- 80
climate_factor_2100 <- 1 + CLIMATE_RATE * years_to_2100

burn_comparison_df <- results$empirical_burn_fractions %>%
  mutate(period = "Historic") %>%
  bind_rows(
    results$empirical_burn_fractions %>%
      mutate(
        burn_fraction = burn_fraction * climate_factor_2100,
        period = "Projected 2100"
      )
  ) %>%
  mutate(
    period = factor(period, levels = c("Historic", "Projected 2100")),
    region_short = case_when(
      geo_label == "United States - California" ~ "California",
      geo_label == "Brazil - Mato Grosso" ~ "Mato Grosso",
      geo_label == "Indonesia - Papua" ~ "Papua"
    )
  )

results$regrowth <- results$regrowth %>%
  mutate(region_short = case_when(
    geo_label == "United States - California" ~ "California",
    geo_label == "Brazil - Mato Grosso" ~ "Mato Grosso",
    geo_label == "Indonesia - Papua" ~ "Papua"
  ))

panel_a <- ggplot() +
  geom_boxplot(
    data = burn_comparison_df,
    aes(x = region_short, y = burn_fraction, fill = period, color = period),
    position = position_dodge(width = 0.8),
    width = 0.7, alpha = 0.7
  ) +
  geom_segment(
    data = results$regrowth,
    aes(x = as.numeric(factor(region_short)) - 0.4,
        xend = as.numeric(factor(region_short)) + 0.4,
        y = regrowth_rate, yend = regrowth_rate),
    linewidth = 0.8, color = "darkgreen"
  ) +
  geom_hline(yintercept = 0, color = "gray50") +
  scale_fill_manual(
    values = c("Historic" = "steelblue", "Projected 2100" = "firebrick"),
    name = NULL
  ) +
  scale_color_manual(
    values = c("Historic" = "steelblue", "Projected 2100" = "firebrick"),
    name = NULL
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0, 0.18)) +
  labs(title = NULL, x = NULL, y = "Annual Rate") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.97, 0.88),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    # legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    legend.margin = margin(2, 4, 2, 4),
    axis.text.x = element_text(size = 7)
  ) +
  # Add regrowth line+label inside the legend box using a dummy aesthetic
  annotate("segment", x = 1.93, xend = 2.1,
           y = 0.175, yend = 0.175,
           linewidth = 0.8, color = "darkgreen") +
  annotate("text", x = 2.22, y = 0.175,
           label = "Regrowth", hjust = 0, size = 2.5)

# Panel (b): CaR curves by region ---------------------------------------------
car_results_short <- results$car_results %>%
  mutate(
    car_kg_tonne = car_95 * 1000,
    region_short = case_when(
      geo_label == "United States - California" ~ "California",
      geo_label == "Brazil - Mato Grosso" ~ "Mato Grosso",
      geo_label == "Indonesia - Papua" ~ "Papua"
    )
  )

panel_b <- ggplot(car_results_short, aes(x = time_horizon, y = car_kg_tonne, color = region_short)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.8) +
  scale_y_continuous(limits = c(0, 1050)) +
  scale_color_manual(
    values = c("California" = CA_COLOR, "Mato Grosso" = "#7F8C8D", "Papua" = "#95A5A6")
  ) +
  labs(
    title = NULL,
    x = "Time Horizon (years)",
    y = "CaR (kg/tonne)",
    color = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Panel (c): Correlation vs distance (California) -----------------------------
ca_pairs <- within_region$California$pairs

panel_c <- plot_correlation_vs_distance(
  ca_pairs,
  rho0 = ca_rho0,
  lambda_km = ca_lambda_km,
  use_ribbon = TRUE,
  bin_width = 100,
  ca_color = CA_COLOR
) +
  labs(title = NULL) +
  theme_classic(base_size = 8) +
  theme(legend.text = element_text(size = 7))

# Panel (d): PDF comparison ---------------------------------------------------
rho_values <- c(0, 0.1, 0.25)
pdf_data <- map_dfr(rho_values, function(rho) {
  cat(sprintf("Running simulations for rho = %.2f...\n", rho))
  sim_result <- run_diversification_analysis(
    gdf_gadm, regrowth_rates = REGROWTH_RATES,
    correlation = rho, n_simulations = 1000, return_raw = TRUE,
    rescale_firesize = RESCALE_FIRESIZE,
    cache_dir = EFFIS_CACHE
  )
  bind_rows(
    tibble(
      carbon_loss = sim_result$results[['1']]$mc_hectares[, 101],
      rho = rho,
      projects = "K = 1"
    ),
    tibble(
      carbon_loss = sim_result$results[['100']]$mc_hectares[, 101],
      rho = rho,
      projects = "K = 100"
    )
  )
}) %>%
  mutate(
    rho_label = factor(
      sprintf("rho == %.2f", rho),
      levels = sprintf("rho == %.2f", rho_values)
    ),
    carbon_removed = 1000 - carbon_loss
  )

car_data_e <- pdf_data %>%
  group_by(rho_label, projects, rho) %>%
  summarise(
    p5_removed = quantile(carbon_removed, 0.05),
    car = 1000 - quantile(carbon_removed, 0.05),
    .groups = "drop"
  )

div_benefit_e <- car_data_e %>%
  select(rho_label, rho, projects, car) %>%
  pivot_wider(names_from = projects, values_from = car) %>%
  mutate(div_benefit = `K = 1` - `K = 100`)

panel_d <- ggplot(pdf_data, aes(x = carbon_removed, fill = projects)) +
  geom_density(alpha = 0.6, color = NA) +
  geom_vline(xintercept = 1000, linetype = "solid", linewidth = 0.2, color = "darkgreen") +
  geom_vline(
    data = car_data_e,
    aes(xintercept = p5_removed, color = projects),
    linetype = "solid", linewidth = 0.3
  ) +
  facet_wrap(~rho_label, nrow = 1, labeller = label_parsed) +
  scale_fill_manual(values = c("K = 1" = CA_COLOR, "K = 100" = "#6BAED6")) +
  scale_color_manual(values = c("K = 1" = "#B7410E", "K = 100" = "#2171B5"), guide = "none") +
  labs(
    title = NULL,
    x = "Carbon Removed (kg/tonne)",
    y = "Density",
    fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    strip.text = element_text(size = 8)
  ) +
  geom_segment(
    data = car_data_e %>% filter(rho == 0, projects == "K = 1"),
    aes(x = 1000, xend = p5_removed, y = 0.028, yend = 0.028),
    arrow = arrow(length = unit(0.07, "cm"), ends = "both", type = "closed"),
    color = "#B7410E", linewidth = 0.3, inherit.aes = FALSE
  ) +
  geom_text(
    data = car_data_e %>% filter(rho == 0, projects == "K = 1"),
    aes(x = (1000 + p5_removed) / 2, y = 0.035, label = paste0("CaR (K=1)\n", round(car), " kg/t")),
    color = "#B7410E", size = 2.2, inherit.aes = FALSE
  ) +
  geom_segment(
    data = car_data_e %>% filter(rho == 0, projects == "K = 100"),
    aes(x = 1000, xend = p5_removed, y = 0.015, yend = 0.015),
    arrow = arrow(length = unit(0.07, "cm"), ends = "both", type = "closed"),
    color = "#2171B5", linewidth = 0.3, inherit.aes = FALSE
  ) +
  geom_text(
    data = car_data_e %>% filter(rho == 0, projects == "K = 100"),
    aes(x = (1000 + p5_removed) / 2, y = 0.01, label = paste0("CaR (K=100)\n", round(car), " kg/t")),
    color = "#2171B5", size = 2.2, inherit.aes = FALSE
  ) +
  geom_segment(
    data = div_benefit_e %>% filter(rho == 0),
    aes(x = 1000 - `K = 1`, xend = 1000 - `K = 100`, y = 0.05, yend = 0.05),
    arrow = arrow(length = unit(0.07, "cm"), ends = "both", type = "closed"),
    color = "gray30", linewidth = 0.3, inherit.aes = FALSE
  ) +
  geom_text(
    data = div_benefit_e %>% filter(rho == 0),
    aes(x = 1000 - `K = 100` + 30, y = 0.05,
        label = paste0("Div. benefit: ", round(div_benefit), " kg/t")),
    color = "gray30", size = 2.2, hjust = 0, inherit.aes = FALSE
  )

# Panel (e): Benefit vs rho ---------------------------------------------------
cat("Running correlation sweep...\n")
sweep_df <- map_dfr(seq(0, 1, 0.05), function(rho) {
  run_diversification_analysis(gdf_gadm, regrowth_rates = REGROWTH_RATES,
                               correlation = rho, n_simulations = 1000,
                               rescale_firesize = RESCALE_FIRESIZE,
                               cache_dir = EFFIS_CACHE)$summary %>%
    mutate(rho = !!rho)
}) %>%
  filter(Horizon == 100) %>%
  mutate(reduction_pct = 100 * (CaR_K100 - CaR_K1) / CaR_K1)

panel_e <- ggplot(sweep_df, aes(x = rho, y = reduction_pct)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, color = "gray30", span = 0.4) +
  geom_point(size = 0.8, color = "gray30") +
  geom_point(
    data = filter(sweep_df, rho %in% c(0, 0.1, 0.25)),
    size = 2, color = CA_COLOR
  ) +
  geom_text(
    data = filter(sweep_df, rho %in% c(0, 0.1, 0.25)) %>%
      mutate(lbl = case_when(
        rho == 0    ~ "\"Independent\"",
        rho == 0.1  ~ "rho==0.10",
        rho == 0.25 ~ "rho==0.25"
      )),
    aes(label = lbl),
    hjust = 0, nudge_x = 0.03, size = 2.2, color = CA_COLOR, parse = TRUE
  ) +
  labs(
    title = NULL,
    x = expression("Correlation between projects ("*rho*")"),
    y = "CaR reduction (%)\nfrom K=1 to K=100"
  ) +
  theme_classic(base_size = 8)

# Assemble figure -------------------------------------------------------------
cat("Assembling figure...\n")

fig2 <- (panel_a | panel_b | panel_c) /
        (panel_d + panel_e + plot_layout(widths = c(2, 1))) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 9)
  )

# Save outputs ----------------------------------------------------------------
print(fig2)
ggsave(file.path(out_dir, "figure2.pdf"), fig2,
       width = FIG_WIDTH_MM, height = FIG_WIDTH_MM * 9 / 14,
       units = "mm")

# Individual panels (subfigures at original readable sizes)
ggsave(file.path(subfig_dir, "figure2_a.pdf"), panel_a, width = 7, height = 5)
ggsave(file.path(subfig_dir, "figure2_b.pdf"), panel_b, width = 6, height = 4)
ggsave(file.path(subfig_dir, "figure2_c.pdf"), panel_c, width = 5, height = 4)
ggsave(file.path(subfig_dir, "figure2_d.pdf"), panel_d, width = 10, height = 5)
ggsave(file.path(subfig_dir, "figure2_e.pdf"), panel_e, width = 5, height = 4)

cat("Figure 2 saved to outputs/main/figure2.pdf\n")
cat("Individual panels saved to outputs/main/subfigs/\n")
