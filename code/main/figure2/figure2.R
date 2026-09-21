# Figure 2: Fire Risk and Diversification
#
# Main text (3 panels, single row):
#   a  Annual burn fractions by region, historic and projected 2100
#   b  95% CaR over 1-200 year horizons, three regions
#   c  CaR reduction from K=1 to K=100 against inter-project correlation
#
# SI (also built here):
#   si_rho_distance.pdf  Pairwise correlation against distance (California)
#   si_density_k.pdf     Delivery densities, K=1 vs K=100, three correlations
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
#   - outputs/main/subfigs/figure2_{a..e}.pdf (individual panels)
#   - outputs/si/si_rho_distance.pdf, si_density_k.pdf

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
CLIMATE_RATE <- 0.005
RESCALE_FIRESIZE <- FALSE
N_SIMULATIONS <- 5000
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

# Every simulation call site is seeded individually rather than relying on one
# seed at the top of the script. Calls otherwise draw from different points in
# the stream, so the same quantity (the K=1 California CaR) came out differently
# in each panel.
set.seed(CAR_SEED)

# Load data -------------------------------------------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)
set.seed(CAR_SEED)
results <- process_selected_geometries(
  gdf_gadm, SELECTED_REGIONS,
  regrowth_rates = REGROWTH_RATES,
  rescale_firesize = RESCALE_FIRESIZE,
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
    y = "95% CaR (kg per tonne contracted)",
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
  set.seed(CAR_SEED)
  sim_result <- run_diversification_analysis(
    gdf_gadm, regrowth_rates = REGROWTH_RATES,
    correlation = rho, n_simulations = N_SIMULATIONS, return_raw = TRUE,
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
    x = "Carbon removed (kg per tonne contracted)",
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
  set.seed(CAR_SEED)
  run_diversification_analysis(gdf_gadm, regrowth_rates = REGROWTH_RATES,
                               correlation = rho, n_simulations = N_SIMULATIONS,
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
    y = "Change in CaR, K=1 to K=100\n(% of single-project CaR)"
  ) +
  theme_classic(base_size = 8)

# SI figures: panels moved out of the main text -------------------------------
si_dir <- "outputs/si"
if (!dir.exists(si_dir)) dir.create(si_dir, recursive = TRUE)

# Rendered at the width they are displayed at in supplement.tex, so the panel
# base_size is the true on-page point size.
SI_LINEWIDTH_IN <- 500.484 / 72.27
SI_FRAC_RHO <- 0.40
SI_FRAC_DENS <- 0.60

ggsave(file.path(si_dir, "si_rho_distance.pdf"),
       panel_c + theme_classic(base_size = 8),
       width = SI_FRAC_RHO * SI_LINEWIDTH_IN,
       height = SI_FRAC_RHO * SI_LINEWIDTH_IN * 0.85)
ggsave(file.path(si_dir, "si_density_k.pdf"), panel_d,
       width = SI_FRAC_DENS * SI_LINEWIDTH_IN,
       height = SI_FRAC_DENS * SI_LINEWIDTH_IN * 0.55)

cat("\n--- CaR vs CVaR at 95%, kg per tonne contracted (SI table) ---\n")
results$car_results %>%
  filter(time_horizon %in% c(100, 200)) %>%
  mutate(
    region = case_when(
      geo_label == "United States - California" ~ "California",
      geo_label == "Brazil - Mato Grosso"       ~ "Mato Grosso",
      geo_label == "Indonesia - Papua"          ~ "Papua"
    ),
    CaR_95  = round(car_95 * 1000),
    CVaR_95 = round(cvar_95 * 1000),
    uplift  = sprintf("%.1f%%", 100 * (cvar_95 / car_95 - 1))
  ) %>%
  select(region, time_horizon, CaR_95, CVaR_95, uplift) %>%
  arrange(region, time_horizon) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat("\n")

cat(sprintf("Correlation decay fit: rho0 = %.4f, lambda = %.1f km\n",
            ca_rho0, ca_lambda_km))
cat(sprintf("California median pairwise rho, all pairs: %.4f\n",
            median(ca_pairs$cor_spearman)))
cat(sprintf("California median pairwise rho, <100 km:  %.4f\n",
            median(ca_pairs$cor_spearman[ca_pairs$distance_km < 100])))
cat(sprintf("California median pairwise rho, >500 km:  %.4f\n",
            median(ca_pairs$cor_spearman[ca_pairs$distance_km > 500])))
cat(sprintf("Mato Grosso median pairwise rho, all pairs: %.4f\n",
            median(within_region$Mato_Grosso$pairs$cor_spearman)))

# Machine-readable export -----------------------------------------------------
# Every forest number quoted in the main text or SI, in one long-format file, so
# transcription into the manuscript can be checked rather than read off the
# figure annotations.

num_dir <- "outputs/intermediate"
if (!dir.exists(num_dir)) dir.create(num_dir, recursive = TRUE)

num_single <- results$car_results %>%
  filter(time_horizon %in% c(30, 100, 200)) %>%
  mutate(region = case_when(
    geo_label == "United States - California" ~ "California",
    geo_label == "Brazil - Mato Grosso"       ~ "Mato Grosso",
    geo_label == "Indonesia - Papua"          ~ "Papua"
  )) %>%
  select(region, horizon_yr = time_horizon, mean_loss = mean_car,
         car_95, cvar_95) %>%
  pivot_longer(c(mean_loss, car_95, cvar_95),
               names_to = "statistic", values_to = "value") %>%
  mutate(quantity = "single_project", rho = NA_real_,
         value = round(value * 1000, 1), unit = "kg_per_tonne")

num_cvar_uplift <- results$car_results %>%
  filter(time_horizon %in% c(30, 100, 200)) %>%
  mutate(region = case_when(
    geo_label == "United States - California" ~ "California",
    geo_label == "Brazil - Mato Grosso"       ~ "Mato Grosso",
    geo_label == "Indonesia - Papua"          ~ "Papua"
  )) %>%
  transmute(quantity = "single_project", region, horizon_yr = time_horizon,
            rho = NA_real_, statistic = "cvar_uplift",
            value = round(100 * (cvar_95 / car_95 - 1), 2), unit = "percent")

# Panel d: the three correlation levels the main text and SI quote.
num_panel_d <- div_benefit_e %>%
  transmute(rho,
            car_K1 = `K = 1`, car_K100 = `K = 100`,
            div_benefit_kg = div_benefit,
            div_benefit_pct = 100 * div_benefit / `K = 1`) %>%
  pivot_longer(-rho, names_to = "statistic", values_to = "value") %>%
  mutate(quantity = "portfolio_panel_d", region = "California",
         horizon_yr = 100,
         unit = if_else(statistic == "div_benefit_pct",
                        "percent", "kg_per_tonne"),
         value = round(value, 2))

# Panel e: the full rho grid behind the sweep.
num_panel_e <- sweep_df %>%
  transmute(rho, car_K1 = CaR_K1, car_K10 = CaR_K10, car_K100 = CaR_K100,
            div_benefit_kg = CaR_K1 - CaR_K100,
            div_benefit_pct = -reduction_pct) %>%
  pivot_longer(-rho, names_to = "statistic", values_to = "value") %>%
  mutate(quantity = "portfolio_panel_e", region = "California",
         horizon_yr = 100,
         unit = if_else(statistic == "div_benefit_pct",
                        "percent", "kg_per_tonne"),
         value = round(value, 2))

figure2_numbers <- bind_rows(num_single, num_cvar_uplift,
                             num_panel_d, num_panel_e) %>%
  select(quantity, region, horizon_yr, rho, statistic, value, unit)

write.csv(figure2_numbers,
          file.path(num_dir, "figure2_numbers.csv"), row.names = FALSE)
cat("\nWrote", file.path(num_dir, "figure2_numbers.csv"),
    sprintf("(%d rows)\n", nrow(figure2_numbers)))

# The K=1 California CaR at 100 years is quoted in three places and must agree
# across all of them now that each call site is seeded identically.
k1_panel_b <- figure2_numbers %>%
  filter(quantity == "single_project", region == "California",
         horizon_yr == 100, statistic == "car_95") %>% pull(value)
k1_panel_d <- figure2_numbers %>%
  filter(quantity == "portfolio_panel_d", rho == 0,
         statistic == "car_K1") %>% pull(value)
k1_panel_e <- figure2_numbers %>%
  filter(quantity == "portfolio_panel_e", rho == 0,
         statistic == "car_K1") %>% pull(value)
cat(sprintf("K=1 California CaR95 at 100 yr: panel b %.1f | panel d %.1f | panel e %.1f\n",
            k1_panel_b, k1_panel_d, k1_panel_e))
if (max(abs(c(k1_panel_d, k1_panel_e) - k1_panel_b)) > 0.5) {
  warning("K=1 California CaR disagrees across panels; check the per-call-site seeding.")
}

# Assemble figure -------------------------------------------------------------
cat("Assembling figure...\n")

# `patchwork & theme` errors under ggplot2 4.x, so add the tag theme per panel.
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 9))

fig2 <- ((panel_a + tag_theme) | (panel_b + tag_theme) | (panel_e + tag_theme)) +
  plot_annotation(tag_levels = "a")

# Save outputs ----------------------------------------------------------------
print(fig2)
ggsave(file.path(out_dir, "figure2.pdf"), fig2,
       width = FIG_WIDTH_MM, height = FIG_WIDTH_MM * 5 / 14,
       units = "mm")

# Individual panels (subfigures at original readable sizes)
ggsave(file.path(subfig_dir, "figure2_a.pdf"), panel_a, width = 7, height = 5)
ggsave(file.path(subfig_dir, "figure2_b.pdf"), panel_b, width = 6, height = 4)
ggsave(file.path(subfig_dir, "figure2_c.pdf"), panel_c, width = 5, height = 4)
ggsave(file.path(subfig_dir, "figure2_d.pdf"), panel_d, width = 10, height = 5)
ggsave(file.path(subfig_dir, "figure2_e.pdf"), panel_e, width = 5, height = 4)

cat("Figure 2 saved to outputs/main/figure2.pdf\n")
cat("Individual panels saved to outputs/main/subfigs/\n")
