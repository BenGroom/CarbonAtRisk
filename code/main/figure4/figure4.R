# Figure 4: Three procurement rules — how correlation drives diversification
#
# Three procurement rules:
#   Rule 1 "Cheapest tonne":  min cost s.t. Q >= T  (all forest)
#   Rule 2 "Only permanent":  all-DACCS s.t. p5 >= T (all DACCS)
#   Rule 3 "CaR-constrained": min cost s.t. p5 >= T  (diversified mix)
#
# Layout (2 rows x 3 panels):
#   Row 1: (a) Portfolio composition  |  (b) Delivery distributions  |  (c) Mean-SD space
#   Row 2: (d) Budget dual bars       |  (e) Forest share vs rho     |  (f) c_eff vs rho
#
# Also generates SI figures:
#   - si_portfolio_cost_vs_rho.pdf
#   - si_portfolio_eff_price_heatmap.pdf
#
#
# Run from repo root:  Rscript code/main/figure4/figure4.R

library(tidyverse)
library(truncnorm)
library(patchwork)
library(scales)
library(ggrepel)

FIG_WIDTH_MM <- 183
FIG_WIDTH_IN <- FIG_WIDTH_MM / 25.4

theme_set(
  theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 7),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.3, "cm")
    )
)

set.seed(42)
source("code/0_funcs/portfolio_funcs.R")

# "Parameters" ----------

project_costs  <- c(DACS = 500, Forest = 40)
survival_probs <- c(DACS = 0.99, Forest = 0.4)
q         <- 500
z_p5      <- qnorm(0.05)
T_target  <- 30000

scenarios <- list(
  "Independent" = c(DACS = 0, Forest = 0),
  "Correlated"  = c(DACS = 0.5, Forest = 0.2)
)

grid <- expand_grid(
  n_DACS   = seq(0, 500, by = 1),
  n_Forest = seq(0, 500, by = 1)
)

rule_cols <- c(
  "R1: Cheapest"  = "#e41a1c",
  "R2: Permanent" = "#ff7f00",
  "R3: CaR"       = "#377eb8"
)
tech_cols <- c("DACCS" = "#984ea3", "Forest" = "#4daf4a")

out_dir <- "outputs/main"
si_dir  <- "outputs/si"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(si_dir, showWarnings = FALSE, recursive = TRUE)

# "Solve cost-min rules" ----------

solve_costmin <- function(rho_within, scenario_name) {
  stats <- portfolio_stats(
    df = grid, q = q, p = survival_probs, costs = project_costs,
    z = z_p5, rho_within = rho_within, rho_between = 0
  ) %>% mutate(Q = n_total * q)

  r1 <- stats %>% filter(Q >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R1: Cheapest")
  r2 <- stats %>% filter(n_Forest == 0, Q >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R2: Permanent")
  r3 <- stats %>% filter(p5 >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R3: CaR")

  bind_rows(r1, r2, r3) %>% mutate(scenario = scenario_name)
}

results_cm <- map2_dfr(scenarios, names(scenarios), solve_costmin) %>%
  mutate(
    scenario = factor(scenario, levels = c("Independent", "Correlated")),
    rule = factor(rule, levels = names(rule_cols))
  )

# "Derive budget from R3 correlated cost" ----------

B <- results_cm %>%
  filter(scenario == "Correlated", rule == "R3: CaR") %>%
  pull(cost)
cat("Budget B set to R3 correlated cost:", B, "\n")

# "Solve budget rules (correlated only)" ----------

stats_corr <- portfolio_stats(
  df = grid, q = q, p = survival_probs, costs = project_costs,
  z = z_p5, rho_within = scenarios[["Correlated"]], rho_between = 0
) %>% mutate(Q = n_total * q)

affordable <- stats_corr %>% filter(cost <= B, n_total > 0)
budget_r1 <- affordable %>% slice_max(Q, with_ties = FALSE) %>% mutate(rule = "R1: Cheapest")
budget_r2 <- affordable %>% filter(n_Forest == 0) %>%
  slice_max(p5, with_ties = FALSE) %>% mutate(rule = "R2: Permanent")
budget_r3 <- affordable %>% slice_max(p5, with_ties = FALSE) %>% mutate(rule = "R3: CaR")
results_budget <- bind_rows(budget_r1, budget_r2, budget_r3) %>%
  mutate(rule = factor(rule, levels = names(rule_cols)))

cat("Solutions computed\n")

# "Panel (a): Portfolio composition" ----------

comp_df <- results_cm %>%
  select(scenario, rule, n_DACS, n_Forest) %>%
  pivot_longer(c(n_DACS, n_Forest), names_to = "tech", values_to = "n") %>%
  mutate(tech = recode(tech, n_DACS = "DACCS", n_Forest = "Forest"),
         tech = factor(tech, levels = c("DACCS", "Forest")))

pa <- comp_df %>%
  ggplot(aes(x = rule, y = n, fill = tech)) +
  geom_col(position = "stack", width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_stack(vjust = 0.5),
            size = 2, fontface = "bold", color = "white") +
  facet_wrap(~scenario) +
  scale_fill_manual(values = tech_cols, name = NULL) +
  scale_x_discrete(labels = c("R1", "R2", "R3")) +
  labs(x = NULL, y = "Projects") +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 7))

cat("Panel a done\n")

# "Panel (b): Delivery distributions (correlated)" ----------

pdf_cm <- results_cm %>%
  filter(scenario == "Correlated") %>%
  mutate(lower = 0, upper = n_total * q) %>%
  crossing(x = seq(0, max(results_cm$n_total * q), length.out = 2000)) %>%
  mutate(density = dtruncnorm(x, a = lower, b = upper, mean = mu, sd = sd))

x_upper <- max(results_cm %>% filter(scenario == "Correlated") %>% pull(mu)) * 1.3

pb <- pdf_cm %>%
  filter(x > 0) %>%
  ggplot(aes(x = x, y = density, color = rule, fill = rule)) +
  geom_area(alpha = 0.1, position = "identity") +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = rule_cols, name = NULL) +
  scale_fill_manual(values = rule_cols, name = NULL) +
  scale_x_continuous(labels = comma, limits = c(0, x_upper)) +
  labs(x = expression(tCO[2]~delivered), y = "Density") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.3, linewidth = 0.8)),
         color = "none") +
  theme(legend.position = "bottom")

cat("Panel b done\n")

# "Panel (c): Mean-SD space (correlated)" ----------

# Tighter grid: cover solution region with padding
max_n_dacs   <- max(results_cm$n_DACS) * 1.5
max_n_forest <- max(results_cm$n_Forest) * 2
sparse_grid <- expand_grid(
  n_DACS  = seq(0, ceiling(max_n_dacs), by = 1),
  n_Forest = seq(0, ceiling(max_n_forest), by = 1)
)
frontier_corr <- portfolio_stats(
  df = sparse_grid, q = q, p = survival_probs, costs = project_costs,
  z = z_p5, rho_within = scenarios[["Correlated"]], rho_between = 0
) %>% mutate(meets_p5 = p5 >= T_target)

p5_line <- tibble(
  sd = seq(0, max(frontier_corr$sd, na.rm = TRUE), length.out = 200),
  mu = T_target - z_p5 * sd
)

hl <- results_cm %>%
  filter(scenario == "Correlated") %>%
  select(rule, sd, mu, cost, n_DACS, n_Forest) %>%
  mutate(label = paste0("R", match(rule, names(rule_cols)), ": ", n_DACS, "D+", n_Forest, "F"))

pc <- ggplot() +
  geom_point(data = frontier_corr %>% filter(!meets_p5, n_total > 0),
             aes(x = sd, y = mu), color = "grey90", size = 0.1, shape = 16) +
  geom_point(data = frontier_corr %>% filter(meets_p5),
             aes(x = sd, y = mu, color = cost / 1000), size = 0.2, shape = 16) +
  geom_line(data = p5_line, aes(x = sd, y = mu),
            color = "#e41a1c", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = hl, aes(x = sd, y = mu),
             shape = 21, size = 2.5, fill = "white", color = "black", stroke = 0.6) +
  geom_text_repel(data = hl, aes(x = sd, y = mu, label = label),
                  size = 2, fontface = "bold", lineheight = 0.8,
                  min.segment.length = 0, nudge_y = 4000,
                  point.padding = unit(0.3, "lines")) +
  scale_color_viridis_c(name = "Cost ($k)", option = "D", direction = -1,
                        labels = function(x) paste0("$", x, "k")) +
  scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(10000, 60000)) +
  labs(x = expression(SD~(tCO[2])), y = expression(mu~(tCO[2]))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm")) +
  guides(color = guide_colorbar(barwidth = 5, barheight = 0.3))

cat("Panel c done\n")

# "Panel (d): Budget dual bars" ----------

outcome_df <- results_budget %>%
  mutate(Q = n_total * q) %>%
  select(rule, Q, mu, p5) %>%
  pivot_longer(c(Q, mu, p5), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Q", "mu", "p5"),
                         labels = c("Q (contracted)", "mu (expected)", "p5 (95% guaranteed)")))

pd <- outcome_df %>%
  ggplot(aes(x = rule, y = value / 1000, fill = metric)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = T_target / 1000, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Q (contracted)" = "#bdbdbd",
                               "mu (expected)" = "#6baed6",
                               "p5 (95% guaranteed)" = "#d6604d"), name = NULL) +
  scale_x_discrete(labels = c("R1\ncheapest", "R2\npermanent", "R3\nCaR")) +
  labs(x = NULL, y = expression(tCO[2]~"(thousands)")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 6))

cat("Panel d done\n")

# "Panel (e): Forest share vs rho" ----------

rho_grid <- seq(0, 0.4, by = 0.005)

sweep <- map_dfr(rho_grid, function(rho_F) {
  st <- portfolio_stats(
    df = grid, q = q, p = survival_probs, costs = project_costs,
    z = z_p5, rho_within = c(DACS = 0.5, Forest = rho_F), rho_between = 0
  ) %>% mutate(Q = n_total * q)

  r1 <- st %>% filter(Q >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R1: Cheapest")
  r2 <- st %>% filter(n_Forest == 0, Q >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R2: Permanent")
  r3 <- st %>% filter(p5 >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R3: CaR")
  bind_rows(r1, r2, r3) %>% mutate(rho_F = rho_F)
})

pe <- sweep %>%
  mutate(rule = factor(rule, levels = names(rule_cols))) %>%
  ggplot(aes(x = rho_F, y = forest_share, color = rule)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = rule_cols, name = NULL) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
  labs(x = expression(rho[Forest]), y = "Forest share") +
  theme(legend.position = "bottom")

cat("Panel e done\n")

# "Panel (f): Effective cost per guaranteed tonne vs rho" ----------

pf_data <- sweep %>%
  mutate(
    rule  = factor(rule, levels = names(rule_cols)),
    c_eff = if_else(p5 > 0, cost / p5, NA_real_)
  )

# Y-axis cap: 2x max of R2/R3 so R1's spike doesn't squash everything
c_eff_cap <- pf_data %>%
  filter(rule != "R1: Cheapest") %>%
  pull(c_eff) %>% max(na.rm = TRUE) * 2

# R1 solid where defined and below cap, dashed where above cap
pf_r1 <- pf_data %>% filter(rule == "R1: Cheapest", !is.na(c_eff))
pf_rest <- pf_data %>% filter(rule != "R1: Cheapest")

pf <- ggplot() +
  geom_line(data = pf_rest, aes(x = rho_F, y = c_eff, color = rule),
            linewidth = 0.5) +
  geom_line(data = pf_r1 %>% filter(c_eff <= c_eff_cap),
            aes(x = rho_F, y = c_eff, color = rule), linewidth = 0.5) +
  geom_line(data = pf_r1 %>% filter(c_eff > c_eff_cap),
            aes(x = rho_F, y = c_eff, color = rule),
            linewidth = 0.5, linetype = "dashed") +
  scale_color_manual(values = rule_cols, name = NULL) +
  scale_x_continuous(breaks = seq(0, 0.4, 0.1)) +
  coord_cartesian(ylim = c(0, c_eff_cap)) +
  labs(x = expression(rho[Forest]),
       y = expression(c[eff*","*95]~"($/tCO"[2]*")")) +
  theme(legend.position = "none")

cat("Panel f done\n")

# "Assemble main figure" ----------

fig <- (pa + pb + pc) / (pd + pe + pf) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 8))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.margin = margin(2, 6, 2, 2))

ggsave(file.path(out_dir, "figure4.pdf"), fig,
       width = FIG_WIDTH_MM * 1.05, height = FIG_WIDTH_MM * 0.9,
       units = "mm", bg = "white")
cat("Figure 4 saved to", file.path(out_dir, "figure4.pdf"), "\n")

# "SI figures" ----------
cat("Generating SI portfolio figures...\n")

# SI: Min cost vs rho_forest for three rho_DACCS values
rho_dacs_grid <- c(0, 0.2, 0.5)
rho_forest_grid <- seq(0, 0.8, by = 0.01)

sweep_grid <- expand_grid(n_DACS = seq(0, 500, by = 1), n_Forest = seq(0, 500, by = 1))

sweep_results <- expand_grid(
  rho_within_forest = rho_forest_grid,
  rho_within_dacs   = rho_dacs_grid
) %>%
  pmap_dfr(function(rho_within_forest, rho_within_dacs) {
    portfolio_stats(
      df = sweep_grid, q = q, p = survival_probs, costs = project_costs,
      z = z_p5,
      rho_within = c(DACS = rho_within_dacs, Forest = rho_within_forest),
      rho_between = 0
    ) %>%
      pick_min_cost_feasible(target = T_target, target_var = "p5") %>%
      mutate(rho_within_forest = rho_within_forest,
             rho_within_dacs   = rho_within_dacs)
  })

rho_colors <- c("0" = "steelblue", "0.2" = "darkorange", "0.5" = "firebrick")

si_cost_rho <- sweep_results %>%
  ggplot(aes(x = rho_within_forest, y = cost,
             color = factor(rho_within_dacs),
             group = factor(rho_within_dacs))) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = rho_colors, name = expression(rho[DACS])) +
  scale_y_continuous(labels = scales::dollar) +
  labs(x = expression(rho[Forest]), y = "Minimum cost") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(si_dir, "si_portfolio_cost_vs_rho.pdf"),
       si_cost_rho, width = 7, height = 5)
cat("  si_portfolio_cost_vs_rho.pdf saved\n")

# SI: Effective price heatmap
heat_max <- max(results_cm$n_DACS, results_cm$n_Forest) * 2.5
heat_grid <- expand_grid(n_DACS = seq(0, heat_max, by = 2), n_Forest = seq(0, heat_max, by = 2))

heat_df_indep <- portfolio_stats(
  df = heat_grid, q = q, p = survival_probs, costs = project_costs,
  z = z_p5, rho_within = scenarios[["Independent"]], rho_between = 0
) %>% filter(n_total > 0) %>%
  mutate(c_eff_port = cost / p5, feasible = p5 >= T_target)

heat_df_cal <- portfolio_stats(
  df = heat_grid, q = q, p = survival_probs, costs = project_costs,
  z = z_p5, rho_within = scenarios[["Correlated"]], rho_between = 0
) %>% filter(n_total > 0) %>%
  mutate(c_eff_port = cost / p5, feasible = p5 >= T_target)

opt_indep <- heat_df_indep %>% filter(feasible) %>%
  slice_min(order_by = cost, n = 1, with_ties = FALSE)
opt_cal <- heat_df_cal %>% filter(feasible) %>%
  slice_min(order_by = cost, n = 1, with_ties = FALSE)

make_eff_heatmap <- function(heat_data, opt_pt, subtitle) {
  ggplot() +
    geom_tile(data = heat_data %>% filter(!feasible),
              aes(x = n_DACS, y = n_Forest), fill = "grey90") +
    geom_tile(data = heat_data %>% filter(feasible),
              aes(x = n_DACS, y = n_Forest, fill = c_eff_port)) +
    geom_contour(data = heat_data,
                 aes(x = n_DACS, y = n_Forest, z = p5),
                 breaks = T_target, color = "red", linewidth = 0.6, linetype = "dashed") +
    geom_point(data = opt_pt, aes(x = n_DACS, y = n_Forest),
               shape = 4, size = 4, stroke = 2, color = "red") +
    scale_fill_viridis_c(name = expression(c[eff]~("$"/tCO[2])),
                         option = "D", direction = -1,
                         labels = function(x) paste0("$", sprintf("%.3f", x))) +
    labs(title = subtitle, x = "Number of DACCS projects",
         y = "Number of Forest projects") +
    coord_equal() +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "bottom") +
    guides(fill = guide_colorbar(barwidth = 12, barheight = 0.5))
}

si_eff_heatmap <- make_eff_heatmap(heat_df_indep, opt_indep, "Independent") |
  make_eff_heatmap(heat_df_cal, opt_cal, "Calibrated")
ggsave(file.path(si_dir, "si_portfolio_eff_price_heatmap.pdf"),
       si_eff_heatmap, width = 14, height = 6)
cat("  si_portfolio_eff_price_heatmap.pdf saved\n")

# "Print key statistics" ----------
cat("\nKey results (correlated):\n")
results_cm %>%
  filter(scenario == "Correlated") %>%
  select(rule, n_DACS, n_Forest, cost, mu, p5) %>%
  print()

cat("\nDone.\n")
