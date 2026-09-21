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
library(ggrastr)   # panel (a)'s ~150k-point grid is rasterised; axes stay vector

FIG_WIDTH_MM <- 183
FIG_WIDTH_IN <- FIG_WIDTH_MM / 25.4

# theme_set changes ggplot2's default theme for the whole R session, not just
# this script. run_all.R sources every script into one session, so without the
# restore at the end of this file the SI figures that let patchwork supply the
# panel theme (si_deforestation, si_distribution_assumption) render differently
# under run_all.R than they do standalone, and stop matching the manuscript.
default_theme <- theme_get()

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

q         <- 500
z_p5      <- qnorm(0.05)
T_target  <- 30000

# Sticker prices are per tonne CO2; each project delivers q tonnes.
# Only the 12.5:1 ratio enters the optimisation, so rescaling leaves every
# optimal portfolio unchanged and multiplies all reported costs by q.
tonne_prices   <- c(DACS = 500, Forest = 40)
project_costs  <- tonne_prices * q

survival_probs <- c(DACS = 0.99, Forest = 0.4)

scenarios <- list(
  "Independent" = c(DACS = 0, Forest = 0),
  "Correlated"  = c(DACS = 0.5, Forest = 0.2)
)

# Search bounds. n_Forest must reach well past 800: a volume-maximising buyer
# under the budget dual takes ~807 forest projects, and Rule 3's all-forest
# corner reaches ~764 just below the diversification threshold. A bound that
# binds on either manufactures a spurious DACCS holding for Rule 1 and shifts
# the apparent rho_Forest threshold. n_DACS never exceeds ~70 at any solution;
# the `hit` check below reports if either bound is reached.
N_DACS_MAX   <- 200
N_FOREST_MAX <- 2000

grid <- expand_grid(
  n_DACS   = seq(0, N_DACS_MAX, by = 1),
  n_Forest = seq(0, N_FOREST_MAX, by = 1)
)

#' Abort if an optimum sits on the edge of the search grid.
#'
#' Silent boundary solutions are how the rho_Forest threshold and the budget-dual
#' Rule 1 composition were both misreported, so every solve is checked.
check_interior <- function(df, what) {
  hit <- df %>% filter(n_DACS >= N_DACS_MAX | n_Forest >= N_FOREST_MAX)
  if (nrow(hit) > 0) {
    stop("Grid boundary binding in ", what, ": ",
         paste0(hit$n_DACS, "D+", hit$n_Forest, "F", collapse = ", "),
         ". Raise N_DACS_MAX / N_FOREST_MAX.")
  }
  invisible(df)
}

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

  bind_rows(r1, r2, r3) %>% mutate(scenario = scenario_name) %>%
    check_interior(paste("cost-min,", scenario_name))
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
  mutate(rule = factor(rule, levels = names(rule_cols))) %>%
  check_interior("budget dual")

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
             aes(x = sd, y = mu), color = "grey88", size = 0.35, shape = 16) +
  geom_point(data = frontier_corr %>% filter(meets_p5),
             aes(x = sd, y = mu, color = cost / 1e6), size = 0.55, shape = 16) +
  geom_line(data = p5_line, aes(x = sd, y = mu),
            color = "#e41a1c", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = hl, aes(x = sd, y = mu),
             shape = 21, size = 2.5, fill = "white", color = "black", stroke = 0.6) +
  geom_text_repel(data = hl, aes(x = sd, y = mu, label = label),
                  size = 2, fontface = "bold", lineheight = 0.8,
                  min.segment.length = 0, nudge_y = 4000,
                  point.padding = unit(0.3, "lines")) +
  scale_color_viridis_c(name = "Cost ($m)", option = "D", direction = -1,
                        labels = function(x) paste0("$", x, "m")) +
  scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(10000, 60000)) +
  labs(x = expression("Standard deviation of delivery"~(tCO[2])),
       y = expression("Expected delivery"~mu~(tCO[2]))) +
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

# Composition annotation: under a fixed budget each rule maximises its own
# objective, so the labels name the objective rather than the cost-min framing.
comp_lab <- results_budget %>%
  mutate(Q = n_total * q,
         lab = paste0(n_DACS, "D + ", n_Forest, "F"),
         ytop = pmax(Q, mu, p5) / 1000) %>%
  select(rule, lab, ytop)

pd <- outcome_df %>%
  ggplot(aes(x = rule, y = value / 1000, fill = metric)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = T_target / 1000, linetype = "dashed", color = "grey40") +
  geom_text(data = comp_lab, inherit.aes = FALSE,
            aes(x = rule, y = ytop, label = lab),
            vjust = -0.5, size = 1.9, fontface = "bold", color = "grey25") +
  # Plotmath throughout, so D_5 carries its subscript. atop() supplies the
  # two-line stacking that "\n" would otherwise give.
  scale_fill_manual(
    values = c("Q (contracted)" = "#bdbdbd",
               "mu (expected)" = "#6baed6",
               "p5 (95% guaranteed)" = "#d6604d"),
    labels = c(expression(italic(Q)~"(contracted)"),
               expression(mu~"(expected)"),
               expression(italic(D)[5]~"(95% guaranteed)")),
    name = NULL) +
  # atop() rather than "\n" because the third label needs a subscript. Its line
  # gap is wider than "\n" would give and is not adjustable in base plotmath;
  # ggtext would fix it but is not worth a replication-package dependency.
  scale_x_discrete(labels = c(
    expression(atop("R1", "max volume")),
    expression(atop("R2", "max DACCS")),
    expression(atop("R3", "max"~italic(D)[5])))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
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
  bind_rows(r1, r2, r3) %>% mutate(rho_F = rho_F) %>%
    check_interior(paste0("rho_Forest sweep at ", rho_F))
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

# "Price-ratio variants for main-text panel (b)" ----------
# The threshold at which the CaR-constrained buyer leaves the all-forest corner
# is set by the price of the durable technology, not by correlation alone. Sweep
# rho_Forest again at half and twice the central DACCS price to show that.

price_variants <- list(
  list(label = "DACCS $250/tCO2", c_dacs = 250),
  list(label = "DACCS $500/tCO2", c_dacs = 500),
  list(label = "DACCS $1,000/tCO2", c_dacs = 1000)
)
CENTRAL_LABEL <- "DACCS $500/tCO2"

sweep_prices <- map_dfr(price_variants, function(v) {
  costs_v <- c(DACS = v$c_dacs, Forest = tonne_prices[["Forest"]]) * q
  map_dfr(rho_grid, function(rho_F) {
    portfolio_stats(
      df = grid, q = q, p = survival_probs, costs = costs_v,
      z = z_p5, rho_within = c(DACS = 0.5, Forest = rho_F), rho_between = 0
    ) %>%
      filter(p5 >= T_target) %>%
      slice_min(cost, with_ties = FALSE) %>%
      check_interior(paste0(v$label, ", rho_Forest = ", rho_F)) %>%
      mutate(rho_F = rho_F, variant = v$label)
  })
}) %>%
  mutate(variant = factor(variant, levels = vapply(price_variants,
                                                  function(v) v$label, "")))

thresholds <- sweep_prices %>%
  filter(forest_share < 1) %>%
  group_by(variant) %>%
  arrange(rho_F, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

cat("\n--- Diversification threshold by DACCS price ---\n")
thresholds %>%
  mutate(threshold = sprintf("%.3f", rho_F)) %>%
  select(variant, threshold, n_DACS, n_Forest) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat("\n")

# "Main-text panels" ----------
# The main text makes two points and needs two panels. Panel (a) is the
# rules-free version of (c): the feasibility boundary is the subject, since many
# portfolios hit the same permanence target and differ only in cost. Panel (b) is
# (e) restricted to the CaR-constrained optimum, since the three-rule comparison
# now lives in the SI.

cheapest <- hl %>% filter(rule == "R3: CaR")

pa_main <- ggplot() +
  rasterise(
    geom_point(data = frontier_corr %>% filter(!meets_p5, n_total > 0),
               aes(x = sd, y = mu), color = "grey88", size = 0.35, shape = 16),
    dpi = 600
  ) +
  rasterise(
    geom_point(data = frontier_corr %>% filter(meets_p5),
               aes(x = sd, y = mu, color = cost / 1e6), size = 0.55, shape = 16),
    dpi = 600
  ) +
  geom_line(data = p5_line, aes(x = sd, y = mu),
            color = "#e41a1c", linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 7800, y = T_target - z_p5 * 7800 - 3800,
           label = "italic(D)[5] == italic(G)", parse = TRUE,
           size = 2.3, color = "#e41a1c") +
  geom_point(data = cheapest, aes(x = sd, y = mu),
             shape = 21, size = 2.5, fill = "white", color = "black", stroke = 0.6) +
  geom_text_repel(data = cheapest,
                  aes(x = sd, y = mu,
                      label = paste0("cheapest feasible\n", n_DACS, "D + ", n_Forest, "F")),
                  size = 2, fontface = "bold", lineheight = 0.8,
                  min.segment.length = 0, nudge_y = 6000, nudge_x = 1200,
                  point.padding = unit(0.3, "lines")) +
  scale_color_viridis_c(name = "Cost ($m)", option = "D", direction = -1,
                        labels = function(x) paste0("$", x, "m")) +
  scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(10000, 60000)) +
  labs(x = expression("Standard deviation of delivery"~(tCO[2])),
       y = expression("Expected delivery"~mu~(tCO[2]))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm")) +
  guides(color = guide_colorbar(barwidth = 5, barheight = 0.3))

rho_thresh <- sweep %>%
  filter(rule == "R3: CaR") %>%
  arrange(rho_F) %>%
  filter(forest_share < 1) %>%
  slice_head(n = 1)

# Central case solid; the two price variants faint, to show the threshold moving.
# Light / mid / dark blue rather than three alphas of one colour, so the two
# price variants are distinguishable from each other as well as from the centre.
# ColorBrewer Blues 3 / 6 / 9: chosen for separation between all three, not just
# between the variants and the centre.
price_cols <- setNames(
  c("#C6DBEF", "#4292C6", "#084594"),
  vapply(price_variants, function(v) v$label, "")
)

pb_main <- ggplot(sweep_prices, aes(x = rho_F, y = forest_share,
                                    color = variant, linewidth = variant)) +
  geom_vline(xintercept = rho_thresh$rho_F, linetype = "dotted",
             color = "grey45", linewidth = 0.4) +
  geom_line() +
  # Label sits in the top margin, directly above the line it names.
  annotate("text", x = rho_thresh$rho_F, y = 1.06,
           label = sprintf("rho[Forest] == %.2f", rho_thresh$rho_F),
           parse = TRUE, hjust = 0.5, vjust = 0, size = 2.2, color = "grey30") +
  scale_color_manual(values = price_cols, name = NULL) +
  scale_linewidth_manual(values = c(0.45, 0.7, 0.45), name = NULL) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(breaks = seq(0, 0.4, 0.1)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = expression(rho[Forest]),
       y = "Forest share of cost-minimising\nportfolio meeting the target") +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 14, r = 6, b = 2, l = 2))

cat("Main-text panels done\n")

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

pf_r1 <- pf_data %>% filter(rule == "R1: Cheapest", !is.na(c_eff))
pf_rest <- pf_data %>% filter(rule != "R1: Cheapest")

# R1's guarantee collapses before the panel ends: mark where its fifth
# percentile reaches zero, beyond which c_eff is undefined.
rho_r1_zero <- pf_data %>%
  filter(rule == "R1: Cheapest", is.na(c_eff)) %>%
  pull(rho_F) %>% min()

pf <- ggplot() +
  geom_line(data = pf_rest, aes(x = rho_F, y = c_eff, color = rule),
            linewidth = 0.5) +
  geom_line(data = pf_r1, aes(x = rho_F, y = c_eff, color = rule),
            linewidth = 0.5) +
  geom_vline(xintercept = rho_r1_zero, linetype = "dotted",
             color = rule_cols[["R1: Cheapest"]], linewidth = 0.3) +
  annotate("text", x = rho_r1_zero, y = c_eff_cap * 0.72,
           label = "R1: italic(D)[5] == 0", parse = TRUE, angle = 90,
           vjust = -0.5, size = 1.9, color = rule_cols[["R1: Cheapest"]]) +
  scale_color_manual(values = rule_cols, name = NULL,
                     breaks = names(rule_cols)) +
  scale_x_continuous(breaks = seq(0, 0.4, 0.1)) +
  coord_cartesian(ylim = c(0, c_eff_cap)) +
  labs(x = expression(rho[Forest]),
       y = expression(c[eff*","*95]~"($/tCO"[2]*")")) +
  theme(legend.position = "bottom")

cat("Panel f done\n")

# "Assemble main figure" ----------

# `patchwork & theme` errors under ggplot2 4.x, so the tag theme is added to each
# panel individually rather than broadcast across the assembled object.
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 9),
                   plot.margin = margin(2, 6, 2, 2))

fig <- ((pa_main + tag_theme) | (pb_main + tag_theme)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(out_dir, "figure4.pdf"), fig,
       width = FIG_WIDTH_MM, height = FIG_WIDTH_MM * 0.52,
       units = "mm", bg = "white")
cat("Figure 4 saved to", file.path(out_dir, "figure4.pdf"), "\n")

# "SI figures: the panels displaced from the main text" ----------
# Rendered at the width they are displayed at in supplement.tex, so the panel
# base_size is the true on-page point size.
SI_LINEWIDTH_IN <- 500.484 / 72.27

si_panel <- function(p, name, frac, aspect) {
  ggsave(file.path(si_dir, name), p,
         width = frac * SI_LINEWIDTH_IN,
         height = frac * SI_LINEWIDTH_IN * aspect, bg = "white")
  cat("  ", name, "saved\n")
}

si_panel(pa, "si_portfolio_composition.pdf", 0.60, 0.80)
si_panel(pb, "si_portfolio_delivery.pdf",    0.60, 0.70)
si_panel(pd, "si_portfolio_budget_dual.pdf", 0.60, 0.85)
si_panel(pf, "si_portfolio_ceff_rho.pdf",    0.60, 0.85)

# "SI figures" ----------
cat("Generating SI portfolio figures...\n")

# SI: Min cost vs rho_forest for three rho_DACCS values
rho_dacs_grid <- c(0, 0.2, 0.5)
rho_forest_grid <- seq(0, 0.8, by = 0.01)

sweep_grid <- expand_grid(n_DACS = seq(0, N_DACS_MAX, by = 1),
                          n_Forest = seq(0, N_FOREST_MAX, by = 1))

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
      check_interior(paste0("cost-vs-rho sweep at rho_Forest = ", rho_within_forest,
                            ", rho_DACCS = ", rho_within_dacs)) %>%
      mutate(rho_within_forest = rho_within_forest,
             rho_within_dacs   = rho_within_dacs)
  })

rho_colors <- c("0" = "steelblue", "0.2" = "darkorange", "0.5" = "firebrick")

si_cost_rho <- sweep_results %>%
  ggplot(aes(x = rho_within_forest, y = cost,
             color = factor(rho_within_dacs),
             group = factor(rho_within_dacs))) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = rho_colors, name = expression(rho[DACCS])) +
  scale_y_continuous(labels = scales::dollar) +
  labs(x = expression(rho[Forest]), y = "Minimum cost") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(si_dir, "si_portfolio_cost_vs_rho.pdf"),
       si_cost_rho, width = 7, height = 5)
cat("  si_portfolio_cost_vs_rho.pdf saved\n")

# SI: Effective price heatmap
heat_max <- max(results_cm$n_DACS, results_cm$n_Forest) * 2.5
heat_grid <- expand_grid(n_DACS = seq(0, heat_max, by = 1), n_Forest = seq(0, heat_max, by = 1))

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

# The heatmap grid now matches the solver's, so each plotted cross must coincide
# with the fine-grid optimum the SI caption quotes. Assert it rather than assume:
# a mismatch means the caption is describing a point the figure does not show.
opt_indep_fine <- results_cm %>%
  filter(scenario == "Independent", rule == "R3: CaR")
opt_cal_fine <- results_cm %>%
  filter(scenario == "Correlated", rule == "R3: CaR")

check_marker <- function(marked, fine, what) {
  if (marked$n_DACS != fine$n_DACS || marked$n_Forest != fine$n_Forest) {
    stop("Heatmap cross disagrees with the quoted optimum under ", what, ": ",
         "plotted ", marked$n_DACS, "D+", marked$n_Forest, "F against ",
         fine$n_DACS, "D+", fine$n_Forest, "F. The SI caption would be wrong.")
  }
  invisible(TRUE)
}
check_marker(opt_indep, opt_indep_fine, "independence")
check_marker(opt_cal, opt_cal_fine, "calibrated correlations")

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
                         labels = scales::label_dollar(accuracy = 1)) +
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
# Every number the manuscript quotes from this figure, at full precision.

fmt <- function(x) formatC(x, format = "f", big.mark = ",", digits = 1)

cat("\n=========== NUMBERS QUOTED IN THE MANUSCRIPT ===========\n")
cat("Prices: $", tonne_prices["Forest"], "/tCO2 forest, $", tonne_prices["DACS"],
    "/tCO2 DACCS; q = ", q, " tCO2 per project\n", sep = "")
cat("Target G = ", format(T_target, big.mark = ","), " tCO2\n\n", sep = "")

cat("--- Cost-minimising rules ---\n")
results_cm %>%
  mutate(Q = n_total * q, c_eff = cost / p5) %>%
  select(scenario, rule, n_DACS, n_Forest, Q, mu, p5, cost, c_eff) %>%
  arrange(scenario, rule) %>%
  mutate(across(c(Q, mu, p5, cost), fmt), c_eff = round(c_eff, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n--- Budget dual (budget = $", format(B, big.mark = ","),
    ", each rule maximises its own objective) ---\n", sep = "")
results_budget %>%
  mutate(Q = n_total * q) %>%
  select(rule, n_DACS, n_Forest, Q, mu, p5, cost) %>%
  mutate(across(c(Q, mu, p5, cost), fmt)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n--- Rule 3 primal vs dual (should coincide exactly) ---\n")
pr <- results_cm %>% filter(scenario == "Correlated", rule == "R3: CaR")
du <- results_budget %>% filter(rule == "R3: CaR")
cat(sprintf("  primal: %dD + %dF, cost $%s, p5 = %.4f\n",
            pr$n_DACS, pr$n_Forest, format(pr$cost, big.mark = ","), pr$p5))
cat(sprintf("  dual:   %dD + %dF, cost $%s, p5 = %.4f\n",
            du$n_DACS, du$n_Forest, format(du$cost, big.mark = ","), du$p5))
cat(sprintf("  identical portfolio: %s\n",
            identical(c(pr$n_DACS, pr$n_Forest), c(du$n_DACS, du$n_Forest))))

cat("\n--- Diversification threshold (Rule 3 leaves the all-forest corner) ---\n")
thr <- sweep %>%
  filter(rule == "R3: CaR") %>%
  arrange(rho_F) %>%
  filter(forest_share < 1) %>%
  slice_head(n = 1)
cat(sprintf("  rho_Forest = %.3f  (first interior mix: %dD + %dF)\n",
            thr$rho_F, thr$n_DACS, thr$n_Forest))

cat("\n--- Effective cost per guaranteed tonne, c_eff95 ---\n")
sweep %>%
  filter(rho_F %in% c(0, 0.1, 0.2, 0.3, 0.4)) %>%
  mutate(c_eff = if_else(p5 > 0, cost / p5, NA_real_)) %>%
  select(rho_F, rule, c_eff) %>%
  pivot_wider(names_from = rule, values_from = c_eff) %>%
  mutate(across(-rho_F, ~round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n--- Optima marked on the SI heatmap (quote these in the caption) ---\n")
cat(sprintf("  independence: %dD + %dF, cost $%s\n",
            opt_indep$n_DACS, opt_indep$n_Forest,
            format(opt_indep$cost, big.mark = ",")))
cat(sprintf("  calibrated:   %dD + %dF, cost $%s\n",
            opt_cal$n_DACS, opt_cal$n_Forest,
            format(opt_cal$cost, big.mark = ",")))

cat("\n--- SI expected-value procurement rule (mu >= G, forest only) ---\n")
ev <- stats_corr %>% filter(n_DACS == 0, mu >= T_target) %>%
  slice_min(cost, with_ties = FALSE)
cat(sprintf("  %dD + %dF, cost $%s, mu = %s, p5 = %s\n",
            ev$n_DACS, ev$n_Forest,
            formatC(ev$cost, format = "d", big.mark = ","),
            fmt(ev$mu), fmt(ev$p5)))

cat("\n=======================================================\n")
cat("\nDone.\n")

theme_set(default_theme)
