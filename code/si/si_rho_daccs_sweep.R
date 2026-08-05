# SI figure: sensitivity of the portfolio result to the within-DACCS correlation
#
# Section 4 sets rho_DACCS = 0.5 as an illustrative value, since there are no
# observed geological-storage reversals from which to estimate one. This figure
# shows the results across its full range. Two findings:
#
#   (a) Rule 3 remains an interior mix at every value including zero, so the
#       diversification result does not depend on the choice. What drives
#       diversification is rho_Forest, not rho_DACCS.
#   (b) Rule 2 misses the delivery guarantee at every value including zero, so
#       its shortfall is not an artefact of the correlation assumption either.
#
# Run from repo root:  Rscript code/si/si_rho_daccs_sweep.R

library(tidyverse)
library(patchwork)
library(scales)

source("code/0_funcs/portfolio_funcs.R")

FIG_WIDTH_MM <- 183

# Parameters must match figure4.R ---------------------------------------------
q         <- 500
z_p5      <- qnorm(0.05)
T_target  <- 30000
tonne_prices   <- c(DACS = 500, Forest = 40)
project_costs  <- tonne_prices * q
survival_probs <- c(DACS = 0.99, Forest = 0.4)
RHO_FOREST     <- 0.2

N_DACS_MAX   <- 200
N_FOREST_MAX <- 2000

grid <- expand_grid(n_DACS   = seq(0, N_DACS_MAX, by = 1),
                    n_Forest = seq(0, N_FOREST_MAX, by = 1))

out_dir <- "outputs/si"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rule_cols <- c("R2: Permanent" = "#ff7f00", "R3: CaR" = "#377eb8")

# Sweep ------------------------------------------------------------------------
rho_grid <- seq(0, 0.9, by = 0.025)

sweep <- map_dfr(rho_grid, function(rho_D) {
  st <- portfolio_stats(
    df = grid, q = q, p = survival_probs, costs = project_costs,
    z = z_p5, rho_within = c(DACS = rho_D, Forest = RHO_FOREST),
    rho_between = 0
  ) %>% mutate(Q = n_total * q)

  r2 <- st %>% filter(n_Forest == 0, Q >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R2: Permanent")
  r3 <- st %>% filter(p5 >= T_target) %>%
    slice_min(cost, with_ties = FALSE) %>% mutate(rule = "R3: CaR")

  out <- bind_rows(r2, r3) %>% mutate(rho_D = rho_D)
  if (any(out$n_DACS >= N_DACS_MAX | out$n_Forest >= N_FOREST_MAX)) {
    stop("Grid boundary binding at rho_DACCS = ", rho_D)
  }
  out
}) %>% mutate(rule = factor(rule, levels = names(rule_cols)))

# Panel a: Rule 3 composition, by technology -----------------------------------
#
# Plotted as project counts rather than the forest *share*. The share is a ratio
# of two integers whose costs trade at 12.5:1, so it sawtooths: the optimum
# alternates between a 61-DACCS family and a 62-DACCS family, and one extra DACCS
# project displaces about twelve forest projects. Showing both counts makes the
# alternation self-explanatory instead of looking like noise, and it carries the
# actual finding more directly: the DACCS holding is essentially pinned while
# forest does all of the adjusting. Total cost, by contrast, is smooth and
# monotone across the whole range (+6.1%).
comp <- sweep %>%
  filter(rule == "R3: CaR") %>%
  select(rho_D, DACCS = n_DACS, Forest = n_Forest) %>%
  pivot_longer(c(DACCS, Forest), names_to = "tech", values_to = "n") %>%
  mutate(tech = factor(tech, levels = c("DACCS", "Forest")))

# Free y-scales: on a shared axis the 61-to-62 DACCS step is invisible, and the
# forest dips then read as noise rather than as its consequence.
pa <- ggplot(comp, aes(x = rho_D, y = n, colour = tech)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.8) +
  facet_wrap(~tech, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_colour_manual(values = c(DACCS = "#984ea3", Forest = "#4daf4a"),
                      name = NULL, guide = "none") +
  # Project counts are integers; the default breaks give 61.25, 61.50, ...
  scale_y_continuous(breaks = function(lims) {
    unique(round(pretty(lims, n = 4)))
  }) +
  labs(x = expression(rho[DACCS]), y = NULL) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold"))

# Panel b: fifth percentile of delivery, against the target --------------------
pb <- sweep %>%
  ggplot(aes(x = rho_D, y = p5, colour = rule)) +
  geom_hline(yintercept = T_target, linetype = "dashed", colour = "grey40",
             linewidth = 0.4) +
  # R3 binds exactly on the obligation, so the label sits below the line to
  # avoid being overplotted by it.
  annotate("text", x = 0.02, y = T_target, label = "delivery obligation G",
           vjust = 1.7, hjust = 0, size = 2.2, colour = "grey30") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 0.9) +
  scale_colour_manual(values = rule_cols, name = NULL) +
  scale_y_continuous(labels = comma) +
  labs(x = expression(rho[DACCS]), y = expression(D[5]~(tCO[2]))) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom")

tag_theme <- theme(plot.tag = element_text(face = "bold", size = 10))
fig <- ((pa + tag_theme) | (pb + tag_theme)) + plot_annotation(tag_levels = "a")

ggsave(file.path(out_dir, "si_rho_daccs_sweep.pdf"), fig,
       width = FIG_WIDTH_MM, height = FIG_WIDTH_MM * 0.45,
       units = "mm", bg = "white")
cat("Saved", file.path(out_dir, "si_rho_daccs_sweep.pdf"), "\n\n")

# Reportable numbers -----------------------------------------------------------
cat("rho_Forest held at", RHO_FOREST, "\n\n")
sweep %>%
  filter(rho_D %in% c(0, 0.2, 0.5, 0.9)) %>%
  mutate(meets = p5 >= T_target) %>%
  select(rho_D, rule, n_DACS, n_Forest, forest_share, cost, p5, meets) %>%
  mutate(forest_share = round(forest_share, 3),
         p5 = round(p5, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

r3 <- sweep %>% filter(rule == "R3: CaR")
r2 <- sweep %>% filter(rule == "R2: Permanent")
cat("\nRule 3 interior (forest share < 1) at every rho_DACCS: ",
    all(r3$forest_share < 1), "\n", sep = "")
cat("Rule 3 holds forest at every rho_DACCS:                 ",
    all(r3$n_Forest > 0), "\n", sep = "")
cat("Rule 2 misses the guarantee at every rho_DACCS:         ",
    all(r2$p5 < T_target), "\n", sep = "")
cat("Rule 3 forest share range: ",
    sprintf("%.0f%% to %.0f%%", 100 * min(r3$forest_share),
            100 * max(r3$forest_share)), "\n", sep = "")
cat("Rule 3 cost range:         ",
    sprintf("$%s to $%s (+%.1f%%)",
            formatC(min(r3$cost), format = "d", big.mark = ","),
            formatC(max(r3$cost), format = "d", big.mark = ","),
            100 * (max(r3$cost) / min(r3$cost) - 1)), "\n", sep = "")
