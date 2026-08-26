# SI: forest CaR including conversion (deforestation) risk
#
# Extends the fire-and-regrowth model of Section 3.1 with a second, structurally
# different hazard: permanent conversion. Fire loss regrows; converted land does not.
#
# Three pools, U + B + C = 1:
#   U  intact standing carbon
#   B  burned, recoverable at rate r
#   C  converted, permanently lost
# Reported loss is unchanged, L_t = 1 - U_t, so the CaR definition carries over.
#
# Each year:
#   Step 1 (fire and regrowth, exactly as simulate_burn_dynamics):
#     U' = U - U*beta_t + B*r
#     B' = B + U*beta_t - B*r
#   Step 2 (conversion, pro rata across both pools since regrowing land is also
#           clearable):
#     U  = U'(1 - delta_t);  B = B'(1 - delta_t);  C = C + delta_t (U' + B')
#
# Because conversion scales U and B equally, the intact fraction of *surviving*
# land follows the unmodified fire recursion, giving the closed form
#
#     L_t = 1 - (1 - delta_bar)^t * (1 - L*_fire),   L*_fire = beta_bar/(r + beta_bar)
#
# Fire converges to a bounded equilibrium; conversion does not converge at all.
#
# Two specifications, with identical expected loss but very different tails:
#   "gradual"   a fraction delta_t of remaining area is cleared each year.
#               The SI's specification.
#   "absorbing" the project is cleared entirely at tau ~ Geometric(delta_bar).
#               Reported once as a bound.
#
# Data: GFW gadm__tcl__adm1_change v20260424 with the WRI/Google DeepMind 1 km
# dominant-driver attribution of Sims et al. (2025), restricted to five non-fire
# anthropogenic drivers. Wildfire is
# excluded at source, so fire and conversion are disjoint and can be composed
# without double-counting. See si_defor_driver_audit.py for the driver shares.
#
# Run from repo root:  Rscript code/si/si_deforestation.R

library(tidyverse)
library(patchwork)
library(scales)

source("code/0_funcs/fire_funcs.R")
source("code/0_funcs/regrowth_funcs.R")
source("code/si/si_deforestation_common.R")

FIG_WIDTH_MM <- 183
MAX_HORIZON  <- 200
KG_PER_TONNE <- 1000

out_dir <- "outputs/si"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# "Simulator" ------------------------------------------------------------------

#' Simulate cumulative carbon loss under fire, regrowth and permanent conversion.
#'
#' @param burn_fractions Empirical annual burn fractions (EFFIS).
#' @param defor_rates Empirical annual conversion rates (GFW non-fire drivers).
#'   Pass 0 to recover the fire-only model exactly.
#' @param regrowth_rate Annual regrowth rate r.
#' @param climate_rate Annual climate trend gamma; applies to fire only.
#' @param spec "gradual" (fraction cleared each year) or "absorbing"
#'   (whole project cleared at a geometric date).
#' @param max_horizon,n_simulations
#' @param return_pools If TRUE, also return the terminal U/B/C decomposition.
#' @return List with `loss` (n_simulations x (max_horizon+1) cumulative loss,
#'   fraction of original stock) and, optionally, `pools`.
simulate_loss_with_conversion <- function(burn_fractions,
                                          defor_rates,
                                          regrowth_rate,
                                          climate_rate = CLIMATE_RATE,
                                          spec = c("gradual", "absorbing"),
                                          max_horizon = MAX_HORIZON,
                                          n_simulations = N_SIMULATIONS,
                                          return_pools = FALSE) {
  spec <- match.arg(spec)

  loss <- matrix(0.0, nrow = n_simulations, ncol = max_horizon + 1)
  U <- rep(1.0, n_simulations)
  B <- rep(0.0, n_simulations)
  C <- rep(0.0, n_simulations)

  delta_bar <- mean(defor_rates)
  if (spec == "absorbing") {
    # Geometric with support {1, 2, ...}: rgeom counts failures before the first
    # success, so add 1 to get the year of clearing. tau = Inf when delta_bar = 0.
    tau <- if (delta_bar > 0) rgeom(n_simulations, delta_bar) + 1 else rep(Inf, n_simulations)
  }

  for (year in seq_len(max_horizon)) {
    beta <- pmin(sample(burn_fractions, n_simulations, replace = TRUE) *
                   (1 + climate_rate * (year - 1)), 1.0)

    # Step 1: fire and regrowth
    new_burn <- U * beta
    regrow   <- B * regrowth_rate
    U <- U - new_burn + regrow
    B <- B + new_burn - regrow

    # Step 2: conversion
    if (spec == "gradual") {
      delta <- if (delta_bar > 0) {
        sample(defor_rates, n_simulations, replace = TRUE)
      } else {
        rep(0.0, n_simulations)
      }
      C <- C + delta * (U + B)
      U <- U * (1 - delta)
      B <- B * (1 - delta)
    } else {
      hit <- year >= tau
      if (any(hit)) {
        C[hit] <- C[hit] + U[hit] + B[hit]
        U[hit] <- 0
        B[hit] <- 0
      }
    }

    U <- pmin(pmax(U, 0), 1)
    B <- pmin(pmax(B, 0), 1)
    C <- pmin(pmax(C, 0), 1)

    loss[, year + 1] <- 1 - U
  }

  out <- list(loss = loss)
  if (return_pools) out$pools <- list(U = U, B = B, C = C)
  out
}

#' Analytic loss under constant hazards and no climate trend.
closed_form_loss <- function(t, beta_bar, r, delta_bar) {
  L_fire <- beta_bar / (r + beta_bar)
  1 - (1 - delta_bar)^t * (1 - L_fire)
}

# "Validation" -----------------------------------------------------------------
# No results are reported until all six checks pass.

run_validation <- function(defor, burn, regrowth) {
  cat("=== VALIDATION ===\n")
  ok <- TRUE
  chk <- function(label, pass, detail = "") {
    cat(sprintf("  [%s] %s%s\n", if (pass) "PASS" else "FAIL", label,
                if (nzchar(detail)) paste0("  ", detail) else ""))
    pass
  }

  bf <- burn %>% filter(region == "California") %>% pull(burn_fraction)
  dr <- defor %>% filter(region == "California") %>% pull(delta)
  r  <- regrowth[["California"]]

  # 1. Pool invariant U + B + C = 1
  set.seed(101)
  s <- simulate_loss_with_conversion(bf, dr, r, spec = "gradual",
                                     max_horizon = 50, n_simulations = 200,
                                     return_pools = TRUE)
  tot <- s$pools$U + s$pools$B + s$pools$C
  ok <- chk("pool invariant U+B+C=1 (gradual)",
            max(abs(tot - 1)) < 1e-10,
            sprintf("max dev %.2e", max(abs(tot - 1)))) && ok

  set.seed(101)
  s2 <- simulate_loss_with_conversion(bf, dr, r, spec = "absorbing",
                                      max_horizon = 50, n_simulations = 200,
                                      return_pools = TRUE)
  tot2 <- s2$pools$U + s2$pools$B + s2$pools$C
  ok <- chk("pool invariant U+B+C=1 (absorbing)",
            max(abs(tot2 - 1)) < 1e-10,
            sprintf("max dev %.2e", max(abs(tot2 - 1)))) && ok

  # 2. Nesting: delta = 0 must reproduce the existing fire model bit-for-bit
  set.seed(202)
  a <- simulate_loss_with_conversion(bf, 0, r, spec = "gradual",
                                     max_horizon = 100, n_simulations = 300)$loss
  set.seed(202)
  b <- simulate_burn_dynamics(bf, r, CLIMATE_RATE, max_horizon = 100,
                              n_projects = 1, n_simulations = 300, correlation = 0)
  ok <- chk("delta=0 reproduces simulate_burn_dynamics exactly",
            identical(all.equal(a, b, tolerance = 0), TRUE),
            sprintf("max abs diff %.2e", max(abs(a - b)))) && ok

  # 3. Multiplicative separation. Because conversion scales U and B equally, the
  # intact fraction of surviving land follows the unmodified fire recursion:
  #
  #     U_t(fire+conversion) = (1 - delta_bar)^t * U_t(fire only)
  #
  # exactly, at every t. With constant beta and gamma = 0 both runs are
  # deterministic, so this should hold to machine precision. The equilibrium form
  # L* = beta_bar/(r + beta_bar) is only the t -> infinity limit of U_t(fire only),
  # and at t = 100 the fire system has not yet converged, so it is NOT the right
  # target here.
  beta_bar <- mean(bf); delta_bar <- mean(dr)
  set.seed(303)
  both <- simulate_loss_with_conversion(rep(beta_bar, 100), rep(delta_bar, 100), r,
                                        climate_rate = 0, spec = "gradual",
                                        max_horizon = 200, n_simulations = 20)$loss
  set.seed(303)
  fire_only <- simulate_loss_with_conversion(rep(beta_bar, 100), 0, r,
                                             climate_rate = 0, spec = "gradual",
                                             max_horizon = 200, n_simulations = 20)$loss
  horizons <- 0:200
  U_both <- 1 - colMeans(both)
  U_fire <- 1 - colMeans(fire_only)
  sep_err <- max(abs(U_both - (1 - delta_bar)^horizons * U_fire))
  ok <- chk("multiplicative separation U_both = (1-delta)^t * U_fire",
            sep_err < 1e-12,
            sprintf("max dev %.2e; at t=100 loss %.1f kg/t vs asymptotic %.1f",
                    sep_err, KG_PER_TONNE * mean(both[, 101]),
                    KG_PER_TONNE * closed_form_loss(100, beta_bar, r, delta_bar))) && ok

  # 4. Mean equivalence of the two specifications
  set.seed(404)
  g <- simulate_loss_with_conversion(bf, dr, r, spec = "gradual",
                                     n_simulations = 6000)$loss
  set.seed(404)
  ab <- simulate_loss_with_conversion(bf, dr, r, spec = "absorbing",
                                      n_simulations = 6000)$loss
  d100 <- abs(mean(g[, 101]) - mean(ab[, 101]))
  ok <- chk("specifications A and B share a mean at t=100",
            d100 < 0.02,
            sprintf("gradual %.1f vs absorbing %.1f kg/t",
                    KG_PER_TONNE * mean(g[, 101]), KG_PER_TONNE * mean(ab[, 101]))) && ok

  # 5. Absorbing saturation horizon: CaR95 hits total loss once P(cleared) > 5%
  t_star <- ceiling(log(0.95) / log(1 - delta_bar))
  car95_at <- function(m, t) quantile(m[, t + 1], 0.95, names = FALSE)
  ok <- chk("absorbing CaR95 saturates at t* (California)",
            car95_at(ab, t_star + 5) > 0.999,
            sprintf("t* = %d yr, CaR95 at t*+5 = %.1f kg/t",
                    t_star, KG_PER_TONNE * car95_at(ab, t_star + 5))) && ok

  # 6. Fire-only arm reproduces the Stage 1 CaR95 figures
  targets <- c(California = 643, `Mato Grosso` = 507, Papua = 48)
  devs <- purrr::map_dbl(names(targets), function(rg) {
    set.seed(101)
    m <- simulate_loss_with_conversion(
      burn %>% filter(region == rg) %>% pull(burn_fraction),
      0, regrowth[[rg]], spec = "gradual", max_horizon = 100)$loss
    KG_PER_TONNE * quantile(m[, 101], 0.95, names = FALSE) - targets[[rg]]
  })
  ok <- chk("fire-only CaR95 at 100y matches Stage 1",
            max(abs(devs)) < 15,
            sprintf("deviations: %s kg/t",
                    paste(sprintf("%+.0f", devs), collapse = ", "))) && ok

  cat("==================\n\n")
  if (!ok) stop("Validation failed; no results written.")
  invisible(TRUE)
}

# "Run" ------------------------------------------------------------------------

defor    <- load_defor_panel()
burn     <- load_burn_fractions()
regrowth <- load_regrowth_rates()

run_validation(defor, burn, regrowth)

WINDOWS <- list(
  "Full record (2001-25)" = c(2001, 2025),
  "2001-08"               = c(2001, 2008),
  "2009-16"               = c(2009, 2016),
  "2017-25"               = c(2017, 2025)
)

summarise_run <- function(loss, region, spec, window) {
  horizons <- 0:MAX_HORIZON
  tibble(
    region = region, spec = spec, window = window, horizon = horizons,
    mean_loss = KG_PER_TONNE * colMeans(loss),
    car_95    = KG_PER_TONNE * apply(loss, 2, quantile, 0.95, names = FALSE),
    car_99    = KG_PER_TONNE * apply(loss, 2, quantile, 0.99, names = FALSE)
  )
}

results <- purrr::pmap_dfr(REGIONS, function(region, effis_key, iso, adm1, gid_1) {
  bf <- burn  %>% filter(region == !!region) %>% pull(burn_fraction)
  dr <- defor %>% filter(region == !!region) %>% pull(delta)
  r  <- regrowth[[region]]

  # Fire only
  set.seed(101)
  fire <- summarise_run(
    simulate_loss_with_conversion(bf, 0, r, spec = "gradual")$loss,
    region, "fire only", "Full record (2001-25)")

  # Gradual, across resampling windows
  grad <- purrr::imap_dfr(WINDOWS, function(yrs, nm) {
    sub <- defor %>% filter(region == !!region, year >= yrs[1], year <= yrs[2]) %>% pull(delta)
    set.seed(101)
    summarise_run(simulate_loss_with_conversion(bf, sub, r, spec = "gradual")$loss,
                  region, "gradual", nm)
  })

  # Absorbing, full record only
  set.seed(101)
  absb <- summarise_run(
    simulate_loss_with_conversion(bf, dr, r, spec = "absorbing")$loss,
    region, "absorbing", "Full record (2001-25)")

  bind_rows(fire, grad, absb)
}) %>%
  mutate(region = factor(region, levels = REGIONS$region))

readr::write_csv(results, file.path(out_dir, "si_deforestation_results.csv"))
cat("Wrote", file.path(out_dir, "si_deforestation_results.csv"), "\n\n")

# "Figure" ---------------------------------------------------------------------

region_cols <- c(California = "#E67E22", `Mato Grosso` = "#8E44AD", Papua = "#16A085")
base_full <- results %>% filter(window == "Full record (2001-25)")

# (a) fire only vs fire + conversion
pa_dat <- base_full %>% filter(spec %in% c("fire only", "gradual"), horizon >= 1)
pa <- ggplot(pa_dat, aes(x = horizon, y = car_95, colour = region,
                         linetype = spec)) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = region_cols, name = NULL) +
  scale_linetype_manual(values = c("fire only" = "dashed", "gradual" = "solid"),
                        labels = c("Fire only", "Fire + conversion"),
                        name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
  labs(x = "Horizon (years)",
       y = expression(CaR[95]~"(kg per tonne contracted)")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.01, "cm"))

# (b) decomposition of CaR95 at 100 years
pb_dat <- base_full %>%
  filter(horizon == 100, spec %in% c("fire only", "gradual")) %>%
  select(region, spec, car_95) %>%
  pivot_wider(names_from = spec, values_from = car_95) %>%
  mutate(Fire = `fire only`, Conversion = gradual - `fire only`,
         ratio = gradual / `fire only`) %>%
  select(region, Fire, Conversion, ratio)

pb <- pb_dat %>%
  pivot_longer(c(Fire, Conversion), names_to = "source", values_to = "kg") %>%
  mutate(source = factor(source, levels = c("Conversion", "Fire"))) %>%
  ggplot(aes(x = region, y = kg, fill = source)) +
  geom_col(width = 0.62) +
  geom_text(data = pb_dat, inherit.aes = FALSE,
            aes(x = region, y = Fire + Conversion,
                label = sprintf("%.1f×", ratio)),
            vjust = -0.5, size = 2.4, fontface = "bold") +
  scale_fill_manual(values = c(Fire = "#E67E22", Conversion = "#8C6D46"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = expression(CaR[95]~"at 100 years (kg per tonne contracted)")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# (c) robustness: gradual vs absorbing
pc_dat <- base_full %>% filter(spec %in% c("gradual", "absorbing"), horizon >= 1)
pc <- ggplot(pc_dat, aes(x = horizon, y = car_95, colour = region,
                         linetype = spec)) +
  geom_line(linewidth = 0.6) +
  # Line styles annotated in place rather than via a legend, for the same reason
  # as panel b. Both labels sit clear of the curves.
  annotate("text", x = 118, y = 1050, label = "Absorbing (bound)",
           size = 2.1, colour = "grey30", hjust = 0.5) +
  annotate("text", x = 150, y = 700, label = "Gradual", size = 2.1,
           colour = "grey30", hjust = 0.5) +
  scale_colour_manual(values = region_cols, name = NULL, guide = "none") +
  scale_linetype_manual(values = c("gradual" = "solid", "absorbing" = "dotted"),
                        name = NULL, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, NA)) +
  labs(x = "Horizon (years)",
       y = expression(CaR[95]~"(kg per tonne contracted)")) +
  theme_classic(base_size = 9)

tag_theme <- theme(plot.tag = element_text(face = "bold", size = 10))
fig <- ((pa + tag_theme) | (pb + tag_theme) | (pc + tag_theme)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(out_dir, "si_deforestation.pdf"), fig,
       width = FIG_WIDTH_MM, height = FIG_WIDTH_MM * 0.44,
       units = "mm", bg = "white")
cat("Saved", file.path(out_dir, "si_deforestation.pdf"), "\n\n")

# "Reportable numbers" ---------------------------------------------------------

cat("CaR95, kg per tonne, full-record conversion rates:\n")
base_full %>%
  filter(horizon %in% c(30, 100, 200), spec %in% c("fire only", "gradual", "absorbing")) %>%
  select(region, spec, horizon, car_95) %>%
  mutate(car_95 = round(car_95, 0)) %>%
  pivot_wider(names_from = horizon, values_from = car_95, names_prefix = "yr") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nUnderstatement from omitting conversion (CaR95 ratio, gradual / fire only):\n")
base_full %>%
  filter(horizon %in% c(100, 200), spec %in% c("fire only", "gradual")) %>%
  select(region, spec, horizon, car_95) %>%
  pivot_wider(names_from = spec, values_from = car_95) %>%
  mutate(ratio = round(gradual / `fire only`, 2),
         across(c(`fire only`, gradual), ~round(.x, 0))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nResampling-window sensitivity (gradual, CaR95 at 100 years):\n")
results %>%
  filter(spec == "gradual", horizon == 100) %>%
  select(region, window, car_95) %>%
  mutate(car_95 = round(car_95, 0)) %>%
  pivot_wider(names_from = window, values_from = car_95) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nRegional ordering on CaR95 at 100 years:\n")
base_full %>%
  filter(horizon == 100, spec %in% c("fire only", "gradual")) %>%
  select(region, spec, car_95) %>%
  arrange(spec, car_95) %>%
  group_by(spec) %>%
  summarise(order = paste(region, collapse = " < "), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
