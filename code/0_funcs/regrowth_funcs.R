#' Regrowth rate calibration from Zang et al. (2024)
#'
#' Computes post-fire regrowth rates for the CaR simulation using
#' the parametric height-recovery and allometric equations from:
#'
#'   Zang, Y. et al. (2024). "A global dataset of forest regrowth
#'   following wildfires." Scientific Data, 11:1052.
#'   doi:10.1038/s41597-024-03891-5
#'
#' Zang et al. fitted satellite-derived canopy height recovery curves
#' for 8 macro-regions globally. For each region, two equations are
#' published (their Tables 1 and 2):
#'
#'   Height recovery:   H(t) = a * exp(b * t) + c
#'   Height-to-biomass: AGB(H) = alpha * H^beta
#'
#' where t is years since fire, b < 0 (so H rises toward asymptote c),
#' and AGB is in Mg per 900 m^2 (30m x 30m pixel).
#'
#' We convert these to an exponential regrowth rate r by:
#'   1. Computing AGB(t) at several post-fire ages
#'   2. Computing the recovery ratio: AGB(t) / AGB(infinity)
#'   3. Inverting: r(t) = -ln(1 - ratio(t)) / t
#'   4. Averaging r over ages 10-30 years
#'
#' California override: Zang's "North America" macro-region yields
#' r = 8.3%/yr, but this continental average is dominated by fast-
#' regrowing eastern temperate and boreal forests. For California's
#' Mediterranean/montane ecosystems, independent estimates consistently
#' give r = 1.5-3%/yr:
#'   - Cook-Patton et al. (2020): median 1.01 MgC/ha/yr for California,
#'     implying r = 1.5-2.5% depending on mature-forest stock
#'   - Post-fire western US conifer recovery half-times of 25-50 years
#'     (Stevens-Rumann et al. 2018), implying r = 1.4-2.8%/yr
#' We set California r = 2.0%/yr as the midpoint of these estimates.
#' Sensitivity to this choice is reported in Supplementary Section X.

# "Zang et al. (2024) Table 1 parameters" ----------

#' Parameters for the three macro-regions corresponding to our study sites:
#'   - North America  -> California  (montane cordillera)
#'   - South America  -> Mato Grosso (broadleaf evergreen)
#'   - Southeast Asia  -> Papua      (broadleaf evergreen)
#'
#' Height recovery: H(t) = a * exp(b*t) + c   [metres]
#' Allometry:       AGB(H) = alpha * H^beta    [Mg per 900 m^2]

ZANG_PARAMS <- data.frame(
  macro_region = c("North America", "South America", "Southeast Asia"),
  admin_unit   = c("California",    "Mato Grosso",   "Papua"),
  a     = c(-10.70, -24.34, -27.84),
  b     = c(-0.0885, -0.0477, -0.0340),
  c_ht  = c(14.23, 28.03, 30.75),
  alpha = c(0.1372, 0.0541, 0.0194),
  beta  = c(1.5512, 1.894, 2.1604),
  stringsAsFactors = FALSE
)

# "Compute Zang regrowth rate" ----------

#' Compute the implied exponential regrowth rate from Zang et al.
#' parametric equations for a single macro-region.
#'
#' @param a Height model intercept offset (negative)
#' @param b Height model decay rate (negative)
#' @param c_ht Height model asymptote (mature canopy height, metres)
#' @param alpha Allometric coefficient
#' @param beta Allometric exponent
#' @param ages Numeric vector of post-fire ages to evaluate (years)
#' @return Scalar: average implied r across the specified ages
compute_zang_regrowth_rate <- function(a, b, c_ht, alpha, beta,
                                       ages = c(10, 15, 20, 30)) {
  # Asymptotic (mature-forest) AGB
  sbar <- alpha * c_ht^beta

  # AGB at each age
  h_t   <- a * exp(b * ages) + c_ht
  agb_t <- alpha * h_t^beta

  # Recovery ratio (fraction of mature AGB)
  ratio <- agb_t / sbar

  # Invert exponential model: S(t)/Sbar = 1 - exp(-r*t)
  # => r = -ln(1 - ratio) / t
  r_implied <- -log(1 - ratio) / ages

  # Average across ages
  mean(r_implied)
}

# "Get regrowth rates" ----------

#' Return calibrated regrowth rates for California, Mato Grosso, Papua.
#'
#' Computes Zang parametric rates transparently from Table 1 parameters,
#' then applies the California override (see header documentation).
#'
#' @param ca_override Regrowth rate for California (default 0.020 = 2.0%/yr).
#'   Set to NULL to use the raw Zang continental average.
#' @return Named numeric vector of annual regrowth rates (as fractions)
get_regrowth_rates <- function(ca_override = 0.020) {

  rates <- setNames(numeric(nrow(ZANG_PARAMS)), ZANG_PARAMS$admin_unit)

  for (i in seq_len(nrow(ZANG_PARAMS))) {
    p <- ZANG_PARAMS[i, ]
    rates[p$admin_unit] <- compute_zang_regrowth_rate(
      a = p$a, b = p$b, c_ht = p$c_ht,
      alpha = p$alpha, beta = p$beta
    )
  }

  # Report raw Zang values
  cat("=== Zang parametric regrowth rates (raw) ===\n")
  for (nm in names(rates)) {
    cat(sprintf("  %s: %.1f%%/yr\n", nm, rates[nm] * 100))
  }

  # Apply California override
  if (!is.null(ca_override)) {
    cat(sprintf(
      "\nCalifornia override: %.1f%% -> %.1f%% (literature-adjusted)\n",
      rates["California"] * 100, ca_override * 100
    ))
    rates["California"] <- ca_override
  }

  cat("\n=== Final regrowth rates ===\n")
  for (nm in names(rates)) {
    cat(sprintf("  %s: %.1f%%/yr\n", nm, rates[nm] * 100))
  }

  rates
}
