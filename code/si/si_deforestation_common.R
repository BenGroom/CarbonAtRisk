# Shared helpers for the conversion-risk SI analysis.
#
# Sourced by si_deforestation_history.R (descriptive figure) and
# si_deforestation.R (model and results). Holds the region crosswalk, the two
# data loaders, and the parameters that must match figure2.R.

library(dplyr)

# Parameters pinned to figure2.R. RESCALE_FIRESIZE in particular: the default in
# simulate_burn_dynamics is TRUE, and silently taking it produced an SI figure
# ~15% adrift from the main text.
EFFIS_CACHE      <- "data/effis_cache"
DEFOR_PANEL      <- "data/conversion_rates.csv"
CLIMATE_RATE     <- 0.005      # gamma; applies to fire only
RESCALE_FIRESIZE <- FALSE
N_SIMULATIONS    <- 5000
DEFOR_YEARS      <- c(2001, 2025)

# The EFFIS cache key and the GFW (iso, adm1) pair refer to the same GADM unit.
# `adm1` is a bare integer in the panel, so the crosswalk is explicit rather than
# parsed out of the key.
REGIONS <- tibble::tribble(
  ~region,        ~effis_key,   ~iso,   ~adm1,  ~gid_1,
  "California",   "USA_5_1",    "USA",      5,  "USA.5_1",
  "Mato Grosso",  "BRA_12_1",   "BRA",     12,  "BRA.12_1",
  "Papua",        "IDN_23_1",   "IDN",     23,  "IDN.23_1"
)

#' Load the annual conversion-rate panel for the three focal regions.
#'
#' The panel is GFW `gadm__tcl__adm1_change` v20260424 with the WRI/Google DeepMind
#' 1 km dominant-driver attribution of Sims et al. (2025), restricted to five non-fire
#' anthropogenic drivers, so wildfire is excluded at source and the two hazards are
#' disjoint.
#'
#' The file ships with the repository as a three-region extract (USA/5, BRA/12,
#' IDN/23; 2001--2025, 75 rows). Columns are iso, adm1, year, delta, where delta is
#' non-fire anthropogenic loss in the year as a fraction of the unit's 2000 forest
#' area. `si_defor_driver_audit.py` reports the driver shares behind it.
#'
#' @return tibble(region, year, delta)
load_defor_panel <- function(path = DEFOR_PANEL) {
  if (!file.exists(path)) {
    stop("Conversion panel not found at ", path,
         ".\nIt ships with the repository; restore it from version control.")
  }
  raw <- readr::read_csv(path, show_col_types = FALSE)

  out <- REGIONS %>%
    select(region, iso, adm1) %>%
    inner_join(raw, by = c("iso", "adm1")) %>%
    filter(year >= DEFOR_YEARS[1], year <= DEFOR_YEARS[2]) %>%
    select(region, year, delta) %>%
    arrange(region, year)

  # Fail loudly rather than silently modelling two regions, since `adm1` is a bare
  # integer and a mis-join would look like plausible data.
  missing <- setdiff(REGIONS$region, unique(out$region))
  if (length(missing) > 0) {
    stop("No conversion data for: ", paste(missing, collapse = ", "))
  }
  counts <- out %>% count(region)
  if (any(counts$n < 20)) {
    stop("Fewer than 20 years of conversion data for: ",
         paste(counts$region[counts$n < 20], collapse = ", "))
  }
  if (any(out$delta < 0 | out$delta > 1, na.rm = TRUE) || anyNA(out$delta)) {
    stop("Conversion rates outside [0, 1] or missing.")
  }
  out
}

#' Load empirical annual burn fractions for the three focal regions.
#'
#' @return tibble(region, burn_fraction)
load_burn_fractions <- function(cache_dir = EFFIS_CACHE) {
  purrr::pmap_dfr(REGIONS, function(region, effis_key, iso, adm1, gid_1) {
    fires  <- read.csv(file.path(cache_dir, sprintf("effis_fire_%s.csv", effis_key)))
    forest <- readRDS(file.path(cache_dir, sprintf("effis_forest_%s.rds", effis_key)))
    bf <- calculate_burn_fractions(fires, as.numeric(forest$lc1),
                                   rescale_firesize = RESCALE_FIRESIZE)
    tibble::tibble(region = region, burn_fraction = as.numeric(bf))
  })
}

#' Regrowth rates keyed by region name, from the Zang calibration.
load_regrowth_rates <- function() {
  rates <- get_regrowth_rates()
  stats::setNames(as.numeric(rates[REGIONS$region]), REGIONS$region)
}
