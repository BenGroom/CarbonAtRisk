# SI figure: the three phases of the CaR curve, made legible
#
# Main-text Fig 2b plots the CaR curves without marking the regimes, so the three
# phases are not distinguishable. This figure shades the regimes derived
# analytically in the SI (transient 0-30 yr, quasi-equilibrium 30-80 yr, climate
# drift 80+ yr) and overlays the moving equilibrium the system is chasing.
#
# Parameters match figure2.R exactly, including RESCALE_FIRESIZE = FALSE.
#
# Run from repo root:  Rscript code/si/si_car_phases.R

library(tidyverse)
library(patchwork)
library(scales)

source("code/0_funcs/fire_funcs.R")
source("code/0_funcs/regrowth_funcs.R")

# Parameters must match figure2.R ---------------------------------------------
EFFIS_CACHE      <- "data/effis_cache"
PROJECT_AREA     <- 100000
CLIMATE_RATE     <- 0.005
RESCALE_FIRESIZE <- FALSE
N_SIMULATIONS    <- 5000
REGROWTH_RATES   <- get_regrowth_rates()

# Regime boundaries from the SI analytic characterisation
PHASE_BREAKS <- c(0, 30, 80, 200)
PHASE_LABELS <- c("Transient", "Quasi-equilibrium", "Climate drift")
PHASE_FILLS  <- c("#fde0dd", "#e5f5e0", "#deebf7")

REGIONS <- tribble(
  ~key,        ~region,       ~colour,
  "USA_5_1",   "California",  "#E67E22",
  "BRA_12_1",  "Mato Grosso", "#8E44AD",
  "IDN_23_1",  "Papua",       "#16A085"
)
REGROWTH_KEY <- c(California = "USA_5_1", `Mato Grosso` = "BRA_12_1",
                  Papua = "IDN_23_1")

out_dir <- "outputs/si"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(CAR_SEED)

# Simulate ---------------------------------------------------------------------
car_curves <- pmap_dfr(REGIONS, function(key, region, colour) {
  fires  <- read.csv(file.path(EFFIS_CACHE, sprintf("effis_fire_%s.csv", key)))
  forest <- readRDS(file.path(EFFIS_CACHE, sprintf("effis_forest_%s.rds", key)))

  bf <- calculate_burn_fractions(fires, as.numeric(forest$lc1),
                                 project_area = PROJECT_AREA,
                                 rescale_firesize = RESCALE_FIRESIZE)

  r <- REGROWTH_RATES[[region]]

  horizons <- 1:200
  set.seed(CAR_SEED)
  res <- simulate_diversified_car(
    bf, r, CLIMATE_RATE,
    years = horizons, n_simulations = N_SIMULATIONS,
    n_projects = 1, estate_area = 1000, correlation = 0
  )

  # simulate_diversified_car returns only mean_loss and car_95, in the order of
  # `years`; there is no year column to read back.
  tibble(region = region, year = horizons,
         car_95 = res$car_95, mean_loss = res$mean_loss,
         beta_bar = mean(bf), r = r)
}) %>%
  mutate(region = factor(region, levels = REGIONS$region))

region_cols <- setNames(REGIONS$colour, REGIONS$region)

# Reference: the MOVING equilibrium, E[L*_t] = bbar(1+gamma(t-1)) / (r + bbar(1+gamma(t-1))).
#
# The static gamma = 0 fixed point is the wrong reference for these curves: they
# are simulated with gamma = 0.005, so they correctly run away from it and the
# line reads as unexplained. The moving equilibrium is the target the system is
# actually chasing, and the mean loss should track it once the transient is over.
equil_static <- car_curves %>%
  distinct(region, beta_bar, r) %>%
  mutate(equilibrium = 1000 * beta_bar / (r + beta_bar))

equil_moving <- car_curves %>%
  mutate(bb_t = beta_bar * (1 + CLIMATE_RATE * (year - 1)),
         moving_eq = 1000 * bb_t / (r + bb_t)) %>%
  select(region, year, moving_eq)

bands <- tibble(
  xmin  = head(PHASE_BREAKS, -1),
  xmax  = tail(PHASE_BREAKS, -1),
  phase = factor(PHASE_LABELS, levels = PHASE_LABELS)
)

# Panel a: full horizon with regimes shaded ------------------------------------
pa <- ggplot() +
  geom_rect(data = bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.55) +
  # Mean loss is shown alongside CaR95 because the moving equilibrium is an
  # expected loss: it is the level the thin curve tracks, not the thick one.
  geom_line(data = car_curves,
            aes(x = year, y = mean_loss, colour = region),
            linewidth = 0.3, alpha = 0.75) +
  geom_line(data = car_curves,
            aes(x = year, y = car_95, colour = region), linewidth = 0.7) +
  geom_line(data = equil_moving,
            aes(x = year, y = moving_eq, colour = region),
            linetype = "dotted", linewidth = 0.45, show.legend = FALSE) +
  scale_fill_manual(values = setNames(PHASE_FILLS, PHASE_LABELS), name = NULL) +
  scale_colour_manual(values = region_cols, name = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 30, 80, 120, 160, 200)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(xlim = c(0, 200), ylim = c(0, NA)) +
  labs(x = "Horizon (years)",
       y = expression(CaR[95]~"(kg per tonne CO"[2]*"e)")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(0.02, "cm"))

# Single panel. An expanded 0-30 year inset was dropped: once the regimes are
# shaded, the transient is legible in the full-horizon panel, which is what
# main-text Fig 2b lacked.
ggsave(file.path(out_dir, "si_car_phases.pdf"), pa,
       width = 130, height = 100, units = "mm", bg = "white")
cat("Saved", file.path(out_dir, "si_car_phases.pdf"), "\n\n")

# Reportable numbers -----------------------------------------------------------
cat("Analytic equilibrium loss (gamma = 0), kg per tonne:\n")
equil_static %>%
  mutate(beta_bar_pct = round(100 * beta_bar, 3),
         r_pct = round(100 * r, 2),
         equilibrium = round(equilibrium, 1)) %>%
  select(region, r_pct, beta_bar_pct, equilibrium) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nCaR95 at the regime boundaries, kg per tonne:\n")
car_curves %>%
  filter(year %in% c(10, 30, 80, 100, 200)) %>%
  select(region, year, car_95) %>%
  mutate(car_95 = round(car_95, 1)) %>%
  pivot_wider(names_from = year, values_from = car_95,
              names_prefix = "yr") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nShare of the 200-year CaR reached by year 30:\n")
car_curves %>%
  group_by(region) %>%
  summarise(share_at_30 = car_95[year == 30] / car_95[year == 200],
            .groups = "drop") %>%
  mutate(share_at_30 = sprintf("%.0f%%", 100 * share_at_30)) %>%
  as.data.frame() %>% print(row.names = FALSE)
