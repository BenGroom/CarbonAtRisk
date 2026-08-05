# SI: Fire history for the three focal regions
# Faceted bar chart of the annual burn fraction for California, Mato Grosso and
# Papua (2002-2023), from EFFIS.
#
# Plots the burn *fraction* rather than burned area in hectares. Hectares do not
# normalise for region size, so the three panels could not be compared with each
# other, nor with the conversion record in si_deforestation_history.pdf. The
# fraction is also the quantity the simulation actually resamples.
#
# Run from repo root:
#   Rscript code/si/si_fire_history.R
#
# Output:
#   - outputs/si/si_fire_history.pdf

library(tidyverse)

source("code/0_funcs/fire_funcs.R")
source("code/si/si_deforestation_common.R")

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load data --------------------------------------------------------------------
# Same crosswalk, cache and parameters as figure2.R and the conversion analysis,
# so all three hazard figures describe the same underlying quantities.
fire_hist <- pmap_dfr(REGIONS, function(region, effis_key, iso, adm1, gid_1) {
  fires  <- read.csv(file.path(EFFIS_CACHE, sprintf("effis_fire_%s.csv", effis_key)))
  forest <- readRDS(file.path(EFFIS_CACHE, sprintf("effis_forest_%s.rds", effis_key)))
  bf <- calculate_burn_fractions(fires, as.numeric(forest$lc1),
                                 project_area = PROJECT_AREA,
                                 rescale_firesize = RESCALE_FIRESIZE)
  tibble(region = region, year = fires$year, burn_fraction = as.numeric(bf))
}) %>%
  mutate(region = factor(region, levels = REGIONS$region))

mean_bf <- fire_hist %>%
  group_by(region) %>%
  summarise(beta_bar = mean(burn_fraction), .groups = "drop")

# Plot -------------------------------------------------------------------------
p <- ggplot(fire_hist, aes(x = year, y = 100 * burn_fraction)) +
  geom_col(fill = "#E67E22", width = 0.75) +
  geom_hline(data = mean_bf, aes(yintercept = 100 * beta_bar),
             linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  facet_wrap(~region, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(2002, 2022, 4)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Annual burn fraction (% of forest area per year)") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9, hjust = 0))

ggsave(file.path(out_dir, "si_fire_history.pdf"), p,
       width = 110, height = 130, units = "mm", bg = "white")
cat("Saved", file.path(out_dir, "si_fire_history.pdf"), "\n\n")

cat("Annual burn fraction, %/yr:\n")
fire_hist %>%
  group_by(region) %>%
  summarise(mean = 100 * mean(burn_fraction),
            max  = 100 * max(burn_fraction),
            peak_year = year[which.max(burn_fraction)],
            n_years = n(), .groups = "drop") %>%
  mutate(across(c(mean, max), ~round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
