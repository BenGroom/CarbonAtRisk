# SI figure: the conversion record for the three focal regions
#
# Companion to si_fire_history.R. Describes the second hazard's raw record before it
# is used, mirroring how the fire history precedes the fire results.
#
# Panel a: annual conversion rate by region, with full-record and period means.
# Panel b: mean fire hazard against mean conversion hazard, per region.
#
# Run from repo root:  Rscript code/si/si_deforestation_history.R

library(tidyverse)
library(patchwork)
library(scales)

source("code/0_funcs/fire_funcs.R")
source("code/si/si_deforestation_common.R")

FIG_WIDTH_MM <- 183
PERIODS <- tribble(
  ~period,     ~from, ~to,
  "2001-08",    2001, 2008,
  "2009-16",    2009, 2016,
  "2017-25",    2017, 2025
)

out_dir <- "outputs/si"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Data -------------------------------------------------------------------------
defor <- load_defor_panel()
burn  <- load_burn_fractions()

region_levels <- REGIONS$region
defor <- defor %>% mutate(region = factor(region, levels = region_levels))

full_mean <- defor %>%
  group_by(region) %>%
  summarise(delta_bar = mean(delta), .groups = "drop")

period_mean <- defor %>%
  cross_join(PERIODS) %>%
  filter(year >= from, year <= to) %>%
  group_by(region, period, from, to) %>%
  summarise(delta_bar = mean(delta), .groups = "drop")

# Panel a: the annual record ----------------------------------------------------
pa <- ggplot(defor, aes(x = year, y = 100 * delta)) +
  geom_col(fill = "#8C6D46", width = 0.75) +
  geom_hline(data = full_mean, aes(yintercept = 100 * delta_bar),
             linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_segment(data = period_mean,
               aes(x = from - 0.4, xend = to + 0.4,
                   y = 100 * delta_bar, yend = 100 * delta_bar),
               colour = "#C0392B", linewidth = 0.8) +
  facet_wrap(~region, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(2002, 2024, 4)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Annual conversion rate (% of forest area per year)") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9, hjust = 0))

# Panel b: the two hazards side by side ------------------------------------------
haz <- full_mean %>%
  left_join(burn %>%
              group_by(region) %>%
              summarise(beta_bar = mean(burn_fraction), .groups = "drop"),
            by = "region") %>%
  pivot_longer(c(beta_bar, delta_bar), names_to = "hazard", values_to = "rate") %>%
  mutate(hazard = factor(hazard, levels = c("beta_bar", "delta_bar"),
                         labels = c("Fire (modelled)", "Conversion (omitted)")),
         region = factor(region, levels = region_levels))

# A dot plot rather than bars: the axis is logarithmic, so bars would have no
# meaningful baseline and their lengths would be an artefact of the axis limits.
haz_seg <- haz %>%
  select(region, hazard, rate) %>%
  pivot_wider(names_from = hazard, values_from = rate) %>%
  rename(fire = `Fire (modelled)`, conv = `Conversion (omitted)`)

pb <- ggplot() +
  geom_segment(data = haz_seg,
               aes(x = region, xend = region, y = 100 * fire, yend = 100 * conv),
               colour = "grey65", linewidth = 0.5) +
  geom_point(data = haz, aes(x = region, y = 100 * rate, colour = hazard),
             size = 2.6) +
  geom_text(data = haz, aes(x = region, y = 100 * rate,
                            label = sprintf("%.3f", 100 * rate)),
            hjust = -0.35, size = 2.1) +
  scale_colour_manual(values = c("Fire (modelled)" = "#E67E22",
                                 "Conversion (omitted)" = "#8C6D46"), name = NULL) +
  scale_y_log10(limits = c(0.03, 6),
                breaks = c(0.05, 0.1, 0.3, 1, 3),
                labels = c("0.05", "0.1", "0.3", "1", "3")) +
  labs(x = NULL, y = "Mean annual hazard (%/yr, log scale)") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))

tag_theme <- theme(plot.tag = element_text(face = "bold", size = 10))
fig <- ((pa + tag_theme) | (pb + tag_theme)) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(out_dir, "si_deforestation_history.pdf"), fig,
       width = FIG_WIDTH_MM, height = FIG_WIDTH_MM * 0.62,
       units = "mm", bg = "white")
cat("Saved", file.path(out_dir, "si_deforestation_history.pdf"), "\n\n")

# Reportable numbers -------------------------------------------------------------
cat("Mean annual conversion rate, %/yr:\n")
period_mean %>%
  select(region, period, delta_bar) %>%
  mutate(delta_bar = round(100 * delta_bar, 3)) %>%
  pivot_wider(names_from = period, values_from = delta_bar) %>%
  left_join(full_mean %>% mutate(full = round(100 * delta_bar, 3)) %>%
              select(region, full), by = "region") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nMean annual hazard comparison, %/yr:\n")
haz %>%
  mutate(rate = round(100 * rate, 4)) %>%
  pivot_wider(names_from = hazard, values_from = rate) %>%
  mutate(ratio_fire_to_conversion =
           round(`Fire (modelled)` / `Conversion (omitted)`, 2),
         dominant = if_else(`Fire (modelled)` > `Conversion (omitted)`,
                            "fire", "CONVERSION")) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nPeak conversion year per region:\n")
defor %>%
  group_by(region) %>%
  slice_max(delta, n = 1) %>%
  mutate(delta = round(100 * delta, 3)) %>%
  select(region, year, delta) %>%
  as.data.frame() %>% print(row.names = FALSE)
