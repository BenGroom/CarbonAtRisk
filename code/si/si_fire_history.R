# SI: Fire History Bar Charts for Three Regions
# Produces a faceted bar chart of annual burned area (ha) for
# California, Mato Grosso, and Papua (2002-2023).
#
# Run from repo root:
#   Rscript code/si/si_fire_history.R
#
# Output:
#   - outputs/si/si_fire_history.pdf

library(ggplot2)
library(dplyr)
library(sf)
library(httr)
library(jsonlite)

source("code/0_funcs/fire_funcs.R")

# Configuration ---------------------------------------------------------------
GPKG_PATH <- "data/admin_regrowth_with_gpp.gpkg"
EFFIS_CACHE <- "data/effis_cache"

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

regions <- list(
  list(country = "United States", subcountry = "California"),
  list(country = "Brazil",        subcountry = "Mato Grosso"),
  list(country = "Indonesia",     subcountry = "Papua")
)

# Load data --------------------------------------------------------------------
gdf_gadm <- st_read(GPKG_PATH, quiet = TRUE)

all_fires <- data.frame()

for (reg in regions) {
  row <- gdf_gadm %>%
    filter(NAME_0 == reg$country, NAME_1 == reg$subcountry) %>%
    slice(1)

  geo_id <- paste0(row$country_code, ".", row$ID_1, "_1")
  fires_df <- fetch_fire_data(geo_id, cache_dir = EFFIS_CACHE)

  fires_df$region <- reg$subcountry
  all_fires <- bind_rows(all_fires, fires_df)
}

# Plot -------------------------------------------------------------------------
p <- all_fires %>%
  arrange(year) %>%
  ggplot(aes(x = factor(year), y = ba_area_ha)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  facet_wrap(~ region, scales = "free_y", ncol = 1) +
  labs(
    x = "Year",
    y = "Burned Area (ha)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave(file.path(out_dir, "si_fire_history.pdf"), p, width = 10, height = 10)

cat("Saved to outputs/si/si_fire_history.pdf\n")
