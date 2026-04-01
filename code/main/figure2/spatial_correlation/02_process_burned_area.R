# Process MODIS burned area into cell-level annual burn fractions
#
# Prerequisites: Run 01_download_modis.R first
#
# Run from repo root:
#   Rscript code/main/figure2/spatial_correlation/02_process_burned_area.R
#   Rscript code/main/figure2/spatial_correlation/02_process_burned_area.R --cell_size 0.5
#
# Output:
#   - outputs/intermediate/cell_burn_fractions/cell_{size}/{region}.rds
#   - outputs/intermediate/cache/annual_burned_{region}_{year}.tif
#   - outputs/intermediate/cache/forest_mask_{region}.tif

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(parallel)

source("code/main/figure2/spatial_correlation/config.R")
source("code/main/figure2/spatial_correlation/helpers.R")

ensure_dirs()

# ── Load admin boundaries ──────────────────────────────────────────────────

gdf <- st_read(get_path("gpkg"), quiet = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Helper: Get or compute annual burned raster
# ══════════════════════════════════════════════════════════════════════════════

annual_cache_path <- function(region, year) {
  file.path(get_path("cache"), sprintf("annual_burned_%s_%d.tif", region, year))
}

get_annual_burned <- function(region, tiles, year, poly) {
  cache_path <- annual_cache_path(region, year)
  if (file.exists(cache_path)) {
    return(rast(cache_path))
  }

  mcd64a1_dir <- get_path("mcd64a1_dir")
  monthly_burned <- list()
  poly_vect <- NULL

  for (month in 1:12) {
    tile_rasters <- list()

    for (tile in tiles) {
      hdf <- find_hdf(tile, year, month, mcd64a1_dir)
      if (is.null(hdf)) next
      r <- tryCatch(load_burn_date(hdf), error = function(e) NULL)
      if (!is.null(r)) tile_rasters <- c(tile_rasters, list(r))
    }

    if (length(tile_rasters) == 0) next

    burn_date <- if (length(tile_rasters) == 1) {
      tile_rasters[[1]]
    } else {
      do.call(merge, tile_rasters)
    }

    if (is.null(poly_vect)) {
      poly_vect <- vect(st_transform(poly, crs(burn_date)))
    }

    burn_date <- crop(burn_date, poly_vect)
    burned <- burn_date > 0
    monthly_burned <- c(monthly_burned, list(burned))
  }

  if (length(monthly_burned) == 0) return(NULL)

  annual_burned <- if (length(monthly_burned) == 1) {
    monthly_burned[[1]]
  } else {
    app(rast(monthly_burned), max)
  }

  writeRaster(annual_burned, cache_path, overwrite = TRUE)
  annual_burned
}

# ══════════════════════════════════════════════════════════════════════════════
# Helper: Get or compute forest mask
# ══════════════════════════════════════════════════════════════════════════════

forest_mask_path <- function(region) {
  file.path(get_path("cache"), sprintf("forest_mask_%s.tif", region))
}

get_forest_mask <- function(region, tiles, poly) {
  cache_path <- forest_mask_path(region)
  if (file.exists(cache_path)) {
    return(rast(cache_path))
  }

  mcd12q1_dir <- get_path("mcd12q1_dir")
  tile_rasters <- list()

  for (tile in tiles) {
    hdf <- find_hdf_for_year(tile, PARAMS$lc_year, mcd12q1_dir, "MCD12Q1")
    if (is.null(hdf)) next
    r <- tryCatch(load_land_cover(hdf), error = function(e) NULL)
    if (!is.null(r)) tile_rasters <- c(tile_rasters, list(r))
  }

  if (length(tile_rasters) == 0) return(NULL)

  lc <- if (length(tile_rasters) == 1) {
    tile_rasters[[1]]
  } else {
    do.call(merge, tile_rasters)
  }

  poly_vect <- vect(st_transform(poly, crs(lc)))
  lc <- crop(lc, poly_vect)
  lc <- mask(lc, poly_vect)

  forest_mask <- lc %in% c(1:10, 15, 16)

  writeRaster(forest_mask, cache_path, overwrite = TRUE)
  forest_mask
}

# ══════════════════════════════════════════════════════════════════════════════
# Process each region
# ══════════════════════════════════════════════════════════════════════════════

tile_csv <- get_path("tile_lookup")
output_dir <- get_cell_path("cell_burn_fractions")
diag_dir <- get_cell_path("diagnostics")

for (reg in REGIONS) {
  cat(sprintf("\n═══ %s ═══\n", reg$name))

  tiles <- tryCatch(
    get_tiles_for_region(reg$country, reg$adm1, tile_csv),
    error = function(e) {
      cat(sprintf("  Error getting tiles: %s\n", e$message))
      NULL
    }
  )
  if (is.null(tiles)) next

  cat(sprintf("  Tiles: %s\n", paste(tiles, collapse = ", ")))

  poly <- gdf %>%
    filter(NAME_0 == reg$country, NAME_1 == reg$adm1) %>%
    st_make_valid()

  if (nrow(poly) == 0) {
    cat("  ADM1 polygon not found, skipping\n")
    next
  }
  poly <- st_transform(poly, 4326)

  # ── Create grid cells ────────────────────────────────────────────────────

  old_s2 <- sf_use_s2()
  sf_use_s2(FALSE)
  poly <- st_make_valid(poly)

  grid <- st_make_grid(poly, cellsize = PARAMS$cell_size_deg, what = "polygons")
  grid <- st_sf(cell_id = seq_along(grid), geometry = grid)
  grid <- st_intersection(grid, st_geometry(poly))
  sf_use_s2(old_s2)

  grid$cell_id <- seq_len(nrow(grid))
  cat(sprintf("  %d grid cells (%.2f deg)\n", nrow(grid), PARAMS$cell_size_deg))

  centroids <- st_centroid(grid)
  centroid_coords <- st_coordinates(centroids)

  grid_proj <- st_transform(grid, 6933)
  avg_cell_km2 <- mean(as.numeric(st_area(grid_proj))) / 1e6
  min_cell_km2 <- min(as.numeric(st_area(grid_proj))) / 1e6
  max_cell_km2 <- max(as.numeric(st_area(grid_proj))) / 1e6

  p_grid <- ggplot() +
    geom_sf(data = poly, fill = NA, color = "black") +
    geom_sf(data = grid, fill = NA) +
    geom_sf(data = centroids, color = "red", size = 0.5) +
    labs(
      title = sprintf("Grid cells for %s", reg$name),
      subtitle = sprintf("Avg: %.0f km2 | Min: %.0f km2 | Max: %.0f km2",
                         avg_cell_km2, min_cell_km2, max_cell_km2)
    ) +
    theme_minimal()

  ggsave(file.path(diag_dir, sprintf("grid_%s_%s.pdf", reg$name, PARAMS$cell_size_label)),
         p_grid, width = 6, height = 6)

  # ── Create forest mask ───────────────────────────────────────────────────

  cat("  Creating forest mask...\n")
  forest_mask <- get_forest_mask(reg$name, tiles, poly)
  if (is.null(forest_mask)) {
    cat("  Could not create forest mask\n")
  } else {
    cat(sprintf("  Forest mask saved to %s\n", forest_mask_path(reg$name)))
  }

  # ── Process annual burned area ───────────────────────────────────────────

  process_year <- function(year) {
    cat(sprintf("  %d: ", year))
    annual_burned <- get_annual_burned(reg$name, tiles, year, poly)
    if (is.null(annual_burned)) {
      cat("no data\n")
      return(NULL)
    }

    grid_vect <- vect(st_transform(grid, crs(annual_burned)))
    burn_frac <- terra::extract(annual_burned, grid_vect, fun = mean, na.rm = TRUE)
    cat(sprintf("%.4f mean\n", mean(burn_frac[, 2], na.rm = TRUE)))

    data.frame(
      cell_id = grid$cell_id,
      year = year,
      burn_fraction = burn_frac[, 2]
    )
  }

  results <- if (.Platform$OS.type == "unix" && PARAMS$n_workers > 1) {
    mclapply(PARAMS$years, process_year, mc.cores = PARAMS$n_workers)
  } else {
    lapply(PARAMS$years, process_year)
  }
  results <- Filter(function(x) is.data.frame(x), results)

  if (length(results) == 0) {
    cat("  No data for any year, skipping\n")
    next
  }

  cell_burns <- bind_rows(results)
  cell_burns$region <- reg$name

  attr(cell_burns, "centroids") <- centroid_coords
  attr(cell_burns, "n_cells") <- nrow(grid)

  out_path <- file.path(output_dir, sprintf("%s.rds", reg$name))
  saveRDS(cell_burns, out_path)
  cat(sprintf("  Saved to %s (%d rows)\n", out_path, nrow(cell_burns)))
}

cat("\n═══ Processing complete ═══\n")
