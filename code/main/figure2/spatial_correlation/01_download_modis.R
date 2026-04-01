# Download MODIS burned area (MCD64A1) and land cover (MCD12Q1) tiles
#
# Prerequisites:
#   1. NASA Earthdata account (free: https://urs.earthdata.nasa.gov)
#   2. ~/.netrc with: machine urs.earthdata.nasa.gov login YOUR_USERNAME password YOUR_PASSWORD
#   3. ~/.urs_cookies file (touch ~/.urs_cookies)
#
# Run from repo root:
#   Rscript code/main/figure2/spatial_correlation/01_download_modis.R
#
# Output:
#   - outputs/intermediate/modis_mcd64a1/{tile}/*.hdf
#   - outputs/intermediate/modis_mcd12q1/{tile}/*.hdf
#   - outputs/intermediate/tile_lookup.csv

library(sf)
library(httr)
library(jsonlite)
library(parallel)
library(ggplot2)

source("code/main/figure2/spatial_correlation/config.R")
source("code/main/figure2/spatial_correlation/helpers.R")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Determine MODIS tiles for each region
# ══════════════════════════════════════════════════════════════════════════════

cat("═══ STEP 1: Determining MODIS tiles for each region ═══\n\n")

gdf <- st_read(get_path("gpkg"), quiet = TRUE)
if (is.na(st_crs(gdf))) {
  warning("ADM1 layer has no CRS; assuming EPSG:4326")
  st_crs(gdf) <- 4326
}
gdf_ll <- st_transform(gdf, 4326)

target <- data.frame(
  NAME_0 = vapply(REGIONS, function(r) r$country, character(1)),
  NAME_1 = vapply(REGIONS, function(r) r$adm1, character(1)),
  stringsAsFactors = FALSE
)
target$key <- paste(target$NAME_0, target$NAME_1, sep = " || ")

if (all(c("NAME_0", "NAME_1") %in% names(gdf_ll))) {
  gdf_key <- paste(gdf_ll$NAME_0, gdf_ll$NAME_1, sep = " || ")
  gdf_ll <- gdf_ll[gdf_key %in% target$key, , drop = FALSE]
}

hv <- expand.grid(h = 0:35, v = 0:17)
geoms <- mapply(make_tile_poly, hv$h, hv$v,
                MoreArgs = list(modis_config = MODIS), SIMPLIFY = FALSE)
tiles_sinu <- st_sf(
  h = hv$h,
  v = hv$v,
  geometry = st_sfc(geoms, crs = MODIS$sinu_proj)
)

old_s2 <- sf_use_s2()
sf_use_s2(FALSE)
on.exit(sf_use_s2(old_s2), add = TRUE)

if ("st_make_valid" %in% getNamespaceExports("sf")) {
  gdf_ll <- st_make_valid(gdf_ll)
}

gdf_sinu <- st_transform(gdf_ll, MODIS$sinu_proj)
hits <- st_intersects(gdf_sinu, tiles_sinu)

fmt_h <- function(x) sprintf("%02d", x)
fmt_v <- function(x) sprintf("%02d", x)

adm1_id <- paste(gdf_ll$NAME_0, gdf_ll$NAME_1, sep = " || ")

rows <- vector("list", length(hits))
for (i in seq_along(hits)) {
  idx <- hits[[i]]
  if (length(idx) == 0) {
    rows[[i]] <- data.frame(adm1_key = adm1_id[i], tile = character(0),
                            stringsAsFactors = FALSE)
  } else {
    tiles_i <- paste0("h", fmt_h(tiles_sinu$h[idx]), "v", fmt_v(tiles_sinu$v[idx]))
    rows[[i]] <- data.frame(adm1_key = adm1_id[i], tile = tiles_i,
                            stringsAsFactors = FALSE)
  }
}

tile_df <- do.call(rbind, rows)

tile_csv <- get_path("tile_lookup")
dir.create(dirname(tile_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(tile_df, tile_csv, row.names = FALSE)
cat(sprintf("Wrote %d rows to %s\n", nrow(tile_df), tile_csv))

TILES <- sort(unique(tile_df$tile))
cat(sprintf("Unique tiles to download: %s\n\n", paste(TILES, collapse = ", ")))

# Diagnostic plots
diag_dir <- get_path("diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(gdf_ll))) {
  region_key <- paste(gdf_ll$NAME_0[i], gdf_ll$NAME_1[i], sep = " || ")
  tiles_i <- tile_df[tile_df$adm1_key == region_key, , drop = FALSE]
  if (nrow(tiles_i) == 0) next

  tile_ids <- tiles_i$tile
  tile_idx <- which(paste0("h", fmt_h(tiles_sinu$h), "v", fmt_v(tiles_sinu$v)) %in% tile_ids)
  if (length(tile_idx) == 0) next

  poly <- gdf_ll[i, ]
  poly_bbox <- st_bbox(poly)

  pad_x <- (poly_bbox$xmax - poly_bbox$xmin) * 0.5
  pad_y <- (poly_bbox$ymax - poly_bbox$ymin) * 0.5
  xlim <- c(poly_bbox$xmin - pad_x, poly_bbox$xmax + pad_x)
  ylim <- c(poly_bbox$ymin - pad_y, poly_bbox$ymax + pad_y)

  tiles_plot <- st_transform(tiles_sinu[tile_idx, ], 4326)
  poly_plot <- st_transform(poly, 4326)

  cent <- st_centroid(tiles_sinu[tile_idx, ])
  cent_plot <- st_transform(cent, 4326)
  coords <- st_coordinates(cent_plot)
  labels <- paste0("h", fmt_h(tiles_sinu$h[tile_idx]), "v", fmt_v(tiles_sinu$v[tile_idx]))
  label_df <- data.frame(x = coords[, 1], y = coords[, 2], label = labels)

  p <- ggplot() +
    geom_sf(data = tiles_plot, fill = NA, color = "grey60", linewidth = 0.6) +
    geom_sf(data = poly_plot, fill = "steelblue", color = "steelblue", alpha = 0.3) +
    geom_text(data = label_df, aes(x = x, y = y, label = label), size = 3) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = paste("MODIS Tiles:", region_key)) +
    theme_minimal(base_size = 12)

  pdf_path <- file.path(diag_dir, sprintf("tiles_%s_%s.pdf",
                        gsub("\\s+", "_", gdf_ll$NAME_0[i]),
                        gsub("\\s+", "_", gdf_ll$NAME_1[i])))
  ggsave(pdf_path, plot = p, width = 10, height = 7)
  cat(sprintf("Wrote %s\n", pdf_path))
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Download MCD64A1 (Burned Area)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n═══ STEP 2: Downloading MCD64A1 (Burned Area) ═══\n\n")

mcd64a1_dir <- get_path("mcd64a1_dir")

TILE_BBOX <- setNames(
  lapply(TILES, tile_bbox_lonlat, modis_config = MODIS),
  TILES
)

query_cmr_mcd64a1 <- function(tile, year) {
  temporal <- sprintf("%d-01-01T00:00:00Z,%d-12-31T23:59:59Z", year, year)
  bbox <- TILE_BBOX[[tile]]

  all_urls <- character(0)
  all_names <- character(0)
  page <- 1

  repeat {
    resp <- httr::GET(MODIS$cmr_url, query = list(
      short_name = "MCD64A1",
      version = "061",
      provider = "LPCLOUD",
      bounding_box = bbox,
      temporal = temporal,
      page_size = 500,
      page_num = page
    ))

    if (httr::status_code(resp) != 200) {
      warning(sprintf("CMR query failed for %s/%d (HTTP %d)",
                      tile, year, httr::status_code(resp)))
      break
    }

    body <- httr::content(resp, as = "parsed", simplifyVector = TRUE)
    entries <- body$feed$entry

    if (is.null(entries) || length(entries) == 0 ||
        (is.data.frame(entries) && nrow(entries) == 0)) break

    for (i in seq_len(nrow(entries))) {
      title <- entries$title[i]
      if (!grepl(tile, title, fixed = TRUE)) next

      links <- entries$links[[i]]
      hdf_links <- links$href[grepl("\\.hdf$", links$href) &
                               grepl("^https://", links$href)]
      if (length(hdf_links) > 0) {
        all_urls <- c(all_urls, hdf_links[1])
        all_names <- c(all_names, title)
      }
    }

    if (nrow(entries) < 500) break
    page <- page + 1
  }

  data.frame(name = all_names, url = all_urls, stringsAsFactors = FALSE)
}

for (tile in TILES) {
  dir.create(file.path(mcd64a1_dir, tile), recursive = TRUE, showWarnings = FALSE)
}

n_download <- n_skip <- n_fail <- 0

for (tile in TILES) {
  for (year in PARAMS$years) {
    cat(sprintf("[MCD64A1 %s / %d] Querying CMR... ", tile, year))

    granules <- tryCatch(query_cmr_mcd64a1(tile, year), error = function(e) {
      cat(sprintf("ERROR: %s\n", e$message))
      NULL
    })

    if (is.null(granules) || nrow(granules) == 0) {
      cat("no granules found\n")
      next
    }

    cat(sprintf("%d granules\n", nrow(granules)))

    tasks <- lapply(seq_len(nrow(granules)), function(j) {
      fname <- paste0(granules$name[j], ".hdf")
      dest <- file.path(mcd64a1_dir, tile, fname)
      list(url = granules$url[j], dest = dest)
    })

    counts <- download_many(tasks, PARAMS$n_workers, AUTH)
    n_download <- n_download + counts$ok
    n_skip <- n_skip + counts$skip
    n_fail <- n_fail + counts$fail
  }
}

cat(sprintf("\nMCD64A1 done. Downloaded: %d, Skipped: %d, Failed: %d\n",
            n_download, n_skip, n_fail))

# Retry corrupted files
corrupt <- list.files(mcd64a1_dir, pattern = "\\.hdf$", recursive = TRUE, full.names = TRUE)
corrupt <- corrupt[file.size(corrupt) == 0]

if (length(corrupt) > 0) {
  cat(sprintf("\nFound %d corrupt file(s). Retrying...\n", length(corrupt)))
  file.remove(corrupt)

  for (attempt in seq_len(PARAMS$max_retries)) {
    cat(sprintf("--- Retry %d/%d ---\n", attempt, PARAMS$max_retries))
    n_fixed <- 0

    for (tile in TILES) {
      for (year in PARAMS$years) {
        granules <- tryCatch(query_cmr_mcd64a1(tile, year), error = function(e) NULL)
        if (is.null(granules) || nrow(granules) == 0) next

        tasks <- lapply(seq_len(nrow(granules)), function(j) {
          fname <- paste0(granules$name[j], ".hdf")
          dest <- file.path(mcd64a1_dir, tile, fname)
          if (file.exists(dest) && file.size(dest) > 0) NULL
          else list(url = granules$url[j], dest = dest)
        })
        tasks <- Filter(Negate(is.null), tasks)
        if (length(tasks) == 0) next

        counts <- download_many(tasks, PARAMS$n_workers, AUTH)
        n_fixed <- n_fixed + counts$ok
      }
    }

    remaining <- list.files(mcd64a1_dir, pattern = "\\.hdf$", recursive = TRUE, full.names = TRUE)
    remaining <- remaining[file.size(remaining) == 0]
    if (length(remaining) == 0) {
      cat("All MCD64A1 files OK.\n")
      break
    } else {
      file.remove(remaining)
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Download MCD12Q1 (Land Cover)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n═══ STEP 3: Downloading MCD12Q1 (Land Cover) ═══\n\n")

mcd12q1_dir <- get_path("mcd12q1_dir")

query_cmr_mcd12q1 <- function(tile, year) {
  temporal <- sprintf("%d-01-01T00:00:00Z,%d-12-31T23:59:59Z", year, year)
  bbox <- TILE_BBOX[[tile]]

  all_urls <- character(0)
  all_names <- character(0)
  page <- 1

  repeat {
    resp <- httr::GET(MODIS$cmr_url, query = list(
      short_name = "MCD12Q1",
      version = "061",
      provider = "LPCLOUD",
      bounding_box = bbox,
      temporal = temporal,
      page_size = 500,
      page_num = page
    ))

    if (httr::status_code(resp) != 200) {
      warning(sprintf("CMR query failed for %s/%d (HTTP %d)",
                      tile, year, httr::status_code(resp)))
      break
    }

    body <- httr::content(resp, as = "parsed", simplifyVector = TRUE)
    entries <- body$feed$entry

    if (is.null(entries) || length(entries) == 0 ||
        (is.data.frame(entries) && nrow(entries) == 0)) break

    for (i in seq_len(nrow(entries))) {
      title <- entries$title[i]
      if (!grepl(tile, title, fixed = TRUE)) next

      links <- entries$links[[i]]
      hdf_links <- links$href[grepl("\\.hdf$", links$href) &
                               grepl("^https://", links$href)]
      if (length(hdf_links) > 0) {
        all_urls <- c(all_urls, hdf_links[1])
        all_names <- c(all_names, title)
      }
    }

    if (nrow(entries) < 500) break
    page <- page + 1
  }

  data.frame(name = all_names, url = all_urls, stringsAsFactors = FALSE)
}

for (tile in TILES) {
  dir.create(file.path(mcd12q1_dir, tile), recursive = TRUE, showWarnings = FALSE)
}

n_download <- n_skip <- n_fail <- 0

for (tile in TILES) {
  cat(sprintf("[MCD12Q1 %s / %d] Querying CMR... ", tile, PARAMS$lc_year))

  granules <- tryCatch(query_cmr_mcd12q1(tile, PARAMS$lc_year), error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    NULL
  })

  if (is.null(granules) || nrow(granules) == 0) {
    cat("no granules found\n")
    next
  }

  cat(sprintf("%d granules\n", nrow(granules)))

  tasks <- lapply(seq_len(nrow(granules)), function(j) {
    fname <- paste0(granules$name[j], ".hdf")
    dest <- file.path(mcd12q1_dir, tile, fname)
    list(url = granules$url[j], dest = dest)
  })

  counts <- download_many(tasks, PARAMS$n_workers, AUTH)
  n_download <- n_download + counts$ok
  n_skip <- n_skip + counts$skip
  n_fail <- n_fail + counts$fail
}

cat(sprintf("\nMCD12Q1 done. Downloaded: %d, Skipped: %d, Failed: %d\n",
            n_download, n_skip, n_fail))

cat("\n═══ Download complete ═══\n")
