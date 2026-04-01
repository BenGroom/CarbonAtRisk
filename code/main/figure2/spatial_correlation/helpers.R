# Shared helper functions for MODIS spatial correlation pipeline

# ── Tile Lookup ────────────────────────────────────────────────────────────

#' Get MODIS tiles that intersect a region
#' @param country Country name (NAME_0)
#' @param adm1 ADM1 name (NAME_1)
#' @param tile_csv Path to tile lookup CSV
#' @return Character vector of tile IDs (e.g., "h08v04")
get_tiles_for_region <- function(country, adm1, tile_csv) {
  if (!file.exists(tile_csv)) {
    stop(sprintf("Tile CSV not found at %s. Run 01_download_modis.R first.", tile_csv))
  }

  df <- read.csv(tile_csv, stringsAsFactors = FALSE)
  if (!"tile" %in% names(df)) stop("Tile CSV missing 'tile' column.")

  if ("adm1_key" %in% names(df)) {
    key <- paste(country, adm1, sep = " || ")
    tiles <- df$tile[df$adm1_key == key]
  } else if (all(c("NAME_0", "NAME_1") %in% names(df))) {
    tiles <- df$tile[df$NAME_0 == country & df$NAME_1 == adm1]
  } else {
    stop("Tile CSV missing adm1_key or NAME_0/NAME_1 columns.")
  }

  tiles <- sort(unique(tiles))
  if (length(tiles) == 0) {
    stop(sprintf("No tiles found for %s / %s in %s", country, adm1, tile_csv))
  }
  tiles
}

# ── Tile Geometry ──────────────────────────────────────────────────────────

#' Get bounding box for a MODIS tile in lon/lat
tile_bbox_lonlat <- function(tile_id, modis_config) {
  h <- as.integer(sub("^h(\\d{2})v\\d{2}$", "\\1", tile_id))
  v <- as.integer(sub("^h\\d{2}v(\\d{2})$", "\\1", tile_id))
  if (is.na(h) || is.na(v)) stop(sprintf("Bad tile id: %s", tile_id))

  x0 <- modis_config$x_min + (h * modis_config$tile_size)
  x1 <- x0 + modis_config$tile_size
  y1 <- modis_config$y_max - (v * modis_config$tile_size)
  y0 <- y1 - modis_config$tile_size

  poly <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)
  ))), crs = modis_config$sinu_proj)

  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  poly_ll <- sf::st_transform(poly, 4326)
  bb <- sf::st_bbox(poly_ll)
  sprintf("%f,%f,%f,%f", bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])
}

#' Create MODIS tile polygon in sinusoidal projection
make_tile_poly <- function(h, v, modis_config) {
  x0 <- modis_config$x_min + (h * modis_config$tile_size)
  x1 <- x0 + modis_config$tile_size
  y1 <- modis_config$y_max - (v * modis_config$tile_size)
  y0 <- y1 - modis_config$tile_size
  sf::st_polygon(list(rbind(
    c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)
  )))
}

# ── Download Functions ─────────────────────────────────────────────────────

#' Download a file from NASA Earthdata
download_earthdata <- function(url, dest, auth) {
  httr::GET(
    url,
    httr::config(
      netrc = TRUE,
      netrc_file = auth$netrc,
      cookiefile = auth$cookies,
      cookiejar = auth$cookies,
      followlocation = TRUE
    ),
    httr::write_disk(dest, overwrite = FALSE)
  )
}

#' Download a single file with status tracking
download_one <- function(task, auth) {
  dest <- task$dest
  if (file.exists(dest) && file.size(dest) > 0) {
    return(list(status = "skip", dest = dest))
  }

  dl <- tryCatch(download_earthdata(task$url, dest, auth), error = function(e) NULL)

  if (!is.null(dl) && httr::status_code(dl) == 200 &&
      file.exists(dest) && file.size(dest) > 0) {
    list(status = "ok", dest = dest)
  } else {
    if (file.exists(dest)) file.remove(dest)
    list(status = "fail", dest = dest)
  }
}

#' Download multiple files in parallel
download_many <- function(tasks, workers, auth) {
  if (length(tasks) == 0) return(list(ok = 0, skip = 0, fail = 0))

  fn <- function(task) download_one(task, auth)

  if (.Platform$OS.type == "unix") {
    res <- parallel::mclapply(tasks, fn, mc.cores = workers)
  } else {
    res <- lapply(tasks, fn)
  }

  tab <- table(vapply(res, function(x) x$status, character(1)))
  list(
    ok = if ("ok" %in% names(tab)) unname(tab["ok"]) else 0,
    skip = if ("skip" %in% names(tab)) unname(tab["skip"]) else 0,
    fail = if ("fail" %in% names(tab)) unname(tab["fail"]) else 0
  )
}

# ── HDF Loading Functions ──────────────────────────────────────────────────

#' Find HDF file for a tile/year/month
find_hdf <- function(tile, year, month, data_dir, max_doy_diff = 20) {
  tile_dir <- file.path(data_dir, tile)
  if (!dir.exists(tile_dir)) return(NULL)

  files <- list.files(tile_dir, pattern = "\\.hdf$", full.names = TRUE)
  year_files <- grep(sprintf("MCD64A1\\.A%d", year), files, value = TRUE)
  if (length(year_files) == 0) return(NULL)

  month_start_doy <- as.integer(format(as.Date(sprintf("%d-%02d-01", year, month)), "%j"))
  doys <- as.integer(sub(".*\\.A\\d{4}(\\d{3})\\..*", "\\1", basename(year_files)))
  idx <- which.min(abs(doys - month_start_doy))

  if (!is.na(max_doy_diff) && abs(doys[idx] - month_start_doy) > max_doy_diff) return(NULL)
  year_files[idx]
}

#' Find HDF file for a tile/year (any DOY)
find_hdf_for_year <- function(tile, year, data_dir, product_prefix) {
  tile_dir <- file.path(data_dir, tile)
  if (!dir.exists(tile_dir)) return(NULL)

  files <- list.files(tile_dir, pattern = "\\.hdf$", full.names = TRUE)
  year_files <- grep(sprintf("%s\\.A%d", product_prefix, year), files, value = TRUE)
  if (length(year_files) == 0) return(NULL)
  year_files[1]
}

#' Load burn date layer from MCD64A1 HDF
load_burn_date <- function(hdf_path) {
  sds_name <- sprintf(
    'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Monthly_500m_DB_BA:"Burn Date"',
    hdf_path
  )
  terra::rast(sds_name)
}

#' Load land cover type 1 (IGBP) from MCD12Q1 HDF
load_land_cover <- function(hdf_path) {
  sds_name <- sprintf(
    'HDF4_EOS:EOS_GRID:"%s":MCD12Q1:LC_Type1',
    hdf_path
  )
  terra::rast(sds_name)
}
