# Shared configuration for MODIS spatial correlation pipeline
# All paths are relative to project root (run from repo root)

# ── Region Definitions ─────────────────────────────────────────────────────
REGIONS <- list(
  list(
    name = "California",
    country = "United States",
    adm1 = "California",
    geo_id = "USA.5_1"
  ),
  list(
    name = "Mato_Grosso",
    country = "Brazil",
    adm1 = "Mato Grosso",
    geo_id = "BRA.14_1"
  ),
  list(
    name = "Papua",
    country = "Indonesia",
    adm1 = "Papua",
    geo_id = "IDN.24_1"
  )
)

# ── File Paths ─────────────────────────────────────────────────────────────
# All paths relative to project root, routed through 
PATHS <- list(
  # Input data
  gpkg = "data/admin_regrowth_with_gpp.gpkg",

  # MODIS data directories (downloaded artifacts)
  mcd64a1_dir = "outputs/intermediate/modis_mcd64a1",
  mcd12q1_dir = "outputs/intermediate/modis_mcd12q1",

  # Output directories
  tile_lookup = "outputs/intermediate/tile_lookup.csv",
  cell_burn_fractions = "outputs/intermediate/cell_burn_fractions",
  correlation_results = "outputs/intermediate/correlation_results",
  cache = "outputs/intermediate/cache",
  diagnostics = "outputs/intermediate/diagnostics"
)

# ── MODIS Grid Constants ───────────────────────────────────────────────────
MODIS <- list(
  sinu_proj = "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext",
  tile_size = 1111950,  # meters
  x_min = -20015109.354,
  y_max = 10007554.677,
  cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.json"
)

# ── Simulation Parameters ──────────────────────────────────────────────────
PARAMS <- list(
  years = 2002:2023,
  cell_size_deg = 1,  # ~100 km at equator
  lc_year = 2020,
  min_fire_years = 5,
  n_workers = 6,
  max_retries = 3
)

# ── Command-Line Override for Cell Size ───────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
cell_idx <- which(args == "--cell_size")
if (length(cell_idx) > 0 && cell_idx < length(args)) {
  PARAMS$cell_size_deg <- as.numeric(args[cell_idx + 1])
  cat(sprintf("Cell size overridden to: %.2f degrees\n", PARAMS$cell_size_deg))
}

PARAMS$cell_size_label <- sprintf("%.2fdeg", PARAMS$cell_size_deg)

# ── Earthdata Authentication ───────────────────────────────────────────────
AUTH <- list(
  netrc = path.expand("~/.netrc"),
  cookies = path.expand("~/.urs_cookies")
)

# ── Helper to get full path ────────────────────────────────────────────────
get_path <- function(key) {
  path <- PATHS[[key]]
  if (is.null(path)) stop(sprintf("Unknown path key: %s", key))
  path
}

# ── Helper to get cell-size-specific path ──────────────────────────────────
get_cell_path <- function(key) {
  base <- get_path(key)
  file.path(base, sprintf("cell_%s", PARAMS$cell_size_label))
}

# ── Ensure output directories exist ────────────────────────────────────────
ensure_dirs <- function() {
  dir.create(get_path("cache"), recursive = TRUE, showWarnings = FALSE)
  cell_dirs <- c("cell_burn_fractions", "correlation_results", "diagnostics")
  for (d in cell_dirs) {
    dir.create(get_cell_path(d), recursive = TRUE, showWarnings = FALSE)
  }
}
