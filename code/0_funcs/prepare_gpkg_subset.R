# Prepare the region-subset GeoPackage for the ncc replication package
#
# The full admin_regrowth_with_gpp.gpkg file (175 MB, 3,388 admin units
# worldwide) was created by joining GADM Level 1 administrative boundaries
# with mean Gross Primary Productivity (GPP) from the Bi et al. (2022)
# global GPP product (doi:10.1038/s41597-022-01309-2). The GPP values
# (columns gpp_mean, gpp_median, gpp_min, gpp_max) represent the average
# annual GPP (g C/m^2/yr) within each admin unit over 1992-2020.
#
# For the replication package, we include only the three regions used in
# the paper (California, Mato Grosso, Papua), reducing the file from
# 175 MB to ~800 KB.
#
# This script extracts the subset from the full file. It is provided for
# transparency but does not need to be run for replication, since the
# subset file is already included in data/.
#
# To run (from repo root, requires the full file at data/admin_regrowth_with_gpp.gpkg):
#   Rscript code/0_funcs/prepare_gpkg_subset.R

library(sf)

full_path   <- "data/admin_regrowth_with_gpp.gpkg"
subset_path <- "data/admin_regrowth_with_gpp.gpkg"

if (!file.exists(full_path)) {
  stop("Full GeoPackage not found at ", full_path,
       "\nThis script requires the full file (not included in the replication package).",
       "\nThe subset file at ", subset_path, " is already included and sufficient for replication.")
}

gpkg <- st_read(full_path, quiet = TRUE)
cat("Full file:", nrow(gpkg), "admin units\n")

# Extract the three regions used in the paper
subset <- gpkg[(gpkg$ISO == "USA" & gpkg$NAME_1 == "California") |
               (gpkg$ISO == "BRA" & gpkg$NAME_1 == "Mato Grosso") |
               (gpkg$ISO == "IDN" & gpkg$NAME_1 == "Papua"), ]

cat("Subset:", nrow(subset), "regions:",
    paste(subset$NAME_1, collapse = ", "), "\n")
cat("GPP means:", paste(round(subset$gpp_mean, 1), collapse = ", "), "\n")

# Remove existing file (may be a symlink)
if (file.exists(subset_path)) unlink(subset_path)

st_write(subset, subset_path, quiet = TRUE)
cat("Saved to", subset_path,
    sprintf("(%.0f KB)\n", file.size(subset_path) / 1024))
