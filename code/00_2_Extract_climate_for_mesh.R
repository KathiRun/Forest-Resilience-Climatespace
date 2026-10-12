# ==============================================================================
# Extract climate values at the biome point mesh, one file per biome
#
# Description:
#   Extracts the values of all climate layers of the global raster at every point
#   of the biome point mesh and writes one csv per functional biome (grps90).
#   The biome band is excluded from the raster, as the biome ID is already part
#   of the point table. The output files are the input of the PCA step.
#
# Workflow:
#   1. Setup
#   2. Load the raster and the points
#   3. Extract climate values per biome
#
# Inputs:
#   gee_exports/global_raster_clim.tif   multi-band raster of climate layers
#   grid_points_biomes_mesh.csv          point mesh
#       Columns: lon, lat, grps90, ID
#
# Outputs:
#   by_grps90/points_clim_grps90_<g>.csv
#       Point table (lon, lat, grps90, ID) plus one column per climate layer
#
# Required packages (install beforehand):
#   terra, tidyverse, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(terra)       # raster extraction
library(tidyverse)   # readr::write_csv
library(here)        # project-relative paths (run from the project root)

path_raster <- here("gee_exports", "global_raster_clim.tif")
path_points <- here("grid_points_biomes_mesh.csv")
out_dir     <- here("by_grps90")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(memfrac = 0.8)   # let terra use more of the available RAM (default is conservative)


# ---- 2. Load the raster and the points ---------------------------------------

env_clim <- rast(path_raster)
env_clim <- env_clim[[names(env_clim) != "grps90"]]   # drop the biome band, if present

points <- read.csv(path_points)

stopifnot("grps90" %in% names(points))
grps_values <- sort(unique(points$grps90))
message("Found ", length(grps_values), " unique grps90 values: ",
        paste(grps_values, collapse = ", "))


# ---- 3. Extract climate values per biome -------------------------------------
# Processed biome by biome to limit memory use

for (g in grps_values) {

  message("Processing grps90 = ", g)

  pts_sub    <- points[points$grps90 == g, ]
  pts_sub_sp <- vect(pts_sub, geom = c("lon", "lat"), crs = "EPSG:4326")

  # Raster values at the points (same row order as pts_sub)
  extracted <- terra::extract(env_clim, pts_sub_sp, ID = FALSE)
  names(extracted) <- names(env_clim)

  results <- cbind(pts_sub, extracted)

  write_csv(results, file.path(out_dir, paste0("points_clim_grps90_", g, ".csv")))

  rm(pts_sub, pts_sub_sp, extracted, results)
  gc()
}

message("Done. All biome files written to ", out_dir)
