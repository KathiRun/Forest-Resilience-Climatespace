# ==============================================================================
# Create a point mesh of the functional biomes (for upload to Google Earth Engine)
#
# Description:
#   Converts the biome raster (band `grps90`, functional biome ID 1-12) into a
#   table of points, one point per raster cell with a biome value. Cells with
#   value 0 (no biome) and NA cells (masked, e.g. water) are dropped. Each point
#   gets its longitude, latitude, biome ID and a unique numeric ID. The table is
#   used in GEE to extract the climate variables of every point.
#
# Workflow:
#   1. Setup
#   2. Load the raster and select the biome band
#   3. Convert raster cells to points
#   4. Save
#
# Inputs:
#   biome_raster_export.tif
#       Raster with a band named `grps90` (functional biome ID)
#
# Outputs:
#   grid_points_biomes_mesh.csv
#       Columns: lon, lat, grps90, ID
#
# Note: lon/lat are the cell centres, in the coordinate system of the raster
# (geographic coordinates, e.g. WGS84, are assumed).
#
# Required packages (install beforehand):
#   terra, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(terra)   # raster handling
library(here)    # project-relative paths (run from the project root)

path_raster <- here("biome_raster_export.tif")
path_out    <- here("grid_points_biomes_mesh.csv")

band_name <- "grps90"   # raster band holding the biome ID


# ---- 2. Load the raster and select the biome band ----------------------------

env_clim <- rast(path_raster)

if (!band_name %in% names(env_clim)) {
  stop("Band '", band_name, "' not found. Available bands: ",
       paste(names(env_clim), collapse = ", "))
}
r_band <- env_clim[[band_name]]

# Value 0 means "no biome"
r_band[r_band == 0] <- NA


# ---- 3. Convert raster cells to points ---------------------------------------
# NA cells (masked / water / no biome) are dropped

pts <- as.points(r_band, na.rm = TRUE)

df <- data.frame(
  lon    = crds(pts)[, 1],
  lat    = crds(pts)[, 2],
  grps90 = values(pts)[[1]]
)

df$ID <- seq_len(nrow(df))

message("Created ", nrow(df), " points in ", length(unique(df$grps90)), " biomes")


# ---- 4. Save -----------------------------------------------------------------

write.csv(df, path_out, row.names = FALSE)
