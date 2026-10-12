# ==============================================================================
# Assign unique IDs to the sample points
#
# Description:
#   Gives every sample point a unique ID (UUID), extracts longitude and latitude
#   from the point geometry, and writes
#     - the full table of sample points, and
#     - a small table (id, lat, lon) used to extract NDVI in Google Earth Engine.
#   As a check, the distribution of PC1 per biome is plotted.
#
#   IMPORTANT: UUIDs are random. The IDs created here are used in all later
#   steps (e.g. the NDVI extraction in GEE), so run this script only once.
#   It stops if the output already exists (see `overwrite_ids`).
#
# Workflow:
#   1. Setup
#   2. Load data and prepare columns
#   3. Assign IDs and save
#   4. Check the PC1 distribution per biome
#
# Inputs:
#   gee_exports/sample_points_with_composite.csv
#       Columns used: geometry (WKT point), first (forest mask), grps90 (biome ID),
#       PC1_env2
#
# Outputs:
#   gee_exports/sample_points_with_composite_UUID.csv
#   gee_exports/sample_points_for_GEE_UUID.csv
#
# Required packages (install beforehand):
#   uuid, tidyverse, sf, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(uuid)        # unique IDs
library(tidyverse)   # dplyr, ggplot2, ...
library(sf)          # parse the point geometry
library(here)        # project-relative paths (run from the project root)

overwrite_ids <- FALSE   # TRUE = allow replacing existing IDs (breaks links to later steps)

path_input <- here("gee_exports", "sample_points_with_composite.csv")
path_full  <- here("gee_exports", "sample_points_with_composite_UUID.csv")
path_gee   <- here("gee_exports", "sample_points_for_GEE_UUID.csv")

if (file.exists(path_full) && !overwrite_ids) {
  stop("IDs already exist (", path_full, "). Set overwrite_ids <- TRUE only if you ",
       "intend to replace them; later steps depend on the existing IDs.")
}

biome_lookup <- tibble(
  biome_num  = as.character(1:12),
  biome_name = c("Patagonian", "Boreal", "Indomalayan", "Warm Temperate",
                 "Mesoamerican-Andean", "Sahelian-Subequatorial", "East African",
                 "African Rainforest", "Cool Temperate", "Temperate-Mediterranean",
                 "Brazilian Tropical", "Amazon Rainforest")
)


# ---- 2. Load data and prepare columns ----------------------------------------

df <- read.csv(path_input) %>%
  rename(forest_mask = first,
         biome_num   = grps90)

# Longitude and latitude from the point geometry (WKT)
coords <- st_coordinates(st_as_sfc(df$geometry))
df$lon <- coords[, 1]
df$lat <- coords[, 2]


# ---- 3. Assign IDs and save --------------------------------------------------

df$id <- UUIDgenerate(n = nrow(df))
stopifnot(!anyDuplicated(df$id))

# Small table for the NDVI extraction in GEE
df_small <- df %>% select(id, lat, lon)

write.csv(df,       file = path_full, row.names = FALSE)
write.csv(df_small, file = path_gee,  row.names = FALSE)


# ---- 4. Check the PC1 distribution per biome ---------------------------------

df_plot <- df %>%
  mutate(biome_num = as.character(biome_num)) %>%   # same type as in the look-up
  left_join(biome_lookup, by = "biome_num") %>%
  mutate(biome_name = factor(biome_name, levels = biome_lookup$biome_name))

p_check <- ggplot(df_plot, aes(x = PC1_env2, fill = biome_name)) +
  geom_density(alpha = 0.7, color = "black") +
  facet_wrap(~ biome_name, ncol = 4, scales = "free") +
  labs(x = "PC1", y = "Density", title = "Distribution of PC1_env2 per biome") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

print(p_check)
