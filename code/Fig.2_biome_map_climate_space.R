# ==============================================================================
# Figure: Global map of functional biomes shaded by distance from 
# PC1 bioclimatic center
# 
# Description:
#   To generate the map in Figure 2.
#   Samples are aggregated into 0.5-degree grid cells. Each cell is coloured by
#   its functional biome. Opacity scales with the mean absolute value of
#   climatic PC1 (PC1_env2) within the biome: cells close to 0 are transparent,
#   cells far from 0 (in either direction) are opaque. The opacity rescaling is
#   done separately for each biome.
#
# Inputs:
#   PC1_sample10k_allgroups.csv
#   functional_biomes_1_12.shp (Functional biome polygons (attribute grps90 = biome ID))
#
# Outputs:
#   Climatespace_map.jpg
#   Climatespace_map.svg
#
# ==============================================================================

# ---- 1. Setup ----------------------------------------------------------------

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(here)

# Paths (relative to the project root)
path_samples <- here("00_raw_data", "PC1_sample10k_allgroups.csv")
path_biomes  <- here("00_raw_data", "Shp_files_functional_biomes",
                     "functional_biomes_1_12.shp")
output_plot  <- here("03_plots")

if (!dir.exists(output_plot)) dir.create(output_plot, recursive = TRUE)

# Analysis parameters
bin_size    <- 0.5            # grid cell size (degrees)
alpha_range <- c(0.05, 1)     # min/max opacity of grid cells
excluded_id <- 1              # biome ID excluded from the figure (Patagonian)

# ---- 2. Biome look-up table, display order and colours -----------------------

biome_lookup <- tibble(
  grps90     = factor(1:12),
  biome_name = c("Patagonian", "Boreal", "Indomalayan", "Warm Temperate",
                 "Mesoamerican-Andean", "Sahelian-Subequatorial",
                 "East African", "African Rainforest", "Cool Temperate",
                 "Temperate-Mediterranean", "Brazilian Tropical",
                 "Amazon Rainforest")
)

# Legend order (Patagonian is excluded)
biome_order <- c(
  "Boreal", "Cool Temperate", "Temperate-Mediterranean", "Warm Temperate",
  "East African", "Sahelian-Subequatorial", "Brazilian Tropical",
  "Mesoamerican-Andean", "Amazon Rainforest", "Indomalayan",
  "African Rainforest"
)

biome_colors <- c(
  "African Rainforest"      = "#1A6B3C",
  "Amazon Rainforest"       = "#74C476",
  "Indomalayan"             = "green1",
  "Sahelian-Subequatorial"  = "#D62728",
  "East African"            = "#E07B39",
  "Brazilian Tropical"      = "#F5C242",
  "Warm Temperate"          = "#7B4FA6",
  "Temperate-Mediterranean" = "#3A7DC9",
  "Mesoamerican-Andean"     = "#17BECF",
  "Cool Temperate"          = "#006D8F",
  "Boreal"                  = "#5C4033"
)

# ---- 3. Load and prepare sample data -----------------------------------------

df_samples <- read.csv(path_samples) %>%
  mutate(grps90 = factor(grps90, levels = levels(biome_lookup$grps90))) %>%
  left_join(biome_lookup, by = "grps90") %>%
  filter(grps90 != excluded_id)

# Aggregate to grid cells and compute opacity
# (mean PC1 per cell, rescaled to alpha_range within each biome)
df_binned <- df_samples %>%
  mutate(
    lon_bin = round(lon / bin_size) * bin_size,
    lat_bin = round(lat / bin_size) * bin_size,
    abs_pc1 = abs(PC1_env2)
  ) %>%
  group_by(biome_name, lon_bin, lat_bin) %>%
  summarise(mean_abs_pc1 = mean(abs_pc1), .groups = "drop") %>%
  group_by(biome_name) %>%
  mutate(alpha_fill = rescale(mean_abs_pc1, to = alpha_range)) %>%
  ungroup() %>%
  mutate(biome_name = factor(biome_name, levels = biome_order))

# ---- 4. Spatial layers -------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

biome_shp <- st_read(path_biomes, quiet = TRUE) %>%
  filter(grps90 != excluded_id)

# ---- 5. Map ------------------------------------------------------------------

p <- ggplot() +
  # Base map
  geom_sf(data = world, fill = "white", color = "grey90", linewidth = 0.2) +
  # Grid cells: colour = biome, opacity = distance from climatic center
  geom_tile(
    data = df_binned,
    aes(x = lon_bin, y = lat_bin, fill = biome_name, alpha = alpha_fill),
    width  = bin_size,
    height = bin_size
  ) +
  # Biome outlines
  geom_sf(data = biome_shp, fill = NA, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = biome_colors, name = "Biome") +
  scale_alpha_continuous(range = alpha_range, guide = "none") +
  coord_sf(expand = FALSE, ylim = c(-60, 85)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 12),
    legend.key       = element_rect(fill = "white"),
    legend.key.size  = unit(0.6, "cm"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    plot.margin      = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
  ) +
  # Show legend keys at full opacity
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

print(p)

# ---- 6. Save -----------------------------------------------------------------

ggsave(
  filename = file.path(output_plot, "Climatespace_map.jpg"),
  plot = p, width = 12, height = 5, units = "in", dpi = 900, bg = "white"
)

ggsave(
  filename = file.path(output_plot, "Climatespace_map.svg"),
  plot = p, width = 180, height = 130, units = "mm", dpi = 900, bg = "white"
)
