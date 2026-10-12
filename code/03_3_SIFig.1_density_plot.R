# ==============================================================================
# PC1 density across sampling stages, one panel per biome
#
# Description:
#   For each functional biome the density of PC1 is compared across three stages
#   of the sampling pipeline (Patagonian, biome 1, is excluded):
#     grey  = full bioclimatic space: every grid point of the biome
#     green = stratified sampling: sample of the unfiltered space
#             (up to 100 points per equal-width PC1 bin)
#     red   = final dataset after filtering: the analysed subsample
#             (up to 100 points per PC1 bin, up to 2000 per biome)
#   Each density is normalised (area = 1), so the curves are comparable in
#   shape despite the very different numbers of points behind them.
#
# Workflow:
#   1. Setup, parameters and lookups
#   2. Load the final (analysed) dataset
#   3. Collect the PC1 values of each biome for the three stages
#   4. Plot one density panel per biome
#
# Inputs:
#   05_data_outputs/PCA/points_grps90_<g>_PCA_scores.csv
#       full bioclimatic space (PC1_env2)
#   05_data_outputs/PCA_quantiles_sampled/points_grps90_<g>_PC1_sample10k.csv
#       stratified sample (PC1_env2)
#   01_tidy_data/NewSampling/all_data_combined_subsampled_100perbin.csv
#       final analysed subsample (columns biome_num, PC1_env2)
#
# Outputs:
#   03_plots/PCA_quantiles_sampled/PC1_density_comparison_summary_finalSubsample.png
#
# Note: the stratified sample must have been drawn with the same outlier filter
# as applied here to the full space. The loop warns if the sample does not reach
# the lower tail of the full space (sign of an older, stricter filter); in that
# case rerun the equal-width-bin sampling script first.
#
# Required packages (install beforehand):
#   data.table, ggplot2, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(data.table)   # fast csv reading
library(ggplot2)      # plotting
library(here)         # project-relative paths (run from the project root)

plot_font <- "Helvetica"

# Plotting theme (text sizes in pt; larger sizes are set in the plot call)
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = plot_font) +
    theme(
      axis.title        = element_text(size = 7),
      axis.text         = element_text(size = 7, color = "black"),
      axis.line         = element_line(linewidth = 0.4, color = "black"),
      axis.ticks        = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      panel.grid        = element_blank(),
      panel.border      = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.title        = element_text(size = 7, hjust = 0),
      plot.subtitle     = element_text(size = 7, color = "black"),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(1, 1, 1, 1, "mm"),
      legend.position   = "none"
    )
}

# ---- Paths
full_dir        <- here("05_data_outputs", "PCA")          # full space
sample_dir      <- here("05_data_outputs", "PCA_quantiles_sampled")  # stratified sample
final_file      <- here("01_tidy_data", "NewSampling",
                        "all_data_combined_subsampled_100perbin.csv") # analysed subsample
output_plot_dir <- here("03_plots", "PCA_quantiles_sampled")

dir.create(output_plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Parameters
# Patagonian (biome 1) is excluded: too little data to be informative here.
groups_to_run <- 2:12

# Same outlier filter as in the sampling script, applied to the full space
# (biome 2, Boreal, has one point whose PC1 is far off the rest of the cloud).
outlier_filters <- list(
  "2" = function(df) df[PC1_env2 >= -10]
)

# ---- Lookups
biome_lookup <- data.table(
  grps90     = 1:12,
  biome_name = c("Patagonian", "Boreal", "Indomalayan", "Warm Temperate",
                 "Mesoamerican-Andean", "Sahelian-Subequatorial", "East African",
                 "African Rainforest", "Cool Temperate", "Temperate-Mediterranean",
                 "Brazilian Tropical", "Amazon Rainforest")
)

# Panel order
biome_order <- c(
  "Boreal", "Cool Temperate", "Temperate-Mediterranean", "Warm Temperate",
  "East African", "Sahelian-Subequatorial", "Brazilian Tropical",
  "Mesoamerican-Andean", "Amazon Rainforest", "Indomalayan", "African Rainforest"
)

stage_levels <- c("full bioclimatic space", "stratified sampling",
                  "final dataset after filtering")
stage_colors <- c("full bioclimatic space"          = "grey50",
                  "stratified sampling"             = "#4DAF4A",
                  "final dataset after filtering"   = "#E41A1C")


# ---- 2. Load the final dataset -----------------------------------------------
# Loaded once and sliced per biome in the loop below.

final_df <- fread(final_file, select = c("biome_num", "PC1_env2"))
setnames(final_df, "biome_num", "grps90")


# ---- 3. Collect PC1 values per biome and stage -------------------------------

dens_list <- list()

for (g in groups_to_run) {
  
  full_file   <- file.path(full_dir,   paste0("points_grps90_", g, "_PCA_scores.csv"))
  sample_file <- file.path(sample_dir, paste0("points_grps90_", g, "_PC1_sample10k.csv"))
  
  missing_files <- c(full_file, sample_file)[!file.exists(c(full_file, sample_file))]
  if (length(missing_files) > 0) {
    message("Skipping biome ", g, ": missing ", paste(missing_files, collapse = ", "))
    next
  }
  
  # Full bioclimatic space (with outlier filter)
  full_g <- fread(full_file, select = "PC1_env2")[!is.na(PC1_env2)]
  outlier_filter <- outlier_filters[[as.character(g)]]
  if (!is.null(outlier_filter)) full_g <- outlier_filter(full_g)
  
  sample_g <- fread(sample_file, select = "PC1_env2")
  final_g  <- final_df[grps90 == g, .(PC1_env2)]
  
  # The lowest PC1 bin of the sample is open-ended and holds every point below
  # the 0.1st percentile, so the sample's minimum must lie below that percentile.
  # If not, the sample was drawn with a stricter filter than the one applied here.
  if (min(sample_g$PC1_env2) > quantile(full_g$PC1_env2, 0.001)) {
    warning("Biome ", g, ": the stratified sample does not reach the lower tail of ",
            "the full bioclimatic space - it was probably drawn with an older ",
            "outlier filter. Rerun the sampling script before this one.")
  }
  if (nrow(final_g) == 0) warning("Biome ", g, ": no points in ", final_file)
  
  message("Biome ", g, " - n per stage: full bioclimatic space = ", nrow(full_g),
          ", stratified sampling = ", nrow(sample_g),
          ", final dataset after filtering = ", nrow(final_g))
  
  dens_g <- rbindlist(list(
    full_g[,   .(PC1_env2, stage = "full bioclimatic space")],
    sample_g[, .(PC1_env2, stage = "stratified sampling")],
    final_g[,  .(PC1_env2, stage = "final dataset after filtering")]
  ))
  dens_g[, grps90 := g]
  dens_list[[as.character(g)]] <- dens_g
  
  rm(full_g, sample_g, final_g, dens_g)
  gc()
}


# ---- 4. Plot: PC1 density per stage, one panel per biome ---------------------

if (length(dens_list) == 0) {
  message("No biome processed - summary plot skipped.")
} else {
  
  dens_summary <- rbindlist(dens_list)
  dens_summary <- merge(dens_summary, biome_lookup, by = "grps90")
  dens_summary[, stage := factor(stage, levels = stage_levels)]
  dens_summary[, biome_name := factor(biome_name, levels = biome_order)]
  
  p_summary <- ggplot(dens_summary, aes(x = PC1_env2, color = stage)) +
    geom_density(linewidth = 0.6, key_glyph = "rect") +
    scale_color_manual(values = stage_colors) +
    facet_wrap(~ biome_name, scales = "free") +
    labs(x = "Bioclimatic space (PC1)", y = "Density", color = NULL) +
    theme_nature(base_size = 14) +
    theme(legend.position   = "top",
          legend.text       = element_text(size = 16),
          legend.key.size   = unit(8, "mm"),
          legend.key.width  = unit(6, "mm"),
          legend.key.height = unit(6, "mm"),
          strip.background  = element_blank(),
          strip.text        = element_text(size = 14, face = "bold", hjust = 0),
          axis.title        = element_text(size = 14),
          axis.text         = element_text(size = 12, color = "black"))
  
  out_file <- file.path(output_plot_dir, "PC1_density_comparison_summary_finalSubsample.png")
  ggsave(out_file, p_summary, width = 400, height = 300, units = "mm", dpi = 300,
         limitsize = FALSE)
  message("Wrote ", out_file)
}
