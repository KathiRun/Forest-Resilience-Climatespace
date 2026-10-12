# ==============================================================================
# Edge-enriched sampling of points along PC1 (equal-width bins)
#
# Description:
#   For each functional biome (grps90 = 1-12), points are sampled evenly across
#   the biome's PC1 axis. PC1 is divided into equal-width bins and up to
#   `n_per_bin` points are drawn at random from each bin.
#
#   Equal-width bins are used instead of equal-frequency bins (e.g. deciles).
#   Deciles contain the same number of points by construction, so a fixed quota
#   per decile is just a plain random sample and does not enrich the tails.
#   With equal-width bins over a robust core range (0.1st-99.9th percentile),
#   sparse edge bins keep (almost) all of their points, while dense central bins
#   are thinned to the cap. The two outermost bins are open-ended, so extreme
#   points are retained without determining the bin width of the other bins.
#
# Workflow:
#   1. Setup and parameters
#   2. Bin and sample the points of each biome (optional diagnostic plots)
#   3. Combine the samples of all biomes into one file
#
# Inputs:
#   05_data_outputs/PCA/points_grps90_<g>_PCA_scores.csv
#       Columns used: lon, lat, grps90, PC1_env2 (PC scores from the PCA step)
#
# Outputs:
#   05_data_outputs/PCA_quantiles_sampled/points_grps90_<g>_PC1_sample10k.csv
#   05_data_outputs/PCA_quantiles_sampled/PC1_sample10k_allgroups.csv
#   03_plots/PCA_quantiles_sampled/...   (only if plot_per_group = TRUE)
#
# Required packages (install beforehand):
#   data.table, ggplot2, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(data.table)   # fast csv reading and grouped sampling
library(ggplot2)      # diagnostic plots
library(here)         # project-relative paths (run from the project root)

plot_font <- "Helvetica"

# Plotting theme (text sizes in pt)
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

# ---- Parameters
groups_to_run  <- 1:12   # biome IDs (grps90)
n_bins         <- 20     # number of PC1 bins (incl. the two open-ended edge bins)
n_per_bin      <- 100    # maximum points per bin (up to n_bins * n_per_bin per biome)
core_quantiles <- c(0.001, 0.999)   # robust core range used for the bin edges
seed           <- 354

# Diagnostic plots per biome (PC1 histogram with bin edges, points per bin before
# and after sampling). The sampled data are written either way.
plot_per_group <- FALSE

# Biome-specific outlier filters, applied to the raw points before binning.
# Biome 2 (Boreal) has one point whose PC1 is far off the rest of the cloud.
outlier_filters <- list(
  "2" = function(df) df[PC1_env2 >= -10]
)

# ---- Paths
input_dir       <- here("05_data_outputs", "PCA")
output_data_dir <- here("05_data_outputs", "PCA_quantiles_sampled")
output_plot_dir <- here("03_plots", "PCA_quantiles_sampled")

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
if (plot_per_group) dir.create(output_plot_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)


# ---- 2. Bin and sample each biome --------------------------------------------

sampled_list <- list()

for (g in groups_to_run) {
  
  group_file <- file.path(input_dir, paste0("points_grps90_", g, "_PCA_scores.csv"))
  if (!file.exists(group_file)) {
    message("Skipping biome ", g, ": PCA scores file not found (", group_file, ")")
    next
  }
  
  df <- fread(group_file, select = c("lon", "lat", "grps90", "PC1_env2"))
  
  outlier_filter <- outlier_filters[[as.character(g)]]
  if (!is.null(outlier_filter)) df <- outlier_filter(df)
  
  # Equal-width bins over the robust core range. The outermost edges are open
  # (-Inf / Inf), so extreme points fall into the first / last bin.
  core_range <- unname(quantile(df$PC1_env2, core_quantiles))
  bin_cuts   <- c(-Inf, seq(core_range[1], core_range[2], length.out = n_bins - 1), Inf)
  
  df[, pc1_bin := cut(PC1_env2, breaks = bin_cuts, include.lowest = TRUE, labels = FALSE)]
  
  # Stratified sample: up to n_per_bin points per bin (all points if fewer)
  df_sample <- df[, .SD[sample(.N, min(.N, n_per_bin))], by = pc1_bin]
  
  bin_counts <- df[, .N, by = pc1_bin][order(pc1_bin)]
  short_bins <- bin_counts[N < n_per_bin]
  if (nrow(short_bins) > 0) {
    message("Biome ", g, ": bins with fewer than ", n_per_bin,
            " points (all points kept):")
    print(short_bins)
  }
  
  # Diagnostic plots
  if (plot_per_group) {
    
    # Raw PC1 distribution with bin edges
    p_hist <- ggplot(df, aes(x = PC1_env2)) +
      geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "grey70", color = NA) +
      geom_density(color = "#377EB8", linewidth = 0.5) +
      geom_vline(xintercept = bin_cuts[2:n_bins], linetype = "dashed",
                 color = "grey40", linewidth = 0.3) +
      labs(title    = paste("Group", g, "- PC1 distribution with equal-width bin edges"),
           subtitle = paste0(n_bins, " bins, edges from the ",
                             core_quantiles[1] * 100, "th to ",
                             core_quantiles[2] * 100, "th percentile"),
           x = "Bioclimatic space (PC1)", y = "Density") +
      theme_nature()
    
    ggsave(file.path(output_plot_dir, paste0("PC1_bin_histogram_", g, ".png")),
           p_hist, width = 88, height = 68, units = "mm", dpi = 300)
    
    # Points per bin before vs after sampling
    sample_counts <- df_sample[, .N, by = pc1_bin][order(pc1_bin)]
    counts_long <- rbindlist(list(
      bin_counts[,    .(pc1_bin, N, stage = "before sampling")],
      sample_counts[, .(pc1_bin, N, stage = "after sampling")]
    ))
    counts_long[, stage := factor(stage, levels = c("before sampling", "after sampling"))]
    
    p_counts <- ggplot(counts_long, aes(x = factor(pc1_bin), y = N, fill = stage)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_y_log10() +
      scale_fill_manual(values = c("before sampling" = "grey70",
                                   "after sampling"  = "#4DAF4A")) +
      labs(title    = paste("Group", g, "- points per PC1 bin, before vs after sampling"),
           subtitle = paste0("Cap = ", n_per_bin, " per bin (log scale)"),
           x = "PC1 bin (equal width)", y = "Point count") +
      theme_nature() +
      theme(legend.position = "top", legend.title = element_blank())
    
    ggsave(file.path(output_plot_dir, paste0("PC1_bin_sample_counts_", g, ".png")),
           p_counts, width = 120, height = 80, units = "mm", dpi = 300)
  }
  
  # Save the sample of this biome
  out_file <- file.path(output_data_dir, paste0("points_grps90_", g, "_PC1_sample10k.csv"))
  fwrite(df_sample, out_file)
  message("Biome ", g, ": wrote ", nrow(df_sample), " points to ", out_file)
  
  sampled_list[[as.character(g)]] <- df_sample
  
  # Free memory before the next biome
  rm(df, df_sample, bin_counts, short_bins)
  gc()
}


# ---- 3. Combine the samples of all biomes ------------------------------------

df_all <- rbindlist(sampled_list)

combined_file <- file.path(output_data_dir, "PC1_sample10k_allgroups.csv")
fwrite(df_all, combined_file)

message("Combined ", length(sampled_list), " biome sample(s) into ", combined_file,
        " (", nrow(df_all), " points)")
