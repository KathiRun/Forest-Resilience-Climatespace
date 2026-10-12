# ==============================================================================
# PCA of bioclimatic variables per functional biome
#
# Description:
#   For each functional biome (grps90 = 1-12), a PCA (centred and scaled) is run
#   on the climate variables of points covering the whole climate space of the
#   biome. The script saves (i) a PCA biplot per biome, (ii) the fitted PCA
#   model, (iii) the PC1/PC2 scores of every point and (iv) the PC loadings of
#   all biomes in one table.
#
# Workflow:
#   1. Setup, variable sets and display names
#   2. Paths and parameters
#   3. PCA, biplot and scores per biome
#   4. Combine and save PCA loadings
#
# Inputs:
#   by_grps90/points_clim_grps90_<g>.csv   (one file per biome)
#       Columns used: ID, lon, lat, grps90 and the climate variables below
#
# Outputs:
#   03_plots/PCA/PCA_biplot_<g>_<env_set>.png
#   04_table_outputs/PCA/PCA_loadings.csv
#   05_data_outputs/PCA/PCA_model_<g>_<env_set>.rds
#   05_data_outputs/PCA/points_grps90_<g>_PCA_scores.csv
#
#
# Required packages (install beforehand):
#   data.table, ggfortify, ggrepel, tidyverse, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(data.table)   # fast csv reading and writing
library(ggfortify)    # autoplot() for prcomp objects
library(ggrepel)      # non-overlapping loading labels (used by ggfortify)
library(tidyverse)    # dplyr, ggplot2, ...
library(here)         # project-relative paths

plot_font <- "Helvetica"

# Climate variable sets. Each set defines one PCA.
env_sets <- list(
  env1 = c(
    "CHELSA_BIO_Annual_Mean_Temperature",
    "WorldClim2_SolarRadiation_AnnualMean",
    "CHELSA_BIO_Annual_Precipitation",
    "CHELSA_exBIO_GrowingSeasonLength"
  ),
  env2 = c(
    "CHELSA_BIO_Annual_Mean_Temperature",
    "WorldClim2_SolarRadiation_AnnualMean",
    "CHELSA_BIO_Annual_Precipitation",
    "CHELSA_exBIO_GrowingSeasonLength",
    "CHELSA_exBIO_AridityIndex",
    "CHELSA_BIO_Mean_Diurnal_Range",
    "CHELSA_BIO_Precipitation_Seasonality",
    "CHELSA_BIO_Min_Temperature_of_Coldest_Month"
  ),
  env3 = c(
    "CHELSA_BIO_Annual_Mean_Temperature",
    "WorldClim2_SolarRadiation_AnnualMean",
    "CHELSA_BIO_Annual_Precipitation",
    "CHELSA_exBIO_GrowingSeasonLength",
    "CHELSA_exBIO_AridityIndex",
    "CHELSA_BIO_Mean_Diurnal_Range",
    "CHELSA_BIO_Precipitation_Seasonality",
    "CHELSA_BIO_Min_Temperature_of_Coldest_Month",
    "SG_Soil_pH_H2O_030cm"
  )
)

# Short names used in biplots and loadings table
display_names <- c(
  "CHELSA_BIO_Annual_Mean_Temperature"          = "MAT",
  "WorldClim2_SolarRadiation_AnnualMean"        = "SR",
  "CHELSA_BIO_Annual_Precipitation"             = "MAP",
  "CHELSA_exBIO_GrowingSeasonLength"            = "GSL",
  "CHELSA_exBIO_AridityIndex"                   = "Aridity",
  "CHELSA_BIO_Mean_Diurnal_Range"               = "MDTR",
  "CHELSA_BIO_Precipitation_Seasonality"        = "PS",
  "CHELSA_BIO_Min_Temperature_of_Coldest_Month" = "MINT",
  "SG_Soil_pH_H2O_030cm"                        = "pH"
)

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
      plot.title        = element_text(size = 7, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = 7, color = "black"),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(1, 1, 1, 0, "mm"),
      legend.position   = "none"
    )
}


# ---- 2. Paths and parameters -------------------------------------------------

# Parameters
groups_to_run   <- 1:12     # biome IDs (grps90)
env_sets_to_run <- "env2"   # which variable sets to run (names of env_sets)

# Paths
input_dir       <- here("by_grps90")
output_plot_PCA <- here("03_plots", "PCA")
output_tab_PCA  <- here("04_table_outputs", "PCA")
output_data_PCA <- here("05_data_outputs", "PCA")

for (d in c(output_plot_PCA, output_tab_PCA, output_data_PCA)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Columns to read from each input file (identifiers + variables of the active sets)
needed_cols <- c("lon", "lat", "grps90", "ID",
                 unique(unlist(env_sets[env_sets_to_run])))


# ---- 3. PCA, biplot and scores per biome -------------------------------------

pca_loadings_list <- list()

for (g in groups_to_run) {
  
  message("Running PCA for biome (grps90): ", g)
  
  group_file <- file.path(input_dir, paste0("points_clim_grps90_", g, ".csv"))
  if (!file.exists(group_file)) stop("Input file not found: ", group_file)
  
  df_group <- fread(group_file, select = needed_cols)
  
  # Table collecting the PC scores of this biome
  group_scores <- df_group[, .(ID, lon, lat, grps90)]
  
  for (env_name in env_sets_to_run) {
    
    env     <- env_sets[[env_name]]
    pc1_col <- paste0("PC1_", env_name)
    pc2_col <- paste0("PC2_", env_name)
    
    # PCA on complete cases only
    df_env    <- as.data.frame(df_group[, ..env])
    rows_used <- complete.cases(df_env)
    
    df_pca <- df_env[rows_used, ]
    colnames(df_pca) <- display_names[env]
    
    pca_model <- prcomp(df_pca, center = TRUE, scale. = TRUE)
    saveRDS(pca_model,
            file = file.path(output_data_PCA, paste0("PCA_model_", g, "_", env_name, ".rds")))
    
    # Scores (NA for points with missing climate data)
    set(group_scores, j = c(pc1_col, pc2_col), value = NA_real_)
    set(group_scores, i = which(rows_used), j = pc1_col, value = pca_model$x[, 1])
    set(group_scores, i = which(rows_used), j = pc2_col, value = pca_model$x[, 2])
    
    # Loadings
    pca_loadings_list[[paste(env_name, g, sep = "_")]] <- data.frame(
      Env_set  = env_name,
      Group    = g,
      Variable = unname(display_names[env]),
      PC1      = pca_model$rotation[, 1],
      PC2      = pca_model$rotation[, 2]
    )
    
    # Biplot
    p <- autoplot(pca_model,
                  loadings              = TRUE,
                  loadings.label        = TRUE,
                  loadings.label.repel  = TRUE,
                  loadings.colour       = "#377EB8",
                  loadings.label.colour = "black",
                  loadings.label.size   = 2.5,
                  loadings.label.family = plot_font,
                  colour = "grey80",
                  alpha  = 0.3,
                  size   = 0.8) +
      ggtitle(paste("Group", g))
    
    # Symmetric y-axis around 0, based on the plotted range
    y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
    y_max   <- max(abs(y_range))
    
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
      scale_y_continuous(limits = c(-y_max, y_max),
                         expand = expansion(mult = 0, add = 0)) +
      coord_fixed() +
      theme_nature()
    
    ggsave(file.path(output_plot_PCA, paste0("PCA_biplot_", g, "_", env_name, ".png")),
           p, width = 88, height = 68, units = "mm", dpi = 900)
    
    rm(df_env, df_pca, pca_model, p)
  }
  
  # Save the scores of this biome
  fwrite(group_scores,
         file.path(output_data_PCA, paste0("points_grps90_", g, "_PCA_scores.csv")))
  
  # Free memory before the next biome
  rm(df_group, group_scores)
  gc()
}


# ---- 4. Combine and save PCA loadings ----------------------------------------

pca_loadings_all <- bind_rows(pca_loadings_list)
rownames(pca_loadings_all) <- NULL

write.csv(pca_loadings_all,
          file = file.path(output_tab_PCA, "PCA_loadings.csv"),
          row.names = FALSE)
