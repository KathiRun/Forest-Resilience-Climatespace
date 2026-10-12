# ==============================================================================
# Relationship between resilience metrics and bioclimatic space (PC1)
#
# Description:
#   For each functional biome, three resilience metrics are related to the
#   position along the biome-specific bioclimatic axis (PC1_env2):
#     - TAC (temporal autocorrelation at lag 1 (plotted as "AC1"))
#     - Variance
#     - Sens (sensitivity)
#   Points are sampled evenly across PC1 space (stratified subsample per PC1
#   bin). Linear and quadratic relationships are tested, correcting for spatial
#   autocorrelation with a residual autocovariate (RAC).
#
# Workflow:
#   1.  Setup, helper functions, lookups, data preparation
#   2.  Stratified subsampling
#   3.  Spatial autocorrelation (correlograms -> RAC distance thresholds)
#   4.  Model fitting and selection (null / linear / quadratic, nested F-tests)
#   5.  Final response-curve plots, RAC-corrected statistics, parabola vertices
#   7.  Supplementary tables
#
# Inputs:
#   all_data_combined.csv
#   PCA_loadings.csv
#
# Outputs:
#   Figure 3-6 and Supplementary Figures 2-12
#   Supplementary Tables 3-6
#
# Required packages (install beforehand):
#   here, tidyverse, patchwork, scales, spdep, gridExtra, ggtext, flextable,
#   ncf, png, sp
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(here)        # project-relative paths
library(tidyverse)   # dplyr, ggplot2, tidyr, purrr, ...
library(patchwork)   # combine ggplots
library(scales)      # custom axis transformations
library(spdep)       # spatial weights, residual autocovariate, Moran's I
library(gridExtra)   # arrange raster grobs (correlogram grid)
library(grid)        # rasterGrob
library(ggtext)      # markdown in ggplot text
library(flextable)   # formatted supplementary tables
# Also used via '::': ncf, png, sp

# ---- Analysis parameters
seed             <- 123    # random seed (subsampling, correlograms)
n_per_bin        <- 100    # points sampled per biome and PC1 bin
corr_cutoff      <- 0.1    # correlation below which spatial autocorrelation is ignored
sig_level        <- 0.05   # significance level for model selection
run_correlograms <- TRUE   # FALSE = re-use the saved threshold distances

plot_font <- "Helvetica"


# ---- 2. Helper functions -----------------------------------------------------

# Significance stars and formatted p-values
get_stars <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

format_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", signif(p, 2))
}

# Vectorised p-value formatting for tables (an exact 0 is printed as "0")
format_p_table <- function(p) {
  ifelse(p == 0, "0", formatC(p, format = "f", digits = 5))
}

# Safe file-name component
safe_name <- function(x) gsub("[^A-Za-z0-9]", "_", x)

# Reversed power transformation (expands the y-axis around high values)
rev_power_trans <- function(power = 2) {
  trans_new(
    name      = paste0("rev_power_", power),
    transform = function(x) -(x^power),
    inverse   = function(x) (-x)^(1 / power),
    domain    = c(0, Inf)
  )
}

# Reversed square-root transformation (expands the y-axis around low values)
rev_sqrt_trans <- function() {
  trans_new(
    name      = "rev_sqrt",
    transform = function(x) -sqrt(x),
    inverse   = function(x) (-x)^2,
    domain    = c(0, Inf)
  )
}

# Evenly spaced y-axis breaks (on the transformed scale if a transformation is given)
get_y_breaks <- function(y_values, n_breaks = 8, digits = NULL, trans = NULL) {
  rng <- range(y_values, na.rm = TRUE)
  
  if (!is.null(trans)) {
    trans_rng    <- trans$transform(rng)
    trans_breaks <- seq(trans_rng[1], trans_rng[2], length.out = n_breaks)
    breaks       <- trans$inverse(trans_breaks)
  } else {
    breaks <- scales::breaks_pretty(n = n_breaks)(rng)
    breaks <- breaks[breaks >= rng[1] & breaks <= rng[2]]
  }
  
  labels <- if (!is.null(digits)) round(breaks, digits) else round(breaks, 6)
  list(breaks = breaks, labels = labels)
}

# Colour palette for a biome. The colour direction depends on the sign of the PC1
# loading of mean annual temperature (cold biomes) or precipitation (other biomes).
get_palette <- function(biome_name, loadings_df) {
  
  is_cold_biome <- biome_name %in% c("Boreal", "Cool Temperate")
  var_needed    <- if (is_cold_biome) "MAT" else "MAP"
  
  loading <- loadings_df %>%
    filter(biome_name == !!biome_name, Variable == var_needed) %>%
    pull(PC1)
  
  if (length(loading) == 0) {
    stop("No PC1 loading found for ", var_needed, " in biome: ", biome_name)
  }
  if (length(loading) > 1) {
    warning("Multiple PC1 loadings found for ", var_needed, " in biome: ",
            biome_name, " - using the first one.")
    loading <- loading[1]
  }
  
  pal <- if (is_cold_biome) {
    c("#313695", "#4575B4", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027", "#A50026")
  } else {
    c("#A50026", "#D73027", "#FDAE61", "#FFFFBF", "#ABD9E9", "#4575B4", "#313695")
  }
  
  if (loading < 0) pal <- rev(pal)
  pal
}

# Common plotting theme (text sizes in pt)
theme_nature <- function(base_size = 7, axis_title = base_size, axis_text = base_size,
                         plot_title = base_size, subtitle = base_size,
                         legend_position = "none") {
  theme_classic(base_size = base_size, base_family = plot_font) +
    theme(
      axis.title        = element_text(size = axis_title),
      axis.text         = element_text(size = axis_text, color = "black"),
      axis.line         = element_line(linewidth = 0.4, color = "black"),
      axis.ticks        = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      panel.grid        = element_blank(),
      panel.border      = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.title        = element_text(size = plot_title, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = subtitle, color = "black"),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(1, 1, 1, 1, "mm"),
      legend.position   = legend_position,
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 6)
    )
}

# Residual autocovariate (RAC): inverse-distance weighted mean of neighbouring
# residuals within `dist` km (great-circle distances)
compute_rac <- function(mod, dat, dist) {
  coords <- cbind(dat$lon, dat$lat)
  autocov_dist(resid(mod), coords,
               nbs         = dist,
               type        = "inverse",
               zero.policy = TRUE,
               style       = "W",
               longlat     = TRUE)
}

# Predictions with 95% confidence interval (t-based)
predict_ci <- function(mod, newdata, level = 0.95) {
  pr    <- predict(mod, newdata = newdata, se.fit = TRUE)
  tcrit <- qt(1 - (1 - level) / 2, df = pr$df)
  newdata$fit <- pr$fit
  newdata$lwr <- pr$fit - tcrit * pr$se.fit
  newdata$upr <- pr$fit + tcrit * pr$se.fit
  newdata
}

# F-test of two nested linear models
nested_anova <- function(m_small, m_large) {
  a <- anova(m_small, m_large)
  list(F = a$F[2], df = a$Df[2], res_df = a$Res.Df[2], p = a$`Pr(>F)`[2])
}

# ggsave() that reports warnings with the file name instead of printing them at the end
ggsave_logged <- function(filename, plot, ...) {
  withCallingHandlers(
    ggsave(filename = filename, plot = plot, ...),
    warning = function(w) {
      message("WARNING while saving ", basename(filename), ": ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}


# ---- 3. Lookups ---------------------------------------------------------------

responses <- c("TAC", "Variance", "Sens")
y_labels  <- c(TAC = "AC1", Variance = "Variance", Sens = "Sensitivity")

# Number of decimals on the y-axis
digits_lookup <- list(Sens = 0, TAC = 2, Variance = 4)

# y-axis transformation per response
y_trans_lookup <- list(
  Sens     = rev_power_trans(power = 2),
  TAC      = rev_sqrt_trans(),
  Variance = rev_sqrt_trans()
)

# Annotations in the response-curve plots: c(upper label, lower label)
direction_labels <- list(
  TAC      = c("\u2191 Recovery faster",     "\u2193 Recovery slower"),
  Variance = c("\u2191 Variability lower",    "\u2193 Variability higher"),
  Sens     = c("\u2191 Compensation higher",  "\u2193 Compensation lower")
)

biome_lookup <- tibble(
  biome_num  = as.character(1:12),
  biome_name = c("Patagonian", "Boreal", "Indomalayan", "Warm Temperate",
                 "Mesoamerican-Andean", "Sahelian-Subequatorial", "East African",
                 "African Rainforest", "Cool Temperate", "Temperate-Mediterranean",
                 "Brazilian Tropical", "Amazon Rainforest")
)

# Plotting order (Patagonian is excluded)
biome_order <- c(
  "Boreal", "Cool Temperate", "Temperate-Mediterranean", "Warm Temperate",
  "East African", "Sahelian-Subequatorial", "Brazilian Tropical",
  "Mesoamerican-Andean", "Amazon Rainforest", "Indomalayan", "African Rainforest"
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

col_linear    <- "#3778c2"
col_quadratic <- "#5c002c"


# ---- 4. Paths ----------------------------------------------------------------

path_data       <- here("01_tidy_data", "all_data_combined.csv")
path_subsample  <- here("01_tidy_data",
                        "all_data_combined_subsampled_100perbin.csv")
path_loadings   <- here("04_table_outputs", "PCA_loadings.csv")

output_plot    <- here("03_plots")
output_tabl    <- here("04_table_outputs")

for (d in c(output_plot, output_tabl,
            file.path(output_plot, responses))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ---- 5. Data preparation -----------------------------------------------------

# PC1 range of each biome (used for x-axis limits and colour scales)
get_pc1_limits <- function(df_biome) {
  data.frame(pc1_min = min(df_biome$PC1_env2, na.rm = TRUE),
             pc1_max = max(df_biome$PC1_env2, na.rm = TRUE))
}

# PCA loadings (for colour palettes), with biome names
loadings_df <- read.csv(path_loadings) %>%
  mutate(Group = as.character(Group)) %>%
  left_join(biome_lookup, by = c("Group" = "biome_num"))

# Metric data: exclude Patagonian, remove NAs, keep non-negative TAC only
df <- read.csv(path_data) %>%
  filter(biome_name != "Patagonian") %>%
  filter(!is.na(PC1_env2), !is.na(TAC), !is.na(Variance), !is.na(Sens)) %>%
  filter(TAC >= 0)

biomes <- biome_order[biome_order %in% unique(df$biome_name)]

# Stratified subsample: n_per_bin points per biome and PC1 bin
# (all points if a bin contains fewer)
set.seed(seed)
df_sub <- df %>%
  group_by(biome_name, pc1_bin) %>%
  group_modify(~ slice_sample(.x, n = min(n_per_bin, nrow(.x)))) %>%
  ungroup()

write.csv(df_sub, file = path_subsample, row.names = FALSE)

# Sanity checks: points per biome, and per biome and bin
df_sub %>% count(biome_name)
df_sub %>% count(biome_name, pc1_bin) %>% print(n = Inf)


# ---- 6. Spatial autocorrelation ----------------------------------------------
# Spline correlograms of the residuals of AC1 ~ PC1 give, per biome, the distance
# (km) at which spatial correlation drops below `corr_cutoff`. This distance is
# used for the RAC of all three metrics.

corr_csv <- file.path(output_tabl, "Correlogram_threshold_distances.csv")

if (run_correlograms) {
  
  set.seed(seed)
  correlog_rows  <- list()
  correlog_grobs <- list()
  
  for (b in biomes) {
    message("Computing correlogram threshold for biome: ", b)
    
    df_biome  <- df_sub %>% filter(biome_name == b)
    mod_ref   <- lm(TAC ~ PC1_env2, data = df_biome)
    
    correlog <- ncf::spline.correlog(
      x      = as.numeric(df_biome$lon),
      y      = as.numeric(df_biome$lat),
      z      = as.numeric(resid(mod_ref)),
      xmax   = 500,
      resamp = 100,
      latlon = TRUE
    )
    
    distances    <- as.vector(correlog$real$predicted$x)
    correlations <- as.vector(correlog$real$predicted$y)
    threshold_distance <- distances[which(correlations < corr_cutoff)[1]]
    
    # Fallback if the correlation never drops below the cutoff
    if (is.na(threshold_distance)) {
      threshold_distance <- max(distances)
      message("  WARNING: correlation never drops below ", corr_cutoff, " for biome ",
              b, ". Using max distance: ", round(threshold_distance, 1), " km")
    }
    message("  Threshold distance: ", round(threshold_distance, 1), " km")
    
    correlog_rows[[b]] <- data.frame(biome = b, threshold_distance = threshold_distance)
    
    # Draw the base-R correlogram to a temporary PNG and keep it as a grob
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = 800, height = 600, res = 150, family = plot_font)
    plot(correlog,
         main = paste0(b, "\n(threshold: ", round(threshold_distance, 1), " km)"),
         cex.main = 0.9, xlab = "Distance (km)", ylab = "Correlation")
    abline(h = corr_cutoff, col = "red",    lty = 2, lwd = 1.5)
    abline(h = 0,           col = "grey40", lty = 1, lwd = 0.8)
    dev.off()
    
    correlog_grobs[[b]] <- rasterGrob(png::readPNG(tmp), interpolate = TRUE)
    unlink(tmp)
  }
  
  correlog_results <- bind_rows(correlog_rows)
  write.csv(correlog_results, file = corr_csv, row.names = FALSE)
  
  grid_plot <- arrangeGrob(grobs = correlog_grobs, ncol = 3, nrow = 4)
  ggsave(file.path(output_plot, "Correlograms.jpg"), plot = grid_plot,
         width = 185, height = 210, units = "mm", dpi = 300)
  
} else {
  correlog_results <- read.csv(corr_csv)
}

threshold_for <- function(b) correlog_results$threshold_distance[correlog_results$biome == b]


# ---- 7. Model fitting and selection ------------------------------------------
# For each response and biome, three RAC-corrected models are fitted:
# null (RAC only), linear (PC1 + RAC) and quadratic (PC1 + PC1^2 + RAC).
# Selection is sequential: the quadratic model is chosen only if it improves
# significantly on the linear model (nested F-test).

# Fit the base model, compute the RAC from its residuals, fit the candidates
fit_candidate_models <- function(df_biome, resp, threshold_distance) {
  mod_base      <- lm(reformulate("PC1_env2", response = resp), data = df_biome)
  df_biome$rac  <- compute_rac(mod_base, df_biome, threshold_distance)
  list(
    df_biome      = df_biome,
    mod_base      = mod_base,
    mod_null      = lm(reformulate("rac", response = resp), data = df_biome),
    mod_linear    = lm(reformulate(c("PC1_env2", "rac"), response = resp),
                       data = df_biome),
    mod_quadratic = lm(reformulate(c("PC1_env2", "I(PC1_env2^2)", "rac"),
                                   response = resp), data = df_biome)
  )
}

# Model-comparison plot: both fits with 95% CI; the selected model is in bold
make_diag_plot <- function(df_biome, resp, b, mod_lin, mod_quad, best_model,
                           y_trans, digits, y_label) {
  
  pc1_range <- range(df_biome$PC1_env2, na.rm = TRUE)
  y_scale   <- get_y_breaks(df_biome[[resp]], digits = digits, trans = y_trans)
  
  newdata <- data.frame(
    PC1_env2 = seq(pc1_range[1], pc1_range[2], length.out = 1000),
    rac      = mean(df_biome$rac, na.rm = TRUE)
  )
  pred_lin  <- predict_ci(mod_lin,  newdata)
  pred_quad <- predict_ci(mod_quad, newdata)
  
  bold_if <- function(txt, is_best) if (is_best) paste0("**", txt, "**") else txt
  subtitle <- paste(
    bold_if(paste0("Linear (blue):    AIC = ", round(AIC(mod_lin), 1)),
            best_model == "linear"),
    bold_if(paste0("Quadratic (purple): AIC = ", round(AIC(mod_quad), 1)),
            best_model == "quadratic"),
    sep = "<br>"
  )
  
  ggplot(df_biome, aes(x = PC1_env2, y = .data[[resp]])) +
    geom_point(alpha = 0.7, colour = "darkgrey", size = 0.8) +
    geom_ribbon(data = pred_lin, aes(x = PC1_env2, ymin = lwr, ymax = upr),
                fill = col_linear, alpha = 0.15, inherit.aes = FALSE) +
    geom_ribbon(data = pred_quad, aes(x = PC1_env2, ymin = lwr, ymax = upr),
                fill = col_quadratic, alpha = 0.15, inherit.aes = FALSE) +
    geom_line(data = pred_lin, aes(x = PC1_env2, y = fit),
              color = col_linear, linewidth = 1, inherit.aes = FALSE) +
    geom_line(data = pred_quad, aes(x = PC1_env2, y = fit),
              color = col_quadratic, linewidth = 1, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5, linetype = "dashed") +
    labs(title = b, subtitle = subtitle, x = "Bioclimatic space (PC1)", y = y_label) +
    scale_x_continuous(limits = pc1_range, oob = scales::censor,
                       expand = expansion(mult = 0.05)) +
    scale_y_continuous(trans = y_trans, breaks = y_scale$breaks,
                       labels = y_scale$labels) +
    theme_nature(base_size = 15, axis_title = 16, axis_text = 15,
                 plot_title = 18, subtitle = 17) +
    theme(plot.subtitle = ggtext::element_markdown(size = 12))
}

fitted_models <- list()   # fitted_models[[response]][[biome]]
anova_all     <- list()   # ANOVA statistics per response
aic_all       <- list()   # AIC / adjusted R2 per response
diag_plots    <- list()   # model-comparison plots

for (resp in responses) {
  message("Fitting models for response: ", resp)
  
  fitted_models[[resp]] <- list()
  diag_plots[[resp]]    <- list()
  anova_rows <- list()
  aic_rows   <- list()
  
  for (b in biomes) {
    
    fit <- fit_candidate_models(df_sub %>% filter(biome_name == b), resp, threshold_for(b))
    fitted_models[[resp]][[b]] <- fit
    
    a_lin  <- nested_anova(fit$mod_null,   fit$mod_linear)      # linear vs null
    a_quad <- nested_anova(fit$mod_linear, fit$mod_quadratic)   # quadratic vs linear
    best_model <- if (a_quad$p < sig_level) "quadratic" else "linear"
    
    anova_rows[[b]] <- data.frame(
      biome            = b,
      n                = nrow(fit$df_biome),
      F_null_lin       = round(a_lin$F, 3),
      df_null_lin      = a_lin$df,
      resdf_null_lin   = a_lin$res_df,
      p_null_lin       = round(a_lin$p, 5),
      stars_null_lin   = get_stars(a_lin$p),
      F_lin_quad       = round(a_quad$F, 3),
      df_lin_quad      = a_quad$df,
      resdf_lin_quad   = a_quad$res_df,
      p_lin_quad       = round(a_quad$p, 5),
      stars_lin_quad   = get_stars(a_quad$p),
      best_model_Anova = best_model
    )
    
    aic_rows[[b]] <- data.frame(
      biome           = b,
      AIC_linear      = AIC(fit$mod_linear),
      AIC_quadratic   = AIC(fit$mod_quadratic),
      AdjR2_linear    = summary(fit$mod_linear)$adj.r.squared,
      AdjR2_quadratic = summary(fit$mod_quadratic)$adj.r.squared
    )
    
    diag_plots[[resp]][[b]] <- make_diag_plot(
      fit$df_biome, resp, b, fit$mod_linear, fit$mod_quadratic, best_model,
      y_trans = y_trans_lookup[[resp]], digits = digits_lookup[[resp]],
      y_label = y_labels[[resp]]
    )
  }
  
  anova_all[[resp]] <- bind_rows(anova_rows)
  aic_all[[resp]]   <- bind_rows(aic_rows) %>%
    mutate(
      best_model_AIC   = if_else(AIC_quadratic < AIC_linear, "quadratic", "linear"),
      best_model_AdjR2 = if_else(AdjR2_quadratic > AdjR2_linear, "quadratic", "linear")
    )
}

# Fitted objects plus the ANOVA-selected model for one response and biome
get_selected_fit <- function(resp, b) {
  fit <- fitted_models[[resp]][[b]]
  fit$selected <- anova_all[[resp]]$best_model_Anova[anova_all[[resp]]$biome == b]
  fit$mod      <- if (fit$selected == "quadratic") fit$mod_quadratic else fit$mod_linear
  fit
}

# Combined result tables
anova_combined <- bind_rows(anova_all, .id = "response")

aic_combined <- bind_rows(aic_all, .id = "response") %>%
  mutate(biome = factor(biome, levels = biome_order),
         response = factor(response, levels = responses)) %>%
  arrange(biome, response) %>%
  mutate(across(c(biome, response), as.character))

write.csv(anova_combined, file.path(output_tabl, "Anova_ClimSpace_ALL.csv"), row.names = FALSE)
write.csv(aic_combined,   file.path(output_tabl, "Best_model_ClimSpace_ALL.csv"), row.names = FALSE)

# Model-comparison panels (one per response, all biomes)
for (resp in responses) {
  
  plots_resp <- diag_plots[[resp]][biomes]
  plots_tagged <- lapply(seq_along(plots_resp), function(i) {
    plots_resp[[i]] +
      labs(tag = letters[i]) +
      theme(plot.tag = element_text(face = "bold", size = 18),
            plot.tag.position = "topleft")
  })
  panel <- wrap_plots(plots_tagged, ncol = 3)
  
  for (ext in c("jpg", "svg")) {
    ggsave_logged(
      file.path(output_plot, resp, paste0(resp, "_all_biomes_diagnostic_panel.", ext)),
      panel, width = 18, height = 22, dpi = 300
    )
  }
}


# ---- 8. Final response curves and RAC-corrected statistics --------------------

# Response-curve plot: points coloured along PC1 (0 at the palette midpoint),
# fitted line with 95% CI, optional parabola vertex and annotations.
make_plot <- function(df_biome, mod, resp, y_label, y_trans, digits,
                            pc1_limits, palette, text_size = 7,
                            vertex = NULL, annotations = NULL,
                            title = "", subtitle = NULL) {
  
  pc1_min    <- pc1_limits$pc1_min
  pc1_max    <- pc1_limits$pc1_max
  y_scale    <- get_y_breaks(df_biome[[resp]], digits = digits, trans = y_trans)
  y_min      <- min(df_biome[[resp]], na.rm = TRUE)
  y_max      <- max(df_biome[[resp]], na.rm = TRUE)
  zero_point <- (0 - pc1_min) / (pc1_max - pc1_min)   # position of PC1 = 0 in [0, 1]
  
  pred_df <- predict_ci(mod, data.frame(
    PC1_env2 = seq(min(df_biome$PC1_env2, na.rm = TRUE),
                   max(df_biome$PC1_env2, na.rm = TRUE), length.out = 1000),
    rac      = mean(df_biome$rac, na.rm = TRUE)
  ))
  
  p <- ggplot(df_biome, aes(x = PC1_env2, y = .data[[resp]], color = PC1_env2)) +
    geom_point(size = 1, alpha = 0.8, shape = 19) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5, linetype = "dashed") +
    # Centre PC1 = 0 on the palette midpoint, spread the other colours evenly
    scale_color_gradientn(
      colors = palette,
      limits = c(pc1_min, pc1_max),
      values = c(0, zero_point / 2, zero_point, (zero_point + 1) / 2, 1)
    ) +
    labs(title = title, subtitle = subtitle, x = "Bioclimatic space (PC1)", y = y_label) +
    scale_x_continuous(limits = c(pc1_min, pc1_max), oob = scales::censor,
                       expand = expansion(mult = 0.05)) +
    scale_y_continuous(trans = y_trans, breaks = y_scale$breaks, labels = y_scale$labels,
                       expand = expansion(mult = 0.05)) +
    coord_cartesian(clip = "off") +
    geom_ribbon(data = pred_df, aes(x = PC1_env2, ymin = lwr, ymax = upr),
                fill = "grey60", alpha = 0.3, inherit.aes = FALSE) +
    geom_line(data = pred_df, aes(x = PC1_env2, y = fit),
              color = "black", linewidth = 0.8, inherit.aes = FALSE) +
    theme_nature(base_size = text_size)
  
  # Parabola vertex: triangle plus a dotted line
  if (!is.null(vertex)) {
    seg_length <- (y_max - y_min) * 0.1
    p <- p +
      geom_segment(
        data = data.frame(x = vertex$x, ystart = vertex$y - seg_length,
                          yend = vertex$y + seg_length),
        aes(x = x, xend = x, y = ystart, yend = yend),
        color = "grey20", linewidth = 0.5, linetype = "dotted", inherit.aes = FALSE
      ) +
      geom_point(data = data.frame(x = vertex$x, y = vertex$y), aes(x = x, y = y),
                 shape = 24, size = 3, fill = "black", color = "black",
                 inherit.aes = FALSE)
  }
  
  # Annotations: overall test statistics (top-left) and direction labels
  if (!is.null(annotations)) {
    txt_size <- text_size / ggplot2::.pt
    p <- p +
      annotate("text", x = pc1_min, y = y_min, label = annotations$p_text,
               hjust = -0.1, vjust = 1, size = txt_size, family = plot_font) +
      annotate("text", x = pc1_min, y = y_min, label = annotations$upper,
               hjust = -0.1, vjust = 6.5, size = txt_size, family = plot_font) +
      annotate("text", x = pc1_min, y = y_max, label = annotations$lower,
               hjust = -0.1, vjust = -1.0, size = txt_size, family = plot_font)
  }
  p
}

# Parabola vertex of a quadratic model. Returns NULL for linear models and for
# parabolas opening downward on the data scale. The y-axes of all three metrics are
# reversed, so an upward-opening parabola appears as a hump in the plots.
get_vertex <- function(mod, df_biome) {
  coefs <- coef(mod)
  if (!"I(PC1_env2^2)" %in% names(coefs) || coefs[["I(PC1_env2^2)"]] <= 0) return(NULL)
  x <- -coefs[["PC1_env2"]] / (2 * coefs[["I(PC1_env2^2)"]])
  y <- as.numeric(predict(mod, newdata = data.frame(
    PC1_env2 = x, rac = mean(df_biome$rac, na.rm = TRUE))))
  list(x = x, y = y)
}

rac_results    <- list()
vertex_results <- list()

for (resp in responses) {
  
  rac_rows    <- list()
  vertex_rows <- list()
  
  for (b in biomes) {
    
    fit      <- get_selected_fit(resp, b)
    df_biome <- fit$df_biome
    mod      <- fit$mod
    b_num    <- biome_lookup$biome_num[biome_lookup$biome_name == b]
    
    # Overall test of the selected model against the RAC-only model
    a_overall <- nested_anova(fit$mod_null, mod)
    
    # Residual spatial autocorrelation after RAC correction (5 nearest neighbours)
    coords      <- cbind(df_biome$lon, df_biome$lat)
    lw          <- nb2listw(knn2nb(knearneigh(coords, k = 5)), style = "W")
    moran_after <- moran.test(resid(mod), lw)
    
    rac_rows[[b]] <- data.frame(
      biome             = b,
      response          = resp,
      selected_model    = fit$selected,
      threshold_dist_km = round(threshold_for(b), 1),
      AdjR2_base        = round(summary(fit$mod_base)$adj.r.squared, 4),
      AdjR2_rac         = round(summary(mod)$adj.r.squared, 4),
      F_overall         = round(a_overall$F, 3),
      df_numerator      = a_overall$df,
      df_denominator    = a_overall$res_df,
      pval_overall      = round(a_overall$p, 5),
      stars_overall     = get_stars(a_overall$p),
      moransI_after     = round(unname(moran_after$estimate["Moran I statistic"]), 4),
      moransI_pval      = round(moran_after$p.value, 5)
    )
    
    vertex <- if (fit$selected == "quadratic") get_vertex(mod, df_biome) else NULL
    if (!is.null(vertex)) {
      vertex_rows[[b]] <- data.frame(
        biome       = b,
        response    = resp,
        PC1_vertex  = round(vertex$x, 3),
        resp_vertex = round(vertex$y, 3),
        PC1_min     = round(min(df_biome$PC1_env2, na.rm = TRUE), 3),
        PC1_max     = round(max(df_biome$PC1_env2, na.rm = TRUE), 3)
      )
    }
    
    p <- make_plot(
      df_biome, mod, resp,
      y_label    = y_labels[[resp]],
      y_trans    = y_trans_lookup[[resp]],
      digits     = digits_lookup[[resp]],
      pc1_limits = get_pc1_limits(df_biome),
      palette    = get_palette(b, loadings_df),
      vertex     = vertex,
      annotations = list(
        p_text = paste0(format_p(a_overall$p), " ", get_stars(a_overall$p)),
        upper  = direction_labels[[resp]][1],
        lower  = direction_labels[[resp]][2]
      )
    )
    
    ggsave_logged(
      file.path(output_plot, resp,
                paste0("ClimateSpace_PC1_", resp, "_", b_num, ".svg")),
      p, width = 88, height = 67, units = "mm", dpi = 900
    )
  }
  
  rac_results[[resp]]    <- bind_rows(rac_rows)
  vertex_results[[resp]] <- bind_rows(vertex_rows)
}

rac_all    <- bind_rows(rac_results)
vertex_all <- bind_rows(vertex_results)

write.csv(rac_all,    file.path(output_tabl, "RAC_corrected_ALL.csv"), row.names = FALSE)
write.csv(vertex_all, file.path(output_tabl, "Vertex_point.csv"),      row.names = FALSE)


# ---- 10. Supplementary tables ------------------------------------------------

# Model selection (nested ANOVA)
supp_table <- anova_combined %>%
  mutate(biome    = factor(biome, levels = biome_order),
         response = factor(response, levels = responses)) %>%
  arrange(biome, response) %>%
  transmute(
    Biome       = as.character(biome),
    Response    = unname(y_labels[as.character(response)]),
    n,
    F_null_lin,
    df_null_lin = paste0(df_null_lin, ", ", resdf_null_lin),
    p_null_lin  = format_p_table(p_null_lin),
    F_lin_quad,
    df_lin_quad = paste0(df_lin_quad, ", ", resdf_lin_quad),
    p_lin_quad  = format_p_table(p_lin_quad),
    Final_model = best_model_Anova
  )

ft_anova <- flextable(supp_table) %>%
  set_header_labels(
    Biome = "Biome", Response = "Response", n = "n",
    F_null_lin = "F", df_null_lin = "df", p_null_lin = "P",
    F_lin_quad = "F", df_lin_quad = "df", p_lin_quad = "P",
    Final_model = "Selected model"
  ) %>%
  add_header_row(
    values    = c("", "", "", "ANOVA: null & linear model",
                  "ANOVA: linear & quadratic model", ""),
    colwidths = c(1, 1, 1, 3, 3, 1)
  ) %>%
  merge_v(j = c("Biome", "n")) %>%   # merge repeated cells vertically
  theme_box() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header") %>%
  autofit()

save_as_docx(ft_anova, path = file.path(output_tabl, "Supplementary_Table_ModelSelection.docx"))
write.csv(supp_table, file.path(output_tabl, "Supplementary_Table_ModelSelection.csv"),
          row.names = FALSE)

# RAC-corrected selected models
rac_table <- rac_all %>%
  mutate(biome    = factor(biome, levels = biome_order),
         response = factor(response, levels = responses)) %>%
  arrange(biome, response) %>%
  transmute(
    Biome                 = as.character(biome),
    Response              = unname(y_labels[as.character(response)]),
    `Selected model`      = selected_model,
    `F(df_num, df_denom)` = paste0(round(F_overall, 3), "(", df_numerator, ", ",
                                   df_denominator, ")"),
    P                     = format_p_table(pval_overall),
    `Adjusted R2`         = round(AdjR2_rac, 4)
  )

ft_rac <- flextable(rac_table) %>%
  merge_v(j = "Biome") %>%
  theme_box() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header") %>%
  italic(j = c("F(df_num, df_denom)", "P"), part = "header") %>%
  autofit()

save_as_docx(ft_rac, path = file.path(output_tabl, "Supplementary_Table_RAC_Results.docx"))
