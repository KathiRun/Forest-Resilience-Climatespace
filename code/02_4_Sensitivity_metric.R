# ==============================================================================
# Vegetation sensitivity, autocorrelation and climate weights per study site
#
# Method adapted from the vegetation sensitivity index of Seddon et al. (2016,
# Nature 531:229-232).
#
# Description:
#   For every study site, the monthly NDVI and climate time series (temperature,
#   water availability, solar radiation) of 2001-2020 are analysed to obtain
#     - climate coefficient weights (betas) and lag-1 autocorrelation of NDVI,
#     - the mean of the NDVI series and the variance of the detrended series
#       (using months with a significant climate weight),
#     - the vegetation sensitivity to each climate variable and in total,
#       unweighted ("NW") and weighted by the climate coefficients ("W").
#   Results are collected in one summary table per study site.
#
# Workflow:
#   1. Setup and parameters
#   2. Run the sensitivity model for every site (parallel)
#   3. Extract climate coefficient weights and autocorrelation
#   4. Extract the mean and variance of the anomalies (variance-mean plots)
#   5. Calculate sensitivity
#   6. Summary table per site
#
# Inputs:
#   functions_sensitivity_metrics.R   (model functions)
#   dataSensitivity/Data_list.RData
#       object `final_list`: one item per site, itself a list of monthly vectors
#       (NDVI and each climate variable) from Jan 2001 to Dec 2020. The item
#       names must match the names used in the functions file. The names of
#       the sites (names(final_list)) are carried into `results` and used as
#       the site ID column of the summary table.
#
# Outputs (under dataSensitivity/):
#   robust_results/VSI2001_2020.RData, badPixels2001_2020.RData
#   robust_results/extractedDF/   coefficients, autocorrelation, means, variances
#   robust_results/sensitivity/   sensitivity vectors and summary.csv
#   Fig.pdf, Fig.RData   variance-mean relationships
#
# Note: site IDs in the summary are the names of `results` (saved in
# VSI2001_2020.RData), so `final_list` must be a named list.
#
# Required packages (install beforehand):
#   MASS, sfsmisc, here  (parallel is part of base R)
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(MASS)      # robust regression (used by the model functions)
library(sfsmisc)   # used by the model functions
library(parallel)  # parallel processing
library(here)      # project-relative paths (run from the project root)

# ---- Parameters
init_year <- 2001    # first year of the time series
end_year  <- 2020    # last year of the time series
ndvi_thresh <- 0.32   # kNDVI threshold used by the model
n_iter    <- 100     # number of iterations of the model
n_cores   <- 8       # number of parallel workers (machine dependent)
run_vsi   <- TRUE    # FALSE = skip step 2 and load the saved results
na_codes  <- c("-8888", "-9999")   # missing-value codes in the csv files

period_us   <- paste0(init_year, "_", end_year)   # used in .RData file names
period_dash <- paste0(init_year, "-", end_year)   # used in csv file names

# ---- Paths
data_dir      <- here("dataSensitivity")
results_dir   <- file.path(data_dir, "robust_results")
extracted_dir <- file.path(results_dir, "extractedDF")
sens_dir      <- file.path(results_dir, "sensitivity")

for (d in c(results_dir, extracted_dir, sens_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Model functions (vsi, vsiSafe, compileCoefTiles, Proj, anomCalcNew, ...)
source(file.path("functions_sensitivity_metric.R"))

# ---- Helper functions
# Read the column `x` of a csv file written by the pipeline
read_x <- function(file, n = NULL) {
  x <- as.numeric(read.csv(file, header = TRUE, na.strings = na_codes)$x)
  if (!is.null(n) && length(x) != n) {
    stop("Unexpected number of rows in ", basename(file), ": ", length(x),
         " (expected ", n, ")")
  }
  x
}

# Rescale to 0-100 based on the full range of the values
rescale_0_100 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 100
}

# Save a sensitivity vector, e.g. save_sens(x, "TempW") -> TempWSensitivity2001-2020.csv
save_sens <- function(x, name) {
  write.csv(x, file = file.path(sens_dir, paste0(name, "Sensitivity", period_dash, ".csv")),
            row.names = FALSE)
}


# ---- 2. Run the sensitivity model for every site -----------------------------
# The model is run per site with vsiSafe(), i.e. vsi() wrapped in tryCatch: a
# degenerate site returns all-NA (with the attribute "vsiError") instead of
# aborting the whole cluster job. NDVI is rescaled inside the model functions
# (see the functions file).

vsi_file <- file.path(results_dir, paste0("VSI", period_us, ".RData"))

if (run_vsi) {
  
  load(file.path(data_dir, "Data_list.RData"))   # provides `final_list`
  site_list <- final_list
  site_list <- site_list[1:10]                 # uncomment for a quick test
  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("vsi", "vsiSafe", "safeRobP", "safeConfint", "det.func",
                      "r2ww", "z.func", "colSds"),
                envir = environment())
  clusterEvalQ(cl, { library(MASS); library(sfsmisc) })
  
  # The parameters are passed as arguments because the workers cannot see the
  # variables of the main session
  results <- tryCatch(
    parLapply(cl, site_list, function(pixel, init_year, end_year, ndvi_thresh, n_iter) {
      vsiSafe(pixel, init.y = init_year, end.y = end_year,
              ndviThresh = ndvi_thresh, niter = n_iter)
    }, init_year = init_year, end_year = end_year,
    ndvi_thresh = ndvi_thresh, n_iter = n_iter),
    finally = stopCluster(cl)
  )
  
  # Sites for which the model failed, and why
  bad_pix <- which(!sapply(results, function(x) is.null(attr(x, "vsiError"))))
  message(length(bad_pix), " sites failed out of ", length(results))
  if (length(bad_pix)) {
    print(table(unlist(lapply(results[bad_pix], function(x) attr(x, "vsiError")))))
  }
  
  badPix <- bad_pix   # object name kept for compatibility with the saved file
  save(badPix, file = file.path(results_dir, paste0("badPixels", period_us, ".RData")))
  save(results, file = vsi_file)
  
} else {
  load(vsi_file)   # provides `results`
}


# ---- 3. Climate coefficient weights and autocorrelation ----------------------

# Climate coefficient weights (betas): raw values, anomalies removed
coef_results <- compileCoefTiles(absVal = TRUE, results = results)
Proj(resultMatrix = coef_results, resultPath = file.path(extracted_dir, "coef_"), raw = TRUE)

# Autocorrelation, anomalies removed
autocor      <- compileAutocor(results = results, absVal = TRUE)
autocor_save <- cutAnomalies(autocor)
write.csv(autocor_save, file = file.path(extracted_dir, "autocor.csv"), row.names = FALSE)


# ---- 4. Mean and variance of the anomalies -----------------------------------
# Mean of the initial time series and variance of the detrended series, for the
# months with a significant climate weight.

# Mean NDVI in each variable
mean_sig <- compileVarTilesAnon(results = results, var = "meanVecSig")
Proj(resultMatrix = mean_sig, fullStand = FALSE,
     resultPath = file.path(extracted_dir, "meanVec_"), raw = TRUE)

# Variance
sd_sig <- compileVarTilesAnon(results = results, var = "sdDetVecSig")
Proj(resultMatrix = sd_sig, fullStand = FALSE,
     resultPath = file.path(extracted_dir, "sdDetVec_"), raw = TRUE)

# Anomalies for each variable (written to extractedDF/anom_<var>.csv)
anom_vars <- c("temp", "wat", "rad", "ndvi")
anom <- setNames(lapply(anom_vars, function(v) {
  anomCalcNew(var = v, namefile = file.path(extracted_dir, "anom_"))
}), anom_vars)

# Variance-mean relationship of one variable
make_anom_plot <- function(anom_var, x_label, panel_label) {
  col_list <- list(
    ndvi1 = rgb(100, 100, 100, 70, maxColorValue = 255),
    temp1 = rgb(100,   0,   0, 70, maxColorValue = 255),
    wat1  = rgb(  0,   0, 100, 70, maxColorValue = 255),
    rad1  = rgb(  0, 100,   0, 70, maxColorValue = 255)
  )
  
  plot(anom_var$x, anom_var$y, col = col_list[[anom_var$varName]], pch = 16,
       xlab = x_label, ylab = "Variance")
  lines(anom_var$xPred, anom_var$yPred)
  
  # Panel label in the top-left corner
  x_range <- range(anom_var$x, na.rm = TRUE)
  y_range <- range(anom_var$y, na.rm = TRUE)
  text(x_range[1] + diff(x_range) * 0.025, y_range[1] + diff(y_range) * 0.975,
       labels = panel_label, cex = 1.1)
}

pdf(file.path(data_dir, "Fig.pdf"))
par(mfrow = c(2, 2))
make_anom_plot(anom$ndvi, paste("Mean NDVI for", period_dash, "period"), "a")
make_anom_plot(anom$temp, paste("Mean temperature (°C) for", period_dash, "period"), "b")
make_anom_plot(anom$wat,  paste("Mean WAT for", period_dash, "period"), "c")
make_anom_plot(anom$rad,  paste("Mean solar radiation for", period_dash, "period"), "d")
dev.off()

dataToSave <- anom[c("temp", "wat", "rad", "ndvi")]
save(dataToSave, file = file.path(data_dir, "Fig.RData"))


# ---- 5. Sensitivity ----------------------------------------------------------
# Sensitivity to a climate variable = log10((NDVI variability + 1) /
# (climate variability + 1)), after rescaling both to 0-100. "NW" = not weighted;
# "W" = weighted by the climate coefficient; "r" = rescaled to 0-100.

# Variability (anomaly) files
ndvi_r <- rescale_0_100(read_x(file.path(extracted_dir, "anom_ndvi.csv")))
wat_r  <- rescale_0_100(read_x(file.path(extracted_dir, "anom_wat.csv")))
temp_r <- rescale_0_100(read_x(file.path(extracted_dir, "anom_temp.csv")))
rad_r  <- rescale_0_100(read_x(file.path(extracted_dir, "anom_rad.csv")))

# Climate coefficients (weights, without t-1)
wat_coef  <- read_x(file.path(extracted_dir, "coef_wat.csv"))
temp_coef <- read_x(file.path(extracted_dir, "coef_temp.csv"))
rad_coef  <- read_x(file.path(extracted_dir, "coef_rad.csv"))

# Sensitivity to each climate variable (not weighted) and total
sens_wat_nw   <- log10((ndvi_r + 1) / (wat_r  + 1))
sens_temp_nw  <- log10((ndvi_r + 1) / (temp_r + 1))
sens_rad_nw   <- log10((ndvi_r + 1) / (rad_r  + 1))
sens_total_nw <- sens_wat_nw + sens_temp_nw + sens_rad_nw

# Weighted sensitivities, each rescaled to 0-100, and weighted total
sens_wat_wr   <- rescale_0_100(sens_wat_nw  * wat_coef)
sens_temp_wr  <- rescale_0_100(sens_temp_nw * temp_coef)
sens_rad_wr   <- rescale_0_100(sens_rad_nw  * rad_coef)
sens_total_w  <- sens_wat_wr + sens_temp_wr + sens_rad_wr

# Save
save_sens(sens_total_w,  "W")
save_sens(sens_total_nw, "NW")

save_sens(sens_wat_wr,  "WATW")
save_sens(sens_temp_wr, "TempW")
save_sens(sens_rad_wr,  "RadW")

save_sens(sens_wat_nw,  "WATNW")
save_sens(sens_temp_nw, "TempNW")
save_sens(sens_rad_nw,  "RadNW")


# ---- 6. Summary table per site -----------------------------------------------
# Sensitivity (weighted, rescaled), autocorrelation and climate weights per site.

site_ids <- names(results)   # site names from VSI2001_2020.RData
if (is.null(site_ids)) stop("`results` has no names: final_list must be a named list")
n_sites <- length(site_ids)

summary_df <- data.frame(
  ID2        = site_ids,
  Sens       = read_x(file.path(sens_dir, paste0("WSensitivity",    period_dash, ".csv")), n_sites),
  TempSens   = read_x(file.path(sens_dir, paste0("TempWSensitivity", period_dash, ".csv")), n_sites),
  WATSens    = read_x(file.path(sens_dir, paste0("WATWSensitivity",  period_dash, ".csv")), n_sites),
  RadSens    = read_x(file.path(sens_dir, paste0("RadWSensitivity",  period_dash, ".csv")), n_sites),
  autocorr   = read_x(file.path(extracted_dir, "autocor.csv"),  n_sites),
  WeightWAT  = read_x(file.path(extracted_dir, "coef_wat.csv"),  n_sites),
  WeightTemp = read_x(file.path(extracted_dir, "coef_temp.csv"), n_sites),
  WeightRad  = read_x(file.path(extracted_dir, "coef_rad.csv"),  n_sites),
  WeightT1   = read_x(file.path(extracted_dir, "coef_t1.csv"),   n_sites)
)

write.csv(summary_df, file = file.path(sens_dir, "summary.csv"), row.names = FALSE)
