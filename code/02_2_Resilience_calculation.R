# ==============================================================================
# Resilience metrics: lag-1 autocorrelation (AC1) and variance
#
# Description:
#   For every sample point, the monthly NDVI series is transformed to kNDVI,
#   gap-filled, and decomposed with STL (seasonal-trend decomposition using
#   loess). Both resilience metrics are computed on the STL remainder:
#     - Variance: variance of the remainder
#     - TAC:      lag-1 temporal autocorrelation (Spearman) of the remainder;
#                 shown as "AC1" in the figures
#
#   Processing steps per point:
#     1. kNDVI transform: tanh(NDVI^2), after rescaling the raw NDVI
#        (stored x 10000) to its natural range (Camps-Valls et al. 2021,
#        Science Advances).
#     2. Linear interpolation of gaps of up to `max_gap` months; leading and
#        trailing gaps are filled with the nearest value.
#     3. STL decomposition of the monthly series (frequency 12).
#   Points that cannot be processed (gap longer than `max_gap` months,
#   all-zero series, STL failure or another error) are not part of the output and
#   are listed with the reason in a separate file.
#
# Workflow:
#   1. Setup and parameters
#   2. Per-series calculation
#   3. Load and check the data
#   4. Calculate AC1 and variance for every point
#   5. Save
#
# Inputs:
#   NDVI_long.csv
#       Long format, one row per point and month, ordered by date, complete
#       monthly series starting in January of `start_year`.
#       Columns: id, date, NDVI (raw scale, x 10000)
#
# Outputs:
#   AC1VAR_results.csv   columns: id, Variance, TAC
#   AC1VAR_errors.csv    columns: id, reason (points without a result)
#
# Required packages (install beforehand):
#   data.table, zoo, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(data.table)   # fast csv reading/writing
library(zoo)          # na.approx() for gap filling
library(here)         # project-relative paths (run from the project root)

# ---- Parameters
start_year   <- 2001    # first year of the monthly series (series starts in January)
ndvi_scale   <- 10000   # raw NDVI is stored multiplied by 10000
max_gap      <- 6       # longest gap (months) that is filled by interpolation
stl_t_window <- 23      # STL trend window (months)
stl_s_window <- 13      # STL seasonal window
stl_l_window <- 13      # STL low-pass window
digits_out   <- 10      # decimals kept in the output

# ---- Paths
path_ndvi   <- here("NDVI_long.csv")
path_out    <- here("AC1VAR_results.csv")
path_errors <- here("AC1VAR_errors.csv")


# ---- 2. Per-series calculation -----------------------------------------------

# Returns a list(Variance, TAC, reason). `reason` is NA if the series was
# processed and describes why the series was skipped otherwise.
compute_ac1_var <- function(ndvi_raw, start_year) {
  
  skip <- function(reason) list(Variance = NA_real_, TAC = NA_real_, reason = reason)
  
  # kNDVI transform and gap filling
  ndvi_k      <- tanh((ndvi_raw / ndvi_scale)^2)
  ndvi_filled <- na.approx(ndvi_k, maxgap = max_gap, na.rm = FALSE, rule = 2)
  
  if (anyNA(ndvi_filled))     return(skip(paste("gap longer than", max_gap, "months")))
  if (!any(ndvi_filled != 0)) return(skip("all-zero series"))
  
  # STL decomposition
  ts_i <- ts(ndvi_filled, frequency = 12, start = c(start_year, 1))
  dett <- tryCatch(
    stl(ts_i, t.window = stl_t_window, s.window = stl_s_window,
        l.window = stl_l_window, na.action = na.pass),
    error = function(e) NULL
  )
  if (is.null(dett)) return(skip("STL failed to converge"))
  
  remainder <- dett$time.series[, "remainder"]
  
  list(
    Variance = var(remainder),
    TAC      = cor(remainder, c(NA, remainder[-length(remainder)]),
                   use = "complete.obs", method = "spearman"),
    reason   = NA_character_
  )
}


# ---- 3. Load and check the data ----------------------------------------------

dt <- fread(path_ndvi, select = c("id", "date", "NDVI"))
setorder(dt, id, date)

# The calculation assumes a complete monthly series starting in January of
# start_year for every point
if (format(min(dt$date), "%Y-%m") != paste0(start_year, "-01")) {
  stop("The series do not start in January ", start_year, " (first date: ", min(dt$date), ")")
}
n_months <- (year(max(dt$date)) - year(min(dt$date))) * 12 +
  month(max(dt$date)) - month(min(dt$date)) + 1
rows_per_id <- dt[, .N, by = id]$N
if (any(rows_per_id != n_months)) {
  stop("Not all points have a complete monthly series of ", n_months, " months")
}


# ---- 4. Calculate AC1 and variance for every point ---------------------------

series <- split(dt$NDVI, dt$id)   # one NDVI vector per point
ids    <- names(series)

timing <- system.time({
  results <- rbindlist(lapply(seq_along(series), function(i) {
    res <- tryCatch(
      compute_ac1_var(series[[i]], start_year),
      error = function(e) {
        list(Variance = NA_real_, TAC = NA_real_,
             reason = paste("error:", conditionMessage(e)))
      }
    )
    data.table(id = ids[i], Variance = res$Variance, TAC = res$TAC, reason = res$reason)
  }))
})

out <- results[is.na(reason),
               .(id, Variance = round(Variance, digits_out), TAC = round(TAC, digits_out))]
err_out <- results[!is.na(reason), .(id, reason)]


# ---- 5. Save -----------------------------------------------------------------

fwrite(out, path_out)
fwrite(err_out, path_errors)

message(sprintf("Processed %d points (%d with results, %d without) in %.1f s (%.4f s per point)",
                length(ids), nrow(out), nrow(err_out),
                timing["elapsed"], timing["elapsed"] / length(ids)))
if (nrow(err_out)) {
  message("Points without a result, by reason:")
  print(err_out[, .N, by = reason])
  message("Written to ", path_errors)
}

print(head(out, 10))
