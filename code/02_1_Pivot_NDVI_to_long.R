# ==============================================================================
# Reshape NDVI time series to long format and apply quality control
#
# Description:
#   The NDVI export has one row per sample point and one column per date and
#   band (column names of the form <date>_<band>, e.g. 2001_01_01_NDVI and
#   2001_01_01_SummaryQA). The script
#     (i)  reshapes it to one row per point and date, with one column per band,
#     (ii) applies the pixel-quality rules below, and
#     (iii) saves the result for the resilience calculations.
#
#   Quality rules (SummaryQA identifies the quality of the measured pixel:
#   0 = good, 1 = marginal, 2 = snow/ice, 3 = cloudy):
#     - SummaryQA is floored to an integer class (it sometimes has decimals).
#     - SummaryQA = 2 (snow/ice): NDVI is set to 0.
#     - SummaryQA = 3 (cloudy):   NDVI is set to NA.
#     - Negative NDVI values are not valid vegetation values and are set to NA.
#
# Workflow:
#   1. Setup
#   2. Load the NDVI export
#   3. Reshape to long format (one row per point and date)
#   4. Quality control
#   5. Save
#
# Inputs:
#   gee_exports/sample_points_with_ndvi.csv
#       Columns: id, <date>_<band> columns (and optionally geometry, batch_index)
#
# Outputs:
#   dataSensitivity/NDVI_long.csv
#       Columns: id, date, NDVI, SummaryQA
#
# Required packages (install beforehand):
#   tidyverse, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(tidyverse)   # dplyr, tidyr, stringr, readr
library(here)        # project-relative paths (run from the project root)

path_ndvi <- here("gee_exports", "sample_points_with_ndvi.csv")
path_out  <- here("NDVI_long.csv")

# Columns that are not <date>_<band> values
id_cols <- c("id", "geometry", "batch_index")

dir.create(dirname(path_out), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Load data ------------------------------------------------------------

# check.names = FALSE keeps the column names exactly as exported
NDVI <- read.csv(path_ndvi, check.names = FALSE)

value_cols <- setdiff(names(NDVI), id_cols)


# ---- 3. Reshape to long format -----------------------------------------------
# The band names (NDVI, SummaryQA) contain no underscores, so the band is
# everything after the LAST underscore of a column name and the date is
# everything before it. Each band becomes its own column.

ndvi_long <- NDVI %>%
  select(id, all_of(value_cols)) %>%
  pivot_longer(
    cols          = all_of(value_cols),
    names_to      = c("date_raw", ".value"),
    names_pattern = "(.*)_([^_]+)$"
  ) %>%
  mutate(
    # Drop a leading "X" (added by R to names starting with a digit), then
    # convert 2001_01_01 -> 2001-01-01
    date = as.Date(str_replace_all(str_remove(date_raw, "^X"), "_", "-"))
  ) %>%
  select(-date_raw) %>%
  relocate(date, .after = id) %>%
  arrange(id, date)


# ---- 4. Quality control ------------------------------------------------------

ndvi_long <- ndvi_long %>%
  mutate(
    SummaryQA = floor(SummaryQA),
    NDVI = case_when(
      SummaryQA == 3 ~ NA_real_,   # cloudy pixel
      SummaryQA == 2 ~ 0,          # snow/ice-covered pixel
      NDVI < 0       ~ NA_real_,   # invalid negative NDVI
      TRUE           ~ NDVI
    )
  )


# ---- 5. Save -----------------------------------------------------------------

write_csv(ndvi_long, path_out)
