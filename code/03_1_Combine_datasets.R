# ==============================================================================
# Combine datasets into one table
#
# Description:
#   Left-joins the resilience metrics (lag-1 autocorrelation and variance) and the
#   vegetation sensitivity onto the table of sample points, matching on the point
#   ID. All rows of the sample-point table are kept; points without a match in
#   the other tables get NA.
#
#   data.table is used because the tables are large. The datasets cannot simply
#   be bound by row (rbind / bind_rows), as the order of the IDs can differ and
#   rows can be missing between files, so they are matched on the ID instead.
#
# Workflow:
#   1. Setup
#   2. Load the datasets
#   3. Join on ID
#   4. Save
#
# Inputs:
#   gee_exports/sample_points_with_composite_UUID.csv         sample points (left table, column id)
#   AC1VAR_results.csv                                        temporal autocorrelation and variance (column id)
#   dataSensitivity/robust_results/sensitivity/summary.csv    sensitivity (columns ID2, Sens)
#                                                   
#
# Outputs:
#   01_tidy_data/all_data_combined.csv
#
# Required packages (install beforehand):
#   data.table, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(data.table)   # fast csv reading/writing and joins
library(here)         # project-relative paths (run from the project root)

path_points      <- here("sample_points_with_composite_UUID.csv")
path_ac1var      <- here("dataSensitivity", "AC1VAR_results.csv")
path_sensitivity <- here("dataSensitivity", "robust_results", "sensitivity", "summary.csv")
path_out         <- here("01_tidy_data", "all_data_combined.csv")

dir.create(dirname(path_out), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Load datasets --------------------------------------------------------

# Left table: the dataset all others are matched onto
points_dt <- fread(path_points)

# Autocorrelation (AC1) and variance
ac1var <- fread(path_ac1var)

# Sensitivity: keep the point ID (renamed from ID2) and the sensitivity value
sensitivity <- fread(path_sensitivity)
setnames(sensitivity, "ID2", "id")
sensitivity <- sensitivity[, .(id, Sens)]

# The IDs of the tables that are joined on must be unique, otherwise the join
# would duplicate rows
stopifnot(!anyDuplicated(ac1var$id), !anyDuplicated(sensitivity$id))


# ---- 3. Join on ID -----------------------------------------------------------

combined <- merge(points_dt, ac1var,      by = "id", all.x = TRUE, sort = FALSE)
combined <- merge(combined,  sensitivity, by = "id", all.x = TRUE, sort = FALSE)

# A left join must keep every row of the sample-point table
stopifnot(nrow(combined) == nrow(points_dt))

message("Combined table: ", nrow(combined), " rows, ", ncol(combined), " columns")


# ---- 4. Save -----------------------------------------------------------------

fwrite(combined, path_out)
