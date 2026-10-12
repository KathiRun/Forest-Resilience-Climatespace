# ==============================================================================
# Build the input list for the sensitivity calculation
#
# Description:
#   Combines the monthly NDVI and climate series of every sample point into the
#   list used by the sensitivity pipeline: one item per point, containing one
#   vector per variable, ordered by date.
#     ndvi.ts  kNDVI = tanh(NDVI^2) after rescaling the raw NDVI (x 10000)
#     temp.ts  mean temperature (degrees C)
#     rad.ts   solar radiation
#     wat.ts   water availability = actual / potential evapotranspiration
#
# Workflow:
#   1. Setup and parameters
#   2. Load data
#   3. Calculate water availability
#   4. Join all variables
#   5. Build and save the list
#
# Inputs (dataSensitivity/):
#   NDVI_long.csv                   columns: id, date, NDVI (raw scale, x 10000)
#   gee_exports/climate_AET.csv     columns: id, date, aet
#   gee_exports/climate_PET.csv     columns: id, date, pet
#   gee_exports/climate_TEMP.csv    columns: id, date, tmean_C
#   gee_exports/climate_RAD.csv     columns: id, date, srad
#
# Outputs:
#   dataSensitivity/Data_list.RData       object `final_list` (named by point id)
#   dataSensitivity/Data_list_test.RData  only if test_n_sites is set
#
# Note: the sensitivity pipeline matches its results to the study-site table by
# position, so the order of the points (ascending id) must equal the row order of
# that table. final_list is named by id to allow checking this.
#
# Required packages (install beforehand):
#   data.table, here
# ==============================================================================


# ---- 1. Setup ----------------------------------------------------------------

library(data.table)   # fast csv reading, joins and grouping
library(here)         # project-relative paths (run from the project root)

# ---- Parameters
ndvi_scale   <- 10000   # raw NDVI is stored multiplied by 10000
test_n_sites <- NULL    # e.g. 100 = build a small test list from 100 random points
seed         <- 123     # random seed for the test subset

# ---- Paths
data_dir <- here("dataSensitivity")   # output list
gee_dir  <- here("gee_exports")       # climate exports
path_out <- file.path(data_dir, if (is.null(test_n_sites)) "Data_list.RData" else "Data_list_test.RData")

# Read the selected columns of a csv file; dates stay as text so that they can be
# matched exactly across files
read_input <- function(dir, file, cols) {
  fread(file.path(dir, file), select = cols, colClasses = c(date = "character"))
}


# ---- 2. Load data ------------------------------------------------------------

aet  <- read_input(gee_dir,  "climate_AET.csv",  c("id", "date", "aet"))
pet  <- read_input(gee_dir,  "climate_PET.csv",  c("id", "date", "pet"))
temp <- read_input(gee_dir,  "climate_TEMP.csv", c("id", "date", "tmean_C"))
rad  <- read_input(gee_dir,  "climate_RAD.csv",  c("id", "date", "srad"))
ndvi <- read_input(here(),   "NDVI_long.csv",    c("id", "date", "NDVI"))

# The join keys must be unique, otherwise the joins would duplicate rows
for (tbl in list(aet, pet, temp, rad, ndvi)) {
  stopifnot(!anyDuplicated(tbl, by = c("id", "date")))
}


# ---- 3. Water availability and harmonised variables --------------------------

# Water availability = AET / PET (only dates present in both tables)
wat <- merge(aet, pet, by = c("id", "date"))
wat[, wat.ts := aet / pet]
wat <- wat[, .(id, date, wat.ts)]

# kNDVI transform
ndvi[, ndvi.ts := tanh((NDVI / ndvi_scale)^2)]
ndvi <- ndvi[, .(id, date, ndvi.ts)]

setnames(temp, "tmean_C", "temp.ts")
rad[, rad.ts := as.numeric(srad)]
rad <- rad[, .(id, date, rad.ts)]


# ---- 4. Join all variables ---------------------------------------------------
# Left joins onto the NDVI table: every NDVI row is kept.

all_dt <- merge(ndvi,   temp, by = c("id", "date"), all.x = TRUE)
all_dt <- merge(all_dt, rad,  by = c("id", "date"), all.x = TRUE)
all_dt <- merge(all_dt, wat,  by = c("id", "date"), all.x = TRUE)

stopifnot(nrow(all_dt) == nrow(ndvi))

# Replace NaN (e.g. 0 / 0 in water availability) with NA in all numeric columns
for (col in names(all_dt)) {
  if (is.numeric(all_dt[[col]])) {
    set(all_dt, which(is.nan(all_dt[[col]])), col, NA)
  }
}

# Checks: share of missing values and infinite values per variable
ts_cols <- c("ndvi.ts", "temp.ts", "rad.ts", "wat.ts")
message("Share of missing values per variable:")
print(all_dt[, lapply(.SD, function(x) round(mean(is.na(x)), 4)), .SDcols = ts_cols])
n_inf <- sapply(ts_cols, function(col) sum(is.infinite(all_dt[[col]])))
if (any(n_inf > 0)) {
  warning("Infinite values found (e.g. PET = 0 in water availability): ",
          paste(names(n_inf)[n_inf > 0], n_inf[n_inf > 0], sep = " = ", collapse = ", "))
}

setorder(all_dt, id, date)

# All points must have a series of the same length
rows_per_id <- all_dt[, .N, by = id]$N
if (length(unique(rows_per_id)) != 1) {
  stop("The points do not all have the same number of monthly values")
}

# Optional small test subset
if (!is.null(test_n_sites)) {
  set.seed(seed)
  sampled_ids <- sample(unique(all_dt$id), test_n_sites)
  all_dt <- all_dt[id %in% sampled_ids]
}


# ---- 5. Build and save the list ----------------------------------------------
# One item per point (ordered by id), each with one vector per variable

split_dt <- split(all_dt, by = "id", keep.by = TRUE)

final_list <- lapply(split_dt, function(dt) list(
  ndvi.ts = dt$ndvi.ts,
  temp.ts = dt$temp.ts,
  rad.ts  = dt$rad.ts,
  wat.ts  = dt$wat.ts
))

save(final_list, file = path_out)
message("Saved ", length(final_list), " points (", unique(rows_per_id),
        " months each) to ", path_out)
