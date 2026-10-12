import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import ee
import geemap
import geopandas as gpd
import pandas as pd
from loguru import logger
from tqdm import tqdm

# Note: after assigning UUIDs to sample_points_with_composite.csv using 01_3_Add_UUID.R
# upload the file to your assets
# Replace with the path to your own Earth Engine asset containing the sample points.
POINTS_ASSET = "projects/YOUR_PROJECT_ID/assets/sample_points_for_GEE_UUID"

CLIMATE_ASSET = "IDAHO_EPSCOR/TERRACLIMATE"

ID_PROPERTY = "id"

START_YEAR = 2001
END_YEAR = 2020

BATCH_SIZE = 1_000
MAX_WORKERS = 6

# variable group -> band(s) it contributes to the final long-format tables
VARIABLES = {
    "PET": ["pet"],
    "AET": ["aet"],
    "RAD": ["srad"],
    "TEMP": ["tmmn", "tmmx", "tmin_C", "tmax_C", "tmean_C", "tmean_K"],
}

BATCH_PARQUET_DIR = Path("tmp_batch_parquets_climate_wide")
OUTPUT_DIR = Path(".")

# Delimiter separating the date prefix from the band name in renamed bands.
# Must not appear inside a band name itself.
DATE_BAND_SEP = "__"


def split_feature_collection(
    fc: ee.FeatureCollection,
    batch_size: int = BATCH_SIZE,
) -> list[ee.FeatureCollection]:
    """Split a FeatureCollection into chunks of `batch_size`."""

    total_points = fc.size().getInfo()
    if total_points == 0:
        return []

    num_batches = math.ceil(total_points / batch_size)
    batches: list[ee.FeatureCollection] = []

    for batch_index in range(num_batches):
        offset = batch_index * batch_size
        current_size = min(batch_size, total_points - offset)
        batch_fc = ee.FeatureCollection(fc.toList(current_size, offset))
        batches.append(batch_fc)

    return batches


def add_derived_temperature_bands(image: ee.Image) -> ee.Image:
    """Add tmin_C, tmax_C, tmean_C, tmean_K bands, mirroring the original JS logic."""

    tmin_c = image.select("tmmn").divide(10).rename("tmin_C")
    tmax_c = image.select("tmmx").divide(10).rename("tmax_C")
    tmean_c = tmin_c.add(tmax_c).divide(2).rename("tmean_C")
    tmean_k = image.select("tmmn").add(image.select("tmmx")).divide(2).rename("tmean_K")

    return image.addBands([tmin_c, tmax_c, tmean_c, tmean_k])


def rename_bands_with_date(image: ee.Image) -> ee.Image:
    """Prefix every band name with its image's date, e.g. 'pet' -> '2001-01-01__pet'.

    This makes the column names produced by .toBands() unambiguous to parse
    back into (date, band) pairs later, regardless of what prefix .toBands()
    itself adds.
    """

    date_str = image.date().format("yyyy-MM-dd")
    band_names = image.bandNames()
    new_names = band_names.map(lambda b: ee.String(date_str).cat(DATE_BAND_SEP).cat(ee.String(b)))
    return image.rename(new_names)


def build_composite() -> tuple[ee.Image, str, list]:
    """Build the single wide multi-band composite covering every variable and date."""

    raw_climate = ee.ImageCollection(CLIMATE_ASSET).filterDate(
        f"{START_YEAR}-01-01", f"{END_YEAR + 1}-01-01"
    )

    climate = ee.ImageCollection(CLIMATE_ASSET).filterDate(
        f"{START_YEAR}-01-01", f"{END_YEAR + 1}-01-01"
    )
    climate = raw_climate.select(["pet", "aet", "srad", "tmmn", "tmmx"])
    climate = climate.map(add_derived_temperature_bands)

    projection = raw_climate.first().select("pet").projection().getInfo()
    print(projection)  # inspect actual keys before assuming 'crs' exists
    
    crs = projection.get("crs") or projection.get("wkt")
    transform = projection["transform"]

    if crs is None:
        raise ValueError(f"Could not determine CRS from projection info: {projection}")

    climate_dated = climate.map(rename_bands_with_date)
    climate_img = climate_dated.toBands()

    return climate_img, crs, transform


def process_batch(
    batch_index: int,
    points_batch: ee.FeatureCollection,
    composite: ee.Image,
    crs: str,
    transform: list,
    batch_parquet_file: Path,
) -> Path | None:
    """Sample one chunk of points against the full wide composite in a single call."""

    sampled = composite.reduceRegions(
        collection=points_batch,
        reducer=ee.Reducer.first(),
        crs=crs,
        crsTransform=transform,
    )

    try:
        batch_df = geemap.ee_to_gdf(sampled)
    except Exception as exc:
        logger.exception(f"Batch {batch_index}: failed to retrieve data: {exc}")
        return None

    if batch_df.empty:
        return None

    batch_df.to_parquet(batch_parquet_file, index=False)
    return batch_parquet_file


def melt_batch_for_variable(
    batch_df: pd.DataFrame,
    id_property: str,
    value_columns: list[str],
) -> pd.DataFrame:
    """Reshape ONE batch's columns for ONE variable into long format.
 
    Doing this per-batch, per-variable keeps every intermediate DataFrame
    small (roughly BATCH_SIZE rows x a handful of columns) instead of ever
    materializing the full points x all-variables x all-dates table at once.
    """

    if id_property not in batch_df.columns:
        raise ValueError(
            f"'{id_property}' not found in batch columns: {list(batch_df.columns)}. "
            "Set ID_PROPERTY at the top of the script to match the actual "
            "property name on your points asset."
        )
 
    cols_for_variable = [
        c
        for c in batch_df.columns
        if DATE_BAND_SEP in c and c.rsplit(DATE_BAND_SEP, 1)[1] in value_columns
    ]
 
    if not cols_for_variable:
        return pd.DataFrame(columns=[id_property, "date", *value_columns])
 
    subset = batch_df[[id_property, *cols_for_variable]]
 
    long_df = subset.melt(id_vars=[id_property], var_name="col", value_name="value")
 
    split_cols = long_df["col"].str.rsplit(DATE_BAND_SEP, n=1, expand=True)
    long_df["band"] = split_cols[1]
    long_df["date"] = split_cols[0].str.extract(r"(\d{4}-\d{2}-\d{2})")
    long_df = long_df.drop(columns="col")
 
    pivoted = long_df.pivot_table(
        index=[id_property, "date"], columns="band", values="value", aggfunc="first"
    ).reset_index()
    pivoted.columns.name = None
 
    return pivoted
 
 
def extract_and_sample(points: ee.FeatureCollection) -> list[Path]:
    """Run (or resume) the batched Earth Engine sampling step. Returns parquet paths."""
 
    BATCH_PARQUET_DIR.mkdir(parents=True, exist_ok=True)
    point_batches = split_feature_collection(points, BATCH_SIZE)
    composite, crs, transform = build_composite()
 
    num_batches = len(point_batches)
    logger.info(
        f"Processing {num_batches} batches in parallel "
        f"(batch_size={BATCH_SIZE}, workers={MAX_WORKERS})"
    )
 
    batch_parquet_files: list[Path] = []
    skipped_batches = 0
 
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {}
 
        for batch_index, points_batch in enumerate(point_batches):
            batch_parquet_file = BATCH_PARQUET_DIR / f"batch_{batch_index:06d}.parquet"
 
            if batch_parquet_file.exists():
                batch_parquet_files.append(batch_parquet_file)
                skipped_batches += 1
                continue
 
            future = executor.submit(
                process_batch,
                batch_index,
                points_batch,
                composite,
                crs,
                transform,
                batch_parquet_file,
            )
            futures[future] = batch_index
 
        logger.info(
            f"Skipped {skipped_batches} batches with existing parquet files. "
            f"Processing {len(futures)} remaining batches."
        )
 
        with tqdm(
            total=num_batches,
            initial=skipped_batches,
            desc="Sampling batches",
            unit="batch",
        ) as progress:
            for future in as_completed(futures):
                batch_index = futures[future]
                try:
                    batch_file = future.result()
                    if batch_file is not None:
                        batch_parquet_files.append(batch_file)
                except ee.EEException as exc:
                    logger.exception(f"Batch {batch_index} failed: {exc}")
                finally:
                    progress.update(1)
 
    return sorted(set(batch_parquet_files))
 
 
def reshape_and_save(unique_batch_files: list[Path]) -> None:
    """Reshape to long format one batch/variable at a time and write per-variable outputs."""
 
    for variable_name, value_columns in VARIABLES.items():
        frames = []
 
        for batch_file in tqdm(unique_batch_files, desc=f"Reshaping {variable_name}"):
            batch_df = gpd.read_parquet(batch_file)
            pivoted = melt_batch_for_variable(batch_df, ID_PROPERTY, value_columns)
            if not pivoted.empty:
                frames.append(pivoted)
 
        if not frames:
            logger.warning(f"[{variable_name}] No data retrieved, skipping output.")
            continue
 
        var_df = pd.concat(frames, ignore_index=True)
        var_df = var_df.rename(columns={ID_PROPERTY: "id"})
        var_df = var_df.sort_values(["id", "date"]).reset_index(drop=True)
 
        out_csv = OUTPUT_DIR / f"climate_{variable_name}.csv"
        out_parquet = OUTPUT_DIR / f"climate_{variable_name}.parquet"
        var_df.to_csv(out_csv, index=False)
        var_df.to_parquet(out_parquet, index=False)
        logger.info(f"[{variable_name}] Wrote {len(var_df)} rows to {out_csv} and {out_parquet}")
 
 
def main() -> None:
    points = ee.FeatureCollection(POINTS_ASSET)
 
    unique_batch_files = extract_and_sample(points)
 
    if not unique_batch_files:
        logger.warning("No data retrieved, nothing to write.")
        return
 
    logger.info("Reshaping to long format (one batch/variable at a time)...")
    reshape_and_save(unique_batch_files)
 
 
if __name__ == "__main__":
    ee.Initialize(
        project="gem-eth-analysis",
        # url="https://earthengine-highvolume.googleapis.com",
    )
    main()