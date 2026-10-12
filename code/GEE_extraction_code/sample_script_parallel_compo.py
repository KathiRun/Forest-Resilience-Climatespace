import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import ee
import geemap
import geopandas as gpd
import pandas as pd
from loguru import logger
from tqdm import tqdm


# Replace with the path to your own Earth Engine asset containing the sample points.
# Upload PC1_sample10k_allgroups.csv to your assets after running the R code 
# 01_1_Climate_space_PCA.R and 01_2_Stratified_sampling_PC1.R
POINTS_ASSET = "projects/YOUR_PROJECT_ID/assets/PC1_sample10k_allgroups"

# STABLE_MASK_ASSET and FOREST_2000_ASSET are created with Build_forest_mask.js
STABLE_MASK_ASSET = "projects/YOUR_PROJECT_ID/assets/stableMask_LC_type"
FOREST_2000_ASSET = "projects/YOUR_PROJECT_ID/assets/forestPerc_2000"

# Note: this is a non-public composite dataset of global climate data, all projected to 1km resolution
# but the datasets extracted here are also publicly available through CHELSA and WorldClim.
# To run all subsequent codes you can use sample_points_with_composite_UUID.csv, which includes
# the necessary climate variables for our dataset
COMPOSITE_ASSET = "projects/PROJECT_ID/CrowtherLab_Composite_30ArcSec"

COMPOSITE_BANDS = [
    "CHELSA_BIO_Annual_Mean_Temperature",
    "CHELSA_BIO_Mean_Diurnal_Range",
    "CHELSA_exBIO_GrowingSeasonLength",
    "CHELSA_BIO_Min_Temperature_of_Coldest_Month",
    "CHELSA_BIO_Annual_Precipitation",
    "CHELSA_BIO_Precipitation_Seasonality",
    "CHELSA_exBIO_AridityIndex",
    "WorldClim2_SolarRadiation_AnnualMean",
    "SG_Soil_pH_H2O_030cm",
]

BATCH_SIZE = 10_000
MAX_WORKERS = 6
OUTPUT_CSV = "sample_points_with_composite.csv"
OUTPUT_PARQUET = "sample_points_with_composite.parquet"
BATCH_PARQUET_DIR = "tmp_batch_parquets"
# Note: after creating sample_points_with_composite.csv we assign UUIDs to the data
# using 01_3_Add_UUID.R 

def batch_parquet_path(batch_dir: Path, batch_index: int) -> Path:
    """Build deterministic parquet path for one batch."""

    return batch_dir / f"batch_{batch_index:06d}.parquet"


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


def process_batch(
    batch_index: int,
    points_batch: ee.FeatureCollection,
    forest_mask: ee.Image,
    composite: ee.Image,
    crs: str,
    transform: list,
    batch_parquet_file: Path,
) -> Path | None:
    """Sample one chunk of points and persist it as a parquet file."""

    points_with_mask = forest_mask.reduceRegions(
        collection=points_batch,
        reducer=ee.Reducer.first(),
        scale=500,
    )

    filtered_batch = points_with_mask.filter(ee.Filter.eq("first", 1))

    sampled = composite.reduceRegions(
        collection=filtered_batch,
        reducer=ee.Reducer.first(),
        crs=crs,
        crsTransform=transform,
    )

    try:
        batch_df = geemap.ee_to_gdf(sampled)

    except Exception as exc:
        filtered_size = filtered_batch.size().getInfo()

        if filtered_size == 0:
            logger.warning(
                f"Batch {batch_index}: no points after forest-mask filtering, skipping"
            )
            return None

        logger.exception(
            f"Batch {batch_index}: failed to retrieve data "
            f"(filtered_size={filtered_size}): {exc}"
        )
        return None

    if batch_df.empty:
        return None

    batch_df["batch_index"] = batch_index
    batch_df.to_parquet(batch_parquet_file, index=False)

    # logger.info(f"Batch {batch_index}: wrote {len(batch_df)} rows")
    return batch_parquet_file


def main() -> None:
    batch_dir = Path(BATCH_PARQUET_DIR)
    batch_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # Load points
    # -------------------------------------------------------------------------
    points = ee.FeatureCollection(POINTS_ASSET)
    # -------------------------------------------------------------------------
    # Forest mask
    #
    # Requirements:
    #   - tree cover in 2000 >= 30%
    #   - stable forest according to stableMask_LC_type
    # -------------------------------------------------------------------------
    forest_2000 = ee.Image(FOREST_2000_ASSET)
    stable_mask = ee.Image(STABLE_MASK_ASSET)

    forest_mask = forest_2000.gte(30).And(stable_mask.eq(1)).rename("forest_mask")

    point_batches = split_feature_collection(points, BATCH_SIZE)

    # -------------------------------------------------------------------------
    # Composite sampling setup
    # -------------------------------------------------------------------------
    composite = ee.Image(COMPOSITE_ASSET).select(COMPOSITE_BANDS)

    projection = composite.projection().getInfo()
    crs = projection["crs"]
    transform = projection["transform"]

    # -------------------------------------------------------------------------
    # Parallel batch processing
    # -------------------------------------------------------------------------
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
            batch_parquet_file = batch_parquet_path(batch_dir, batch_index)

            if batch_parquet_file.exists():
                batch_parquet_files.append(batch_parquet_file)
                skipped_batches += 1
                continue

            future = executor.submit(
                process_batch,
                batch_index,
                points_batch,
                forest_mask,
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
            desc="Processing batches",
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

    unique_batch_files = sorted(set(batch_parquet_files))

    if unique_batch_files:
        batch_frames = [gpd.read_parquet(path) for path in unique_batch_files]
        final_df = gpd.GeoDataFrame(
            pd.concat(batch_frames, ignore_index=True),
            geometry="geometry",
            crs=batch_frames[0].crs,
        )

        if "id" in final_df.columns:
            final_df = final_df.sort_values("id").reset_index(drop=True)
    else:
        final_df = gpd.GeoDataFrame(geometry=[])

    final_df.to_parquet(OUTPUT_PARQUET, index=False)
    final_df.to_csv(OUTPUT_CSV, index=False)
    logger.info(
        f"Done. Wrote {len(final_df)} rows to {OUTPUT_PARQUET} and {OUTPUT_CSV}"
    )


if __name__ == "__main__":
    ee.Initialize(
        project="gem-eth-analysis",
        # url="https://earthengine-highvolume.googleapis.com",
    )
    main()
