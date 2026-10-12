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
POINTS_ASSET = "projects/YOUR_PROJECT_ID/assets/sample_points_for_GEE_UUID"

NDVI_ASSET = "MODIS/061/MOD13A3"

BATCH_SIZE = 10_000
MAX_WORKERS = 6
OUTPUT_CSV = "sample_points_with_ndvi.csv"
OUTPUT_PARQUET = "sample_points_with_ndvi.parquet"
BATCH_PARQUET_DIR = "tmp_batch_parquets"


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
    ndvi_img: ee.Image,
    crs: str,
    transform: list,
    batch_parquet_file: Path,
) -> Path | None:
    """Sample one chunk of points and persist it as a parquet file."""

    sampled = ndvi_img.reduceRegions(
        collection=points_batch,
        reducer=ee.Reducer.first(),
        crs=crs,
        crsTransform=transform,
    )

    try:
        batch_df = geemap.ee_to_gdf(sampled)

    except Exception as exc:
        filtered_size = points_batch.size().getInfo()

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
    
    point_batches = split_feature_collection(points, BATCH_SIZE)

    # -------------------------------------------------------------------------
    # NDVI sampling setup
    # -------------------------------------------------------------------------
    ndvi = ee.ImageCollection(NDVI_ASSET).filterDate("2001-01-01", "2020-12-31")
    ndvi = ndvi.select("NDVI","SummaryQA")

    # Stack every image in the collection into a single multi-band image,
    # one band pair (NDVI, SummaryQA) per date -> wide-format time series per point
    ndvi_img = ndvi.toBands()

    projection = ndvi.first().projection().getInfo()
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
                ndvi_img,
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
