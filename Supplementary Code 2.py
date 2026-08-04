from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi
from skimage import exposure, filters, morphology


DEFAULT_EXPERIMENT_DIR = Path(r"D:\Dropbox\20260617EXP2NEDAY7_ch1Bcatenin_ch2Nuclei")
DEFAULT_OUTPUT_DIR = Path(r"D:\Projects\chunling\20260617EXP2NEDAY7_analysis_results")
SERIES_RE = re.compile(r"^C(?P<channel>[12])-.*\(series (?P<series>\d+)\)\.tif$", re.IGNORECASE)


def get_stack_shape(path: Path) -> tuple[int, int, int]:
    with Image.open(path) as image:
        z_layers = int(getattr(image, "n_frames", 1))
        width, height = image.size
    if z_layers < 1:
        raise ValueError(f"Expected at least one TIFF plane: {path}")
    return z_layers, height, width


def read_plane(image: Image.Image, z_index: int) -> np.ndarray:
    image.seek(z_index)
    arr = np.asarray(image)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D TIFF plane, got shape {arr.shape}")
    return arr.astype(np.float32, copy=False)


def output_stem_from_tiff(path: Path) -> str:
    return path.stem


def segment_nuclear_region(
    nuclei: np.ndarray,
    min_area: int,
    gaussian_sigma: float,
    local_block_size: int,
    adaptive_threshold: bool,
) -> np.ndarray:
    smooth = filters.gaussian(nuclei, sigma=gaussian_sigma, preserve_range=True)
    global_threshold = filters.threshold_otsu(smooth)
    if adaptive_threshold:
        block_size = local_block_size if local_block_size % 2 == 1 else local_block_size + 1
        local_threshold = filters.threshold_local(smooth, block_size=block_size)
        mask = smooth > np.maximum(local_threshold, global_threshold * 0.65)
    else:
        mask = smooth > global_threshold

    mask = morphology.remove_small_objects(mask, min_size=min_area)
    mask = ndi.binary_fill_holes(mask)
    mask = morphology.remove_small_holes(mask, area_threshold=max(64, min_area // 2))
    return mask.astype(bool)


def pair_channel_files(raw_dir: Path) -> dict[int, dict[str, Path]]:
    pairs: dict[int, dict[str, Path]] = {}
    for path in raw_dir.glob("C[12]-*.tif"):
        match = SERIES_RE.match(path.name)
        if not match:
            continue
        series = int(match.group("series"))
        channel = f"C{match.group('channel')}"
        pairs.setdefault(series, {})[channel] = path
    return {series: files for series, files in sorted(pairs.items()) if {"C1", "C2"} <= files.keys()}


def save_qc_overlay(out_path: Path, nuclei: np.ndarray, nuclear_mask: np.ndarray) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nuclei_norm = exposure.rescale_intensity(nuclei, in_range=(np.percentile(nuclei, 1), np.percentile(nuclei, 99)))
    gray = (np.clip(nuclei_norm, 0, 1) * 255).astype(np.uint8)
    rgb = np.dstack([gray, gray, gray])
    edge = morphology.binary_dilation(nuclear_mask) ^ morphology.binary_erosion(nuclear_mask)
    rgb[edge] = [255, 255, 0]
    Image.fromarray(rgb).save(out_path)


def analyze(args: argparse.Namespace) -> pd.DataFrame:
    experiment_dir = Path(args.experiment_dir)
    raw_dir = experiment_dir / "raw"
    output_dir = Path(args.output_dir)
    per_tiff_dir = output_dir / "per_tiff_layer_csv"
    qc_dir = output_dir / "qc_overlays"
    per_tiff_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    pairs = pair_channel_files(raw_dir)
    pairs = {
        series: files
        for series, files in pairs.items()
        if (args.start_series is None or series >= args.start_series)
        and (args.end_series is None or series <= args.end_series)
    }
    if not pairs:
        raise RuntimeError(f"No paired C1/C2 TIFF series found in {raw_dir}")

    output_rows: list[dict[str, float | int | str]] = []

    for series, files in pairs.items():
        c1_shape = get_stack_shape(files["C1"])
        c2_shape = get_stack_shape(files["C2"])
        if c1_shape != c2_shape:
            raise ValueError(f"Shape mismatch in series {series}: C1 {c1_shape}, C2 {c2_shape}")

        z_layers, _, _ = c1_shape
        segmentation_z_index = z_layers // 2
        output_stem = output_stem_from_tiff(files["C1"])

        requested_rows: list[dict[str, float | int]] = []
        with Image.open(files["C1"]) as bcat_tif, Image.open(files["C2"]) as nuclei_tif:
            nuclei_mid = read_plane(nuclei_tif, segmentation_z_index)
            nuclear_mask = segment_nuclear_region(
                nuclei_mid,
                min_area=args.min_area,
                gaussian_sigma=args.gaussian_sigma,
                local_block_size=args.local_block_size,
                adaptive_threshold=args.adaptive_threshold,
            )
            nuclear_area_px = int(nuclear_mask.sum())
            save_qc_overlay(qc_dir / f"{output_stem}_qc.png", nuclei_mid, nuclear_mask)

            for z_index in range(z_layers):
                bcat = read_plane(bcat_tif, z_index)
                if bcat.shape != nuclear_mask.shape:
                    raise ValueError(
                        f"Plane shape mismatch in series {series}, z {z_index}: "
                        f"C1 {bcat.shape}, mask {nuclear_mask.shape}"
                    )

                nuclear_pixels = bcat[nuclear_mask]
                intensity = float(nuclear_pixels.sum()) if nuclear_pixels.size else 0.0
                requested_rows.append(
                    {
                        "stack_number": z_index + 1,
                        "intensity": intensity,
                        "nuclei_number": nuclear_area_px,
                        "target_value": intensity / nuclear_area_px if nuclear_area_px else np.nan,
                    }
                )

        requested_df = pd.DataFrame(requested_rows)
        csv_name = f"{output_stem}_stack_intensity.csv"
        requested_df.to_csv(per_tiff_dir / csv_name, index=False)
        output_rows.append(
            {
                "series": series,
                "csv_file": csv_name,
                "rows": len(requested_df),
                "nuclei_number": nuclear_area_px,
                "segmentation_stack_number": segmentation_z_index + 1,
            }
        )
        print(
            f"series {series:02d}: {nuclear_area_px} nuclear pixels, "
            f"{z_layers} z layers, segmentation z {segmentation_z_index + 1}",
            flush=True,
        )

    return pd.DataFrame(output_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Quantify C1 beta-catenin intensity inside the C2 binary nuclear region.")
    parser.add_argument("--experiment-dir", default=str(DEFAULT_EXPERIMENT_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--min-area", type=int, default=80, help="Minimum nucleus area in pixels.")
    parser.add_argument("--gaussian-sigma", type=float, default=1.4)
    parser.add_argument("--local-block-size", type=int, default=151)
    parser.add_argument("--adaptive-threshold", action="store_true", help="Use slower local thresholding for uneven fields.")
    parser.add_argument("--start-series", type=int)
    parser.add_argument("--end-series", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    analyze(parse_args())
