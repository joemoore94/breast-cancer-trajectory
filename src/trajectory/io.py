"""Data download and loading utilities for GSE161529."""

from __future__ import annotations

import logging
import tarfile
from pathlib import Path

import anndata as ad
import scanpy as sc

from trajectory.constants import SAMPLE_METADATA

logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

GEO_ACCESSION = "GSE161529"
GEO_TAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/suppl/"
    "GSE161529_RAW.tar"
)
GEO_FEATURES_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/suppl/"
    "GSE161529_features.tsv.gz"
)


def download_geo(dest: Path | None = None, force: bool = False) -> Path:
    """Download and extract GSE161529_RAW.tar from GEO.

    Returns the directory containing extracted per-sample folders.
    """
    import urllib.request

    dest = dest or RAW_DIR / GEO_ACCESSION
    dest.mkdir(parents=True, exist_ok=True)

    tar_path = dest / "GSE161529_RAW.tar"
    if tar_path.exists() and not force:
        logger.info("TAR already downloaded: %s", tar_path)
    else:
        logger.info("Downloading %s ...", GEO_TAR_URL)
        urllib.request.urlretrieve(GEO_TAR_URL, tar_path)
        logger.info("Download complete: %s", tar_path)

    extracted_marker = dest / ".extracted"
    if extracted_marker.exists() and not force:
        logger.info("TAR already extracted")
    else:
        logger.info("Extracting %s ...", tar_path)
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(path=dest, filter="data")
        extracted_marker.touch()
        logger.info("Extraction complete")

    features_path = dest / "GSE161529_features.tsv.gz"
    if not features_path.exists() or force:
        logger.info("Downloading features file ...")
        urllib.request.urlretrieve(GEO_FEATURES_URL, features_path)

    return dest


def _prepare_sample_dirs(geo_dir: Path) -> dict[str, Path]:
    """Organize flat GEO files into per-sample directories for sc.read_10x_mtx.

    GSE161529_RAW.tar extracts to flat files like:
        GSM4909253_N-PM0092-Total-barcodes.tsv.gz
        GSM4909253_N-PM0092-Total-matrix.mtx.gz
    with a single shared GSE161529_features.tsv.gz.

    This function creates per-sample subdirectories with standard 10x names
    (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz) so that
    sc.read_10x_mtx() can read them directly.
    """
    features_src = geo_dir / "GSE161529_features.tsv.gz"
    if not features_src.exists():
        raise FileNotFoundError(
            f"Shared features file not found: {features_src}. "
            "Run download_geo.py first."
        )

    matrix_files = sorted(geo_dir.glob("GSM*matrix.mtx.gz"))
    if not matrix_files:
        matrix_files = sorted(geo_dir.glob("GSM*matrix.mtx"))

    sample_dirs: dict[str, Path] = {}
    for mtx_path in matrix_files:
        gsm_id = mtx_path.name.split("_")[0]

        barcode_candidates = list(geo_dir.glob(f"{gsm_id}_*barcodes.tsv.gz"))
        if not barcode_candidates:
            barcode_candidates = list(geo_dir.glob(f"{gsm_id}_*barcodes.tsv"))
        if not barcode_candidates:
            logger.warning("No barcodes file for %s, skipping", gsm_id)
            continue

        sample_dir = geo_dir / gsm_id
        sample_dir.mkdir(exist_ok=True)

        mtx_link = sample_dir / "matrix.mtx.gz"
        bc_link = sample_dir / "barcodes.tsv.gz"
        feat_link = sample_dir / "features.tsv.gz"

        for link, target in [
            (mtx_link, mtx_path),
            (bc_link, barcode_candidates[0]),
            (feat_link, features_src),
        ]:
            if not link.exists():
                link.symlink_to(target.resolve())

        sample_dirs[gsm_id] = sample_dir

    logger.info("Prepared %d sample directories", len(sample_dirs))
    return sample_dirs


def load_samples(geo_dir: Path | None = None) -> ad.AnnData:
    """Load all per-sample MTX files and concatenate into a single AnnData.

    Adds ``sample_id``, ``sample_name``, ``tissue_type``, ``patient_id``, and
    ``fraction`` columns to ``.obs`` from :data:`SAMPLE_METADATA`.
    """
    geo_dir = geo_dir or RAW_DIR / GEO_ACCESSION
    sample_dirs = _prepare_sample_dirs(geo_dir)
    logger.info("Found %d sample directories", len(sample_dirs))

    adatas: list[ad.AnnData] = []
    for gsm_id, sample_dir in sample_dirs.items():
        if gsm_id not in SAMPLE_METADATA:
            logger.warning("Skipping unknown sample: %s", gsm_id)
            continue

        sample_name, tissue_type, patient_id, fraction = SAMPLE_METADATA[gsm_id]
        logger.info("Loading %s (%s) ...", gsm_id, sample_name)

        try:
            adata = sc.read_10x_mtx(sample_dir, var_names="gene_symbols")
        except Exception:
            logger.error("Failed to load %s from %s", gsm_id, sample_dir)
            continue

        adata.obs["sample_id"] = gsm_id
        adata.obs["sample_name"] = sample_name
        adata.obs["tissue_type"] = tissue_type
        adata.obs["patient_id"] = patient_id
        adata.obs["fraction"] = fraction
        adata.obs_names_make_unique()
        adatas.append(adata)

    logger.info("Concatenating %d samples ...", len(adatas))
    # All samples share the same features.tsv.gz so gene sets should match.
    # Use join="inner" to keep only shared genes and avoid dense intermediate.
    combined = ad.concat(adatas, join="inner")
    del adatas
    combined.obs_names_make_unique()
    combined.var_names_make_unique()
    logger.info("Combined shape: %s", combined.shape)
    return combined


def save_adata(adata: ad.AnnData, name: str) -> Path:
    """Save an AnnData object to ``data/processed/{name}.h5ad``."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    path = PROCESSED_DIR / f"{name}.h5ad"
    adata.write_h5ad(path)
    logger.info("Saved %s (%s)", path, adata.shape)
    return path


def load_adata(name: str) -> ad.AnnData:
    """Load an AnnData object from ``data/processed/{name}.h5ad``."""
    path = PROCESSED_DIR / f"{name}.h5ad"
    return sc.read_h5ad(path)
