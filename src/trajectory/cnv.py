"""CNV inference wrappers using infercnvpy."""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
import infercnvpy as cnv

from trajectory.constants import NORMAL_TYPES

logger = logging.getLogger(__name__)


def infer_cnv(
    adata: ad.AnnData,
    *,
    reference_key: str = "tissue_type",
    reference_cat: set[str] | None = None,
    window_size: int = 100,
) -> ad.AnnData:
    """Run infercnvpy CNV inference using normal cells as reference.

    Parameters
    ----------
    adata:
        AnnData with raw counts in ``.raw`` and ``reference_key`` in ``.obs``.
    reference_key:
        Column in ``.obs`` that identifies tissue type.
    reference_cat:
        Values of ``reference_key`` to use as the normal reference.
        Defaults to :data:`NORMAL_TYPES`.
    window_size:
        Genomic window size for smoothing.
    """
    reference_cat = reference_cat or NORMAL_TYPES

    cnv.tl.infercnv(
        adata,
        reference_key=reference_key,
        reference_cat=list(reference_cat),
        window_size=window_size,
    )
    cnv.tl.pca(adata)
    cnv.tl.leiden(adata)
    cnv.tl.cnv_score(adata)

    logger.info(
        "CNV inference complete. Score range: [%.3f, %.3f]",
        adata.obs["cnv_score"].min(),
        adata.obs["cnv_score"].max(),
    )
    return adata


def classify_malignancy(
    adata: ad.AnnData,
    *,
    score_col: str = "cnv_score",
    tissue_col: str = "tissue_type",
    threshold: float | None = None,
) -> ad.AnnData:
    """Label cells as normal, preneoplastic, or malignant.

    Uses a combination of CNV score and known tissue type. Cells from normal
    samples with low CNV scores are labeled normal; cells from BRCA1 samples
    are labeled preneoplastic; cells with high CNV scores from tumor samples
    are labeled malignant.

    If ``threshold`` is None, uses Otsu's method on the CNV score distribution
    of tumor-derived cells.
    """
    labels = np.full(adata.n_obs, "unknown", dtype=object)

    normal_mask = adata.obs[tissue_col].isin(NORMAL_TYPES).values
    brca1_mask = (adata.obs[tissue_col] == "brca1").values
    tumor_mask = ~normal_mask & ~brca1_mask

    if threshold is None:
        tumor_scores = adata.obs.loc[tumor_mask, score_col].values
        from skimage.filters import threshold_otsu

        threshold = float(threshold_otsu(tumor_scores))
        logger.info("Otsu threshold on tumor CNV scores: %.4f", threshold)

    scores = adata.obs[score_col].values

    labels[normal_mask & (scores <= threshold)] = "normal"
    labels[normal_mask & (scores > threshold)] = "normal_high_cnv"
    labels[brca1_mask] = "preneoplastic"
    labels[tumor_mask & (scores > threshold)] = "malignant"
    labels[tumor_mask & (scores <= threshold)] = "tumor_low_cnv"

    adata.obs["malignancy"] = pd.Categorical(
        labels, categories=["normal", "normal_high_cnv", "preneoplastic",
                            "tumor_low_cnv", "malignant", "unknown"]
    )
    logger.info("Malignancy labels:\n%s", adata.obs["malignancy"].value_counts())
    return adata
