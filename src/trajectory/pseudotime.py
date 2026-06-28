"""Pseudotime, diffusion maps, and directionality analysis."""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import pandas as pd

logger = logging.getLogger(__name__)


def find_root_cell(
    adata: ad.AnnData,
    *,
    progenitor_score: str = "score_luminal_progenitor",
    malignancy_col: str = "malignancy",
    normal_label: str = "normal",
) -> str:
    """Find the best root cell: highest progenitor score among normal cells."""
    normal_mask = adata.obs[malignancy_col] == normal_label
    if normal_mask.sum() == 0:
        raise ValueError("No normal cells found for root selection")

    normal_scores = adata.obs.loc[normal_mask, progenitor_score]
    root_idx = normal_scores.idxmax()
    logger.info(
        "Root cell: %s (progenitor score=%.3f)",
        root_idx, normal_scores[root_idx],
    )
    return root_idx


def run_palantir(
    adata: ad.AnnData,
    root_cell: str,
    *,
    n_components: int = 10,
    knn: int = 30,
    use_rep: str = "X_pca_harmony",
    num_waypoints: int = 1200,
) -> ad.AnnData:
    """Run Palantir diffusion map and pseudotime computation.

    Uses the AnnData-native API (palantir >= 1.3). Results are stored
    directly in ``adata.obs`` (``palantir_pseudotime``,
    ``palantir_entropy``) and ``adata.obsm``.
    """
    import palantir

    palantir.utils.run_diffusion_maps(
        adata,
        n_components=n_components,
        knn=knn,
        pca_key=use_rep,
    )
    palantir.utils.determine_multiscale_space(adata)

    palantir.core.run_palantir(
        adata,
        root_cell,
        num_waypoints=num_waypoints,
    )

    logger.info(
        "Palantir pseudotime range: [%.3f, %.3f]",
        adata.obs["palantir_pseudotime"].min(),
        adata.obs["palantir_pseudotime"].max(),
    )
    return adata


def run_cellrank_pseudotime(
    adata: ad.AnnData,
    *,
    time_key: str = "palantir_pseudotime",
) -> Any:
    """Run CellRank 2 PseudotimeKernel using Palantir pseudotime as prior.

    Returns the CellRank estimator for further analysis.
    """
    import cellrank as cr

    pk = cr.kernels.PseudotimeKernel(adata, time_key=time_key)
    pk.compute_transition_matrix(threshold_scheme="soft")

    estimator = cr.estimators.GPCCA(pk)
    estimator.compute_schur(n_components=20)
    estimator.compute_macrostates(n_states=6, cluster_key="leiden")

    logger.info("CellRank PseudotimeKernel macrostates computed")
    return estimator


def run_cytotrace(adata: ad.AnnData) -> ad.AnnData:
    """Run CytoTRACE via CellRank's CytoTRACEKernel.

    ``compute_cytotrace`` writes directly to ``adata.obs``:
    ``ct_score``, ``ct_pseudotime``, ``ct_num_exp_genes``.
    """
    import cellrank as cr

    ctk = cr.kernels.CytoTRACEKernel(adata)
    ctk.compute_cytotrace(layer="X", use_raw=False)

    logger.info(
        "CytoTRACE score range: [%.3f, %.3f]",
        adata.obs["ct_score"].min(), adata.obs["ct_score"].max(),
    )
    return adata


def compare_directionality(adata: ad.AnnData) -> pd.DataFrame:
    """Compute pairwise Spearman correlations between directionality estimates.

    Compares Palantir pseudotime, CytoTRACE pseudotime, and (if available)
    CellRank absorption probabilities.
    """
    from scipy.stats import spearmanr

    cols = []
    if "palantir_pseudotime" in adata.obs:
        cols.append("palantir_pseudotime")
    if "ct_pseudotime" in adata.obs:
        cols.append("ct_pseudotime")

    if len(cols) < 2:
        logger.warning("Need at least 2 pseudotime columns to compare")
        return pd.DataFrame()

    results = []
    for i, c1 in enumerate(cols):
        for c2 in cols[i + 1:]:
            rho, pval = spearmanr(adata.obs[c1], adata.obs[c2])
            results.append({"method_1": c1, "method_2": c2, "spearman_rho": rho, "p_value": pval})

    df = pd.DataFrame(results)
    logger.info("Directionality comparison:\n%s", df.to_string(index=False))
    return df
