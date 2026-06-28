"""QC, normalization, HVG selection, and batch correction."""

from __future__ import annotations

import logging

import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)


def qc_filter(
    adata: ad.AnnData,
    *,
    min_genes: int = 500,
    max_genes: int = 8000,
    max_pct_mito: float = 20.0,
    min_cells_per_gene: int = 10,
) -> ad.AnnData:
    """Apply standard scRNA-seq QC filters.

    Computes mitochondrial QC metrics, flags doublets with Scrublet, and
    filters cells and genes by thresholds.
    """
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    import scrublet as scr

    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    n_before = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mito)
        & (~adata.obs["predicted_doublet"])
    ].copy()

    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    logger.info(
        "QC: %d -> %d cells, %d genes",
        n_before, adata.n_obs, adata.n_vars,
    )
    return adata


def normalize_log(adata: ad.AnnData) -> ad.AnnData:
    """Normalize to 10k counts per cell and log-transform.

    Stores raw counts in ``adata.raw`` before normalization.
    """
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def select_hvg(
    adata: ad.AnnData,
    *,
    n_top_genes: int = 3000,
    batch_key: str | None = "sample_id",
) -> ad.AnnData:
    """Select highly variable genes, optionally batch-aware."""
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        flavor="seurat_v3",
        subset=False,
    )
    logger.info("HVGs: %d selected", adata.var["highly_variable"].sum())
    return adata


def run_pca(adata: ad.AnnData, n_comps: int = 50) -> ad.AnnData:
    """PCA on highly variable genes."""
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)
    return adata


def batch_correct_harmony(
    adata: ad.AnnData,
    batch_key: str = "patient_id",
) -> ad.AnnData:
    """Harmony integration on PCA embeddings.

    Stores corrected embeddings in ``adata.obsm['X_pca_harmony']`` and builds
    the neighbor graph from the corrected space.
    """
    import harmonypy as hm

    harmony = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch_key,
    )
    adata.obsm["X_pca_harmony"] = harmony.Z_corr.T
    logger.info("Harmony integration complete on '%s'", batch_key)
    return adata


def build_neighbors_umap(
    adata: ad.AnnData,
    *,
    use_rep: str = "X_pca_harmony",
    n_neighbors: int = 30,
) -> ad.AnnData:
    """Build kNN graph and compute UMAP from a specified representation."""
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    return adata


def score_cell_types(
    adata: ad.AnnData,
    marker_dict: dict[str, list[str]],
) -> ad.AnnData:
    """Score cells for each cell-type marker panel using sc.tl.score_genes."""
    for ct_name, markers in marker_dict.items():
        present = [g for g in markers if g in adata.var_names]
        if len(present) < 2:
            logger.warning(
                "Skipping %s: only %d/%d markers present",
                ct_name, len(present), len(markers),
            )
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ct_name}")
    return adata


def subset_luminal(adata: ad.AnnData, *, threshold: float = 0.0) -> ad.AnnData:
    """Subset to luminal epithelial cells based on marker scores.

    Keeps cells with positive luminal progenitor or mature luminal scores
    and negative basal scores.
    """
    lp = adata.obs.get("score_luminal_progenitor", 0)
    ml = adata.obs.get("score_mature_luminal", 0)
    basal = adata.obs.get("score_basal", 0)

    luminal_mask = (
        ((lp > threshold) | (ml > threshold))
        & (basal <= threshold)
    )
    sub = adata[luminal_mask].copy()
    logger.info("Luminal subset: %d / %d cells", sub.n_obs, adata.n_obs)
    return sub
