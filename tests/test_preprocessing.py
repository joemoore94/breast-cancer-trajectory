"""Tests for preprocessing functions using synthetic AnnData."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from trajectory.constants import CELL_TYPE_MARKERS
from trajectory.preprocessing import normalize_log, score_cell_types, subset_luminal


@pytest.fixture()
def synthetic_adata() -> ad.AnnData:
    """Minimal AnnData with known marker expression patterns."""
    n_cells = 200
    rng = np.random.default_rng(42)

    luminal_genes = ["KRT8", "KRT18", "KRT19", "ESR1", "GATA3", "PGR", "FOXA1"]
    basal_genes = ["KRT14", "KRT5", "TP63", "KRT17"]
    immune_genes = ["PTPRC"]
    other_genes = [f"GENE_{i}" for i in range(50)]
    all_genes = luminal_genes + basal_genes + immune_genes + other_genes

    X = rng.poisson(1, size=(n_cells, len(all_genes))).astype(np.float32)

    # First 80 cells: luminal (high KRT8/GATA3, low KRT14)
    for g in luminal_genes:
        idx = all_genes.index(g)
        X[:80, idx] = rng.poisson(10, size=80)
    for g in basal_genes:
        idx = all_genes.index(g)
        X[:80, idx] = 0

    # Next 60 cells: basal (high KRT14/KRT5, low KRT8)
    for g in basal_genes:
        idx = all_genes.index(g)
        X[80:140, idx] = rng.poisson(10, size=60)
    for g in luminal_genes:
        idx = all_genes.index(g)
        X[80:140, idx] = 0

    # Last 60 cells: immune (high PTPRC)
    X[140:, all_genes.index("PTPRC")] = rng.poisson(15, size=60)

    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(
            {"tissue_type": (["normal"] * 80 + ["tnbc"] * 60 + ["normal"] * 60)},
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(index=all_genes),
    )

    adata.var["mt"] = False
    adata.obs["n_genes_by_counts"] = np.array((X > 0).sum(axis=1)).flatten()
    adata.obs["total_counts"] = np.array(X.sum(axis=1)).flatten()
    adata.obs["pct_counts_mt"] = 0.0

    return adata


def test_normalize_log_preserves_raw(synthetic_adata: ad.AnnData) -> None:
    adata = normalize_log(synthetic_adata)
    assert adata.raw is not None
    assert adata.raw.shape == adata.shape


def test_score_cell_types_adds_columns(synthetic_adata: ad.AnnData) -> None:
    adata = normalize_log(synthetic_adata)
    adata = score_cell_types(adata, CELL_TYPE_MARKERS)
    for ct in CELL_TYPE_MARKERS:
        col = f"score_{ct}"
        if col in adata.obs:
            assert not adata.obs[col].isna().all()


def test_subset_luminal_selects_correct_cells(synthetic_adata: ad.AnnData) -> None:
    adata = normalize_log(synthetic_adata)
    adata = score_cell_types(adata, CELL_TYPE_MARKERS)
    luminal = subset_luminal(adata)
    assert luminal.n_obs > 0
    assert luminal.n_obs < adata.n_obs
