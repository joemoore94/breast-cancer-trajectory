"""Subset to luminal epithelial cells using marker-based scoring."""

import logging

import scanpy as sc

from trajectory.constants import CELL_TYPE_MARKERS
from trajectory.io import load_adata, save_adata
from trajectory.preprocessing import (
    batch_correct_harmony,
    build_neighbors_umap,
    run_pca,
    score_cell_types,
    select_hvg,
    subset_luminal,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> None:
    adata = load_adata("adata_qc")
    print(f"Loaded: {adata.n_obs} cells")

    adata = score_cell_types(adata, CELL_TYPE_MARKERS)

    sc.tl.leiden(adata, resolution=1.0)

    adata_lum = subset_luminal(adata)
    print(f"Luminal subset: {adata_lum.n_obs} cells")
    print(f"Tissue types:\n{adata_lum.obs['tissue_type'].value_counts()}")

    # Re-run HVG/PCA/Harmony/UMAP on the luminal subset
    adata_lum = select_hvg(adata_lum)
    adata_lum = run_pca(adata_lum)
    adata_lum = batch_correct_harmony(adata_lum)
    adata_lum = build_neighbors_umap(adata_lum)

    save_adata(adata_lum, "adata_luminal")


if __name__ == "__main__":
    main()
