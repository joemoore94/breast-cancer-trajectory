"""QC filtering, normalization, HVG selection, and batch correction."""

import logging

from trajectory.io import load_adata, save_adata
from trajectory.preprocessing import (
    batch_correct_harmony,
    build_neighbors_umap,
    normalize_log,
    qc_filter,
    run_pca,
    select_hvg,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> None:
    adata = load_adata("adata_raw")
    print(f"Loaded: {adata.n_obs} cells")

    adata = qc_filter(adata)
    adata = normalize_log(adata)
    adata = select_hvg(adata)
    adata = run_pca(adata)
    adata = batch_correct_harmony(adata)
    adata = build_neighbors_umap(adata)

    print(f"After QC + integration: {adata.n_obs} cells, {adata.n_vars} genes")
    save_adata(adata, "adata_qc")


if __name__ == "__main__":
    main()
