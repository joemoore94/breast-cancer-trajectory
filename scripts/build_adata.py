"""Load per-sample MTX files from GEO and build a combined AnnData."""

import logging

from trajectory.io import load_samples, save_adata

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> None:
    adata = load_samples()
    print(f"Combined: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"Tissue types:\n{adata.obs['tissue_type'].value_counts()}")
    save_adata(adata, "adata_raw")


if __name__ == "__main__":
    main()
