"""Palantir diffusion map and pseudotime computation."""

import logging

from trajectory.io import load_adata, save_adata
from trajectory.pseudotime import find_root_cell, run_palantir

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> None:
    adata = load_adata("adata_luminal_cnv")
    print(f"Loaded: {adata.n_obs} cells")

    root = find_root_cell(adata)
    print(f"Root cell: {root}")

    adata = run_palantir(adata, root)

    print(f"Pseudotime range: [{adata.obs['palantir_pseudotime'].min():.3f}, "
          f"{adata.obs['palantir_pseudotime'].max():.3f}]")

    save_adata(adata, "adata_trajectory")


if __name__ == "__main__":
    main()
