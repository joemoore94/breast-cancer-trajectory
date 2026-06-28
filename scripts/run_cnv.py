"""Run CNV inference on luminal epithelial cells and classify malignancy."""

import logging

from trajectory.cnv import classify_malignancy, infer_cnv
from trajectory.io import load_adata, save_adata

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> None:
    adata = load_adata("adata_luminal")
    print(f"Loaded: {adata.n_obs} luminal cells")

    adata = infer_cnv(adata)
    adata = classify_malignancy(adata)

    print(f"Malignancy labels:\n{adata.obs['malignancy'].value_counts()}")
    save_adata(adata, "adata_luminal_cnv")


if __name__ == "__main__":
    main()
