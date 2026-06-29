"""CellRank 2 and CytoTRACE directionality comparison."""

import logging
from pathlib import Path

from trajectory.io import load_adata, save_adata
from trajectory.pseudotime import (
    compare_directionality,
    run_cellrank_pseudotime,
    run_cytotrace,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

TABLES_DIR = Path(__file__).resolve().parents[1] / "results" / "tables"


def main() -> None:
    adata = load_adata("adata_trajectory")
    print(f"Loaded: {adata.n_obs} cells")

    adata = run_cytotrace(adata)
    estimator = run_cellrank_pseudotime(adata)

    df = compare_directionality(adata)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    if not df.empty:
        print(df.to_string(index=False))
    df.to_csv(TABLES_DIR / "directionality_comparison.csv", index=False)

    save_adata(adata, "adata_directionality")


if __name__ == "__main__":
    main()
