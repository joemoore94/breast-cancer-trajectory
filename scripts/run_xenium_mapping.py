"""Map trajectory-associated genes onto the 380-gene Xenium breast panel."""

import logging
from pathlib import Path

import pandas as pd

from trajectory.constants import XENIUM_TRAJECTORY_GENES

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

TABLES_DIR = Path(__file__).resolve().parents[1] / "results" / "tables"

SEGBENCH_ROOT = Path(__file__).resolve().parents[2] / "segmentation-benchmark"


def load_xenium_panel() -> list[str]:
    """Load the full 380-gene Xenium panel from segmentation-benchmark.

    Falls back to the curated subset in constants.py if the segmentation-
    benchmark repo is not available.
    """
    panel_file = SEGBENCH_ROOT / "data" / "raw" / "panel.tsv"
    if panel_file.exists():
        df = pd.read_csv(panel_file, sep="\t")
        col = "Name" if "Name" in df.columns else df.columns[0]
        genes = df[col].tolist()
        logger.info("Loaded %d genes from Xenium panel", len(genes))
        return genes

    adata_path = SEGBENCH_ROOT / "data" / "processed" / "roi" / "adata_10x.h5ad"
    if adata_path.exists():
        import scanpy as sc
        adata = sc.read_h5ad(adata_path, backed="r")
        genes = adata.var_names.tolist()
        logger.info("Loaded %d genes from 10x AnnData", len(genes))
        return genes

    logger.warning(
        "segmentation-benchmark not found; using curated subset (%d genes)",
        len(XENIUM_TRAJECTORY_GENES),
    )
    return XENIUM_TRAJECTORY_GENES


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    de_path = TABLES_DIR / "de_early_vs_late.csv"
    if not de_path.exists():
        print(f"DE results not found at {de_path}. Run run_de_programs.py first.")
        return

    de = pd.read_csv(de_path)
    sig = de[
        (de["pvals_adj"] < 0.05)
        & (de["logfoldchanges"].abs() > 0.5)
    ].copy()

    xenium_panel = load_xenium_panel()
    xenium_set = set(xenium_panel)

    sig["in_xenium_panel"] = sig["names"].isin(xenium_set)
    overlap = sig[sig["in_xenium_panel"]]

    print(f"Significant trajectory DE genes: {len(sig)}")
    print(f"  Present in Xenium panel: {len(overlap)}")
    print(f"  Not in panel: {len(sig) - len(overlap)}")

    if not overlap.empty:
        print(f"\nXenium-measurable trajectory genes:")
        for _, row in overlap.sort_values("logfoldchanges", ascending=False).iterrows():
            direction = "UP" if row["logfoldchanges"] > 0 else "DOWN"
            print(f"  {row['names']:12s}  logFC={row['logfoldchanges']:+.2f}  {direction}")

    sig.to_csv(TABLES_DIR / "xenium_panel_trajectory_genes.csv", index=False)
    logger.info("Saved xenium panel mapping to %s", TABLES_DIR)


if __name__ == "__main__":
    main()
