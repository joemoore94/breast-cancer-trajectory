"""Differential expression along pseudotime and gene program identification."""

import logging
from pathlib import Path

import pandas as pd
import scanpy as sc

from trajectory.io import load_adata, save_adata

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

TABLES_DIR = Path(__file__).resolve().parents[1] / "results" / "tables"


def de_along_pseudotime(
    adata: sc.AnnData,
    *,
    n_bins: int = 5,
    pseudotime_col: str = "palantir_pseudotime",
) -> pd.DataFrame:
    """Bin cells by pseudotime and run pairwise Wilcoxon DE."""
    adata.obs["pt_bin"] = pd.cut(
        adata.obs[pseudotime_col],
        bins=n_bins,
        labels=[f"bin_{i}" for i in range(n_bins)],
    )

    sc.tl.rank_genes_groups(
        adata,
        groupby="pt_bin",
        method="wilcoxon",
        use_raw=True,
    )

    results = []
    for group in adata.obs["pt_bin"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=group)
        df["pt_bin"] = group
        results.append(df)

    return pd.concat(results, ignore_index=True)


def early_vs_late_de(
    adata: sc.AnnData,
    *,
    pseudotime_col: str = "palantir_pseudotime",
    quantile: float = 0.25,
) -> pd.DataFrame:
    """DE between earliest and latest pseudotime quartiles."""
    pt = adata.obs[pseudotime_col]
    early_thresh = pt.quantile(quantile)
    late_thresh = pt.quantile(1 - quantile)

    adata.obs["trajectory_phase"] = "middle"
    adata.obs.loc[pt <= early_thresh, "trajectory_phase"] = "early"
    adata.obs.loc[pt >= late_thresh, "trajectory_phase"] = "late"

    sub = adata[adata.obs["trajectory_phase"].isin(["early", "late"])].copy()

    sc.tl.rank_genes_groups(
        sub,
        groupby="trajectory_phase",
        groups=["late"],
        reference="early",
        method="wilcoxon",
        use_raw=True,
    )

    return sc.get.rank_genes_groups_df(sub, group="late")


def main() -> None:
    adata = load_adata("adata_trajectory")
    print(f"Loaded: {adata.n_obs} cells")

    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    df_bins = de_along_pseudotime(adata)
    df_bins.to_csv(TABLES_DIR / "de_pseudotime_bins.csv", index=False)
    logger.info("Pseudotime-binned DE: %d results", len(df_bins))

    df_early_late = early_vs_late_de(adata)
    df_early_late.to_csv(TABLES_DIR / "de_early_vs_late.csv", index=False)
    logger.info("Early vs late DE: %d genes", len(df_early_late))

    sig = df_early_late[
        (df_early_late["pvals_adj"] < 0.05)
        & (df_early_late["logfoldchanges"].abs() > 1.0)
    ]
    print(f"Significant DE genes (early vs late): {len(sig)}")
    print(f"  Upregulated in late: {(sig['logfoldchanges'] > 0).sum()}")
    print(f"  Downregulated in late: {(sig['logfoldchanges'] < 0).sum()}")

    save_adata(adata, "adata_trajectory")


if __name__ == "__main__":
    main()
