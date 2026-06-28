"""Generate all publication-quality figures."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

from trajectory.constants import MALIGNANCY_COLORS, TISSUE_COLORS
from trajectory.io import load_adata
from trajectory.style import apply_style

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

FIG_DIR = Path(__file__).resolve().parents[1] / "results" / "figures"


def save_fig(fig: plt.Figure, name: str, dpi: int = 200) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logger.info("Saved %s", path)


def plot_qc_overview(adata_name: str = "adata_qc") -> None:
    """QC violin plots: genes, counts, mito %."""
    apply_style()
    adata = load_adata(adata_name)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.violin(adata, "n_genes_by_counts", groupby="tissue_type",
                 ax=axes[0], show=False, rotation=45)
    sc.pl.violin(adata, "total_counts", groupby="tissue_type",
                 ax=axes[1], show=False, rotation=45)
    sc.pl.violin(adata, "pct_counts_mt", groupby="tissue_type",
                 ax=axes[2], show=False, rotation=45)
    fig.suptitle("QC metrics by tissue type")
    save_fig(fig, "qc_violins")


def plot_umap_tissue(adata_name: str = "adata_qc") -> None:
    """UMAP colored by tissue type."""
    apply_style(scatter=True)
    adata = load_adata(adata_name)

    fig, ax = plt.subplots(figsize=(10, 8))
    for tt, color in TISSUE_COLORS.items():
        mask = adata.obs["tissue_type"] == tt
        if mask.sum() == 0:
            continue
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=color, label=tt, s=1, alpha=0.5, rasterized=True,
        )
    ax.legend(markerscale=10, frameon=False)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("All cells by tissue type")
    save_fig(fig, "umap_tissue_type")


def plot_umap_pseudotime(adata_name: str = "adata_trajectory") -> None:
    """UMAP colored by Palantir pseudotime."""
    apply_style(scatter=True)
    adata = load_adata(adata_name)

    fig, axes = plt.subplots(1, 3, figsize=(24, 7))

    sc.pl.umap(adata, color="palantir_pseudotime", ax=axes[0],
               show=False, title="Palantir pseudotime", frameon=False)

    for tt, color in TISSUE_COLORS.items():
        mask = adata.obs["tissue_type"] == tt
        if mask.sum() == 0:
            continue
        axes[1].scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=color, label=tt, s=2, alpha=0.5, rasterized=True,
        )
    axes[1].legend(markerscale=8, frameon=False)
    axes[1].set_title("Tissue type")

    for mal, color in MALIGNANCY_COLORS.items():
        mask = adata.obs["malignancy"] == mal
        if mask.sum() == 0:
            continue
        axes[2].scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=color, label=mal, s=2, alpha=0.5, rasterized=True,
        )
    axes[2].legend(markerscale=8, frameon=False)
    axes[2].set_title("Malignancy status")

    save_fig(fig, "umap_pseudotime_panel")


def plot_directionality_comparison(adata_name: str = "adata_trajectory") -> None:
    """Scatter: Palantir pseudotime vs CytoTRACE pseudotime."""
    apply_style()
    adata = load_adata(adata_name)

    if "ct_pseudotime" not in adata.obs:
        logger.warning("CytoTRACE not run yet, skipping")
        return

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        adata.obs["palantir_pseudotime"],
        adata.obs["ct_pseudotime"],
        c=adata.obs["palantir_pseudotime"],
        cmap="viridis", s=1, alpha=0.3, rasterized=True,
    )
    ax.set_xlabel("Palantir pseudotime")
    ax.set_ylabel("CytoTRACE pseudotime")
    ax.set_title("Directionality agreement")

    from scipy.stats import spearmanr
    rho, _ = spearmanr(adata.obs["palantir_pseudotime"], adata.obs["ct_pseudotime"])
    ax.text(0.05, 0.95, f"Spearman ρ = {rho:.3f}",
            transform=ax.transAxes, va="top", fontsize=18)

    save_fig(fig, "directionality_comparison")


def plot_pseudotime_heatmap(adata_name: str = "adata_trajectory") -> None:
    """Gene expression heatmap ordered by pseudotime."""
    apply_style()
    adata = load_adata(adata_name)

    import pandas as pd
    tables_dir = Path(__file__).resolve().parents[1] / "results" / "tables"
    de_path = tables_dir / "de_early_vs_late.csv"
    if not de_path.exists():
        logger.warning("DE results not found, skipping heatmap")
        return

    de = pd.read_csv(de_path)
    top_up = de.nlargest(15, "logfoldchanges")["names"].tolist()
    top_down = de.nsmallest(15, "logfoldchanges")["names"].tolist()
    genes = top_up + top_down
    genes = [g for g in genes if g in adata.var_names]

    if len(genes) < 4:
        logger.warning("Too few DE genes for heatmap")
        return

    ordered = adata[adata.obs.sort_values("palantir_pseudotime").index]

    fig = sc.pl.heatmap(
        ordered, var_names=genes, groupby="pt_bin",
        swap_axes=True, show=False, figsize=(12, 8),
    )
    save_fig(fig["heatmap_ax"].figure, "pseudotime_heatmap")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    plot_qc_overview()
    plot_umap_tissue()
    plot_umap_pseudotime()
    plot_directionality_comparison()
    plot_pseudotime_heatmap()

    print(f"Figures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
