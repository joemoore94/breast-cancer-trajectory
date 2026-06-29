# Normal-to-malignant luminal epithelial trajectory in breast cancer

Can we reconstruct a normal-to-malignant luminal epithelial trajectory from a static scRNA-seq snapshot, and what gene expression programs define progression along that axis?

## Dataset

Pal et al. (2021), GEO accession [GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529). 421,761 cells from 55 patients on 10x Chromium, spanning normal breast tissue (reduction mammoplasties), BRCA1 preneoplastic tissue, and malignant tumors across ER+, HER2+, and TNBC subtypes. The normal and preneoplastic samples are the key advantage over tumor-only atlases - they anchor the healthy end of the trajectory and provide a known high-risk intermediate state (BRCA1 carriers) that should sit between normal and malignant on the inferred axis.

69 samples total: 24 normal (13 patients, some with sorted epithelial and total fractions), 4 BRCA1 preneoplastic, 4 TNBC, 4 TNBC-BRCA1, 6 HER2+, 20 ER+ (including tumor/lymph-node pairs), and 3 male ER+. Sample metadata extracted from the companion [GitHub repo](https://github.com/yunshun/HumanBreast10X) (Chen et al., *Scientific Data*, 2022).

## Approach

The analysis subsets to luminal epithelial cells (KRT8+/KRT14−/GATA3+), infers copy number variation to separate normal from malignant cells, then fits a pseudotime trajectory rooted at the luminal progenitor population from normal tissue.

**Preprocessing.** Per-sample count matrices from GEO are concatenated, QC-filtered (mitochondrial %, gene counts, Scrublet doublet detection), normalized (10k + log1p), and integrated across patients with Harmony on the top 3,000 HVGs.

**Cell type subsetting.** Cells are scored against curated marker panels (luminal progenitor, mature luminal, basal, immune, endothelial, fibroblast) using `sc.tl.score_genes`. Luminal epithelial cells are retained; basal, immune, and stromal populations are excluded. The luminal subset is re-embedded (HVG → PCA → Harmony → UMAP) to resolve finer structure.

**CNV inference.** infercnvpy infers large-scale chromosomal copy number from expression, using cells from normal samples as the reference baseline. An Otsu threshold on CNV scores from tumor-derived cells separates malignant from non-malignant cells. BRCA1-sample cells are labeled preneoplastic regardless of CNV score.

**Pseudotime.** Palantir computes a diffusion map and pseudotime from the Harmony-corrected PCA, rooted at the cell with the highest luminal progenitor score among normal cells (KRT15+, low CNV). This produces a continuous ordering from normal luminal progenitor → mature luminal → malignant.

**Directionality validation.** Two independent checks confirm the trajectory direction without RNA velocity (which requires raw FASTQs not available on GEO). CellRank 2's PseudotimeKernel uses Palantir pseudotime as a prior to compute transition probabilities. CytoTRACE, via CellRank's CytoTRACEKernel, estimates differentiation state from transcriptional diversity alone. Spearman correlation between all three orderings quantifies agreement.

**Gene programs.** Wilcoxon DE between early (bottom quartile) and late (top quartile) pseudotime identifies genes that activate or repress during progression. Pseudotime-binned DE across five bins captures non-monotonic expression dynamics.

**Xenium panel mapping.** Trajectory-associated DE genes are cross-referenced against the 380-gene Xenium breast panel from the companion segmentation-benchmark project, identifying which progression markers are spatially measurable and setting up future spatial trajectory analysis.

## Workflow

The preferred build path is Snakemake. It tracks each derived file, resumes
from the last completed step, and keeps long-running outputs tied to the script
that produced them.

```bash
conda env create -f environment.yml
conda activate trajectory
pip install -e .

# Preview the dependency graph without running anything
snakemake -n

# Run the full analysis
snakemake --cores 4
```

The workflow configuration lives in `configs/workflow.yaml`. Increase or
decrease thread counts there, then rerun the same `snakemake --cores N`
command. Logs for each step are written to `logs/`.

Useful targets:

```bash
# Stop after the integrated all-cell object
snakemake --cores 4 data/processed/adata_qc.h5ad

# Stop after Palantir pseudotime
snakemake --cores 4 data/processed/adata_trajectory.h5ad

# Rebuild only figures from existing processed data/tables
snakemake --cores 1 results/figures/umap_pseudotime_panel.png

# Draw a DAG if graphviz is installed
snakemake --dag | dot -Tpng > results/figures/snakemake_dag.png
```

The full workflow produces:

```text
data/processed/adata_raw.h5ad
data/processed/adata_qc.h5ad
data/processed/adata_luminal.h5ad
data/processed/adata_luminal_cnv.h5ad
data/processed/adata_trajectory.h5ad
data/processed/adata_directionality.h5ad
data/processed/adata_programs.h5ad
results/tables/*.csv
results/figures/*.png
```

## Script Steps

```
scripts/download_geo.py       # Download GSE161529_RAW.tar from GEO
scripts/build_adata.py        # Concatenate per-sample MTX → adata_raw.h5ad
scripts/run_qc.py             # QC, normalize, HVG, PCA, Harmony → adata_qc.h5ad
scripts/run_subset_luminal.py # Marker scoring, luminal subset → adata_luminal.h5ad
scripts/run_cnv.py            # infercnvpy, malignancy classification → adata_luminal_cnv.h5ad
scripts/run_pseudotime.py     # Palantir diffusion map + pseudotime → adata_trajectory.h5ad
scripts/run_directionality.py # CellRank + CytoTRACE comparison → adata_directionality.h5ad
scripts/run_de_programs.py    # DE along pseudotime axis → adata_programs.h5ad
scripts/run_xenium_mapping.py # Cross-reference with Xenium panel
scripts/make_figures.py       # Generate all figures
```

## Setup

Use this only if you want to run scripts manually instead of using Snakemake.

```bash
conda env create -f environment.yml
conda activate trajectory
pip install -e .
```

## Results

*Analysis in progress - figures and findings will be added as each pipeline step completes.*
