configfile: "configs/workflow.yaml"


FIGURES = config["results"]["figures"]
TABLES = config["results"]["tables"]


rule all:
    input:
        "data/processed/adata_programs.h5ad",
        FIGURES,
        TABLES,


rule download_geo:
    output:
        extracted="data/raw/GSE161529/.extracted",
        features="data/raw/GSE161529/GSE161529_features.tsv.gz",
    log:
        "logs/download_geo.log",
    shell:
        "mkdir -p logs && python scripts/download_geo.py > {log} 2>&1"


rule build_raw_adata:
    input:
        "data/raw/GSE161529/.extracted",
        "data/raw/GSE161529/GSE161529_features.tsv.gz",
    output:
        "data/processed/adata_raw.h5ad",
    threads:
        config["threads"]["integration"]
    log:
        "logs/build_adata.log",
    shell:
        "mkdir -p logs && python scripts/build_adata.py > {log} 2>&1"


rule qc:
    input:
        "data/processed/adata_raw.h5ad",
    output:
        "data/processed/adata_qc.h5ad",
    threads:
        config["threads"]["integration"]
    log:
        "logs/run_qc.log",
    shell:
        "mkdir -p logs && python scripts/run_qc.py > {log} 2>&1"


rule subset_luminal:
    input:
        "data/processed/adata_qc.h5ad",
    output:
        "data/processed/adata_luminal.h5ad",
    threads:
        config["threads"]["integration"]
    log:
        "logs/run_subset_luminal.log",
    shell:
        "mkdir -p logs && python scripts/run_subset_luminal.py > {log} 2>&1"


rule cnv:
    input:
        "data/processed/adata_luminal.h5ad",
    output:
        "data/processed/adata_luminal_cnv.h5ad",
    threads:
        config["threads"]["cnv"]
    log:
        "logs/run_cnv.log",
    shell:
        "mkdir -p logs && python scripts/run_cnv.py > {log} 2>&1"


rule pseudotime:
    input:
        "data/processed/adata_luminal_cnv.h5ad",
    output:
        "data/processed/adata_trajectory.h5ad",
    threads:
        config["threads"]["pseudotime"]
    log:
        "logs/run_pseudotime.log",
    shell:
        "mkdir -p logs && python scripts/run_pseudotime.py > {log} 2>&1"


rule directionality:
    input:
        "data/processed/adata_trajectory.h5ad",
    output:
        adata="data/processed/adata_directionality.h5ad",
        table="results/tables/directionality_comparison.csv",
    threads:
        config["threads"]["directionality"]
    log:
        "logs/run_directionality.log",
    shell:
        "mkdir -p logs && python scripts/run_directionality.py > {log} 2>&1"


rule de_programs:
    input:
        "data/processed/adata_directionality.h5ad",
    output:
        adata="data/processed/adata_programs.h5ad",
        bins="results/tables/de_pseudotime_bins.csv",
        early_late="results/tables/de_early_vs_late.csv",
    threads:
        config["threads"]["de"]
    log:
        "logs/run_de_programs.log",
    shell:
        "mkdir -p logs && python scripts/run_de_programs.py > {log} 2>&1"


rule xenium_mapping:
    input:
        "results/tables/de_early_vs_late.csv",
    output:
        "results/tables/xenium_panel_trajectory_genes.csv",
    log:
        "logs/run_xenium_mapping.log",
    shell:
        "mkdir -p logs && python scripts/run_xenium_mapping.py > {log} 2>&1"


rule figures:
    input:
        "data/processed/adata_qc.h5ad",
        "data/processed/adata_programs.h5ad",
        "results/tables/de_early_vs_late.csv",
    output:
        FIGURES,
    threads:
        config["threads"]["figures"]
    log:
        "logs/make_figures.log",
    shell:
        "mkdir -p logs && python scripts/make_figures.py > {log} 2>&1"
