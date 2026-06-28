"""Marker panels, sample metadata, cell-type labels, and color palettes."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Marker panels
# ---------------------------------------------------------------------------

LUMINAL_POS = ["KRT8", "KRT18", "KRT19", "ESR1", "GATA3", "PGR", "FOXA1"]
LUMINAL_NEG = ["KRT14", "KRT5", "TP63"]

BASAL_POS = ["KRT14", "KRT5", "TP63", "KRT17"]
IMMUNE_POS = ["PTPRC"]  # CD45
ENDOTHELIAL_POS = ["PECAM1", "VWF"]  # CD31
FIBROBLAST_POS = ["COL1A1", "DCN", "FAP"]

PROGENITOR_MARKERS = ["KRT15", "KIT", "ALDH1A3"]

LUMINAL_PROGENITOR_POS = ["KRT8", "KRT18", "KIT", "ALDH1A3"]
MATURE_LUMINAL_POS = ["KRT8", "KRT18", "ESR1", "PGR", "FOXA1"]

PROLIFERATION_MARKERS = ["MKI67", "TOP2A", "PCNA"]

# Broad cell-type scoring panels (for sc.tl.score_genes)
CELL_TYPE_MARKERS: dict[str, list[str]] = {
    "luminal_progenitor": LUMINAL_PROGENITOR_POS,
    "mature_luminal": MATURE_LUMINAL_POS,
    "basal": BASAL_POS,
    "immune": IMMUNE_POS,
    "endothelial": ENDOTHELIAL_POS,
    "fibroblast": FIBROBLAST_POS,
}

# ---------------------------------------------------------------------------
# Sample metadata (from SampleStats.txt, GitHub yunshun/HumanBreast10X)
# ---------------------------------------------------------------------------

TISSUE_TYPES = {
    "normal": "Normal",
    "brca1": "BRCA1 preneoplastic",
    "tnbc": "TNBC",
    "tnbc_brca1": "TNBC (BRCA1)",
    "her2": "HER2+",
    "er": "ER+",
    "er_male": "ER+ (male)",
}

# GSM ID -> (sample_name, tissue_type, patient_id, fraction)
# fraction: "total", "epi" (sorted epithelial), "T" (tumor), "LN" (lymph node)
SAMPLE_METADATA: dict[str, tuple[str, str, str, str]] = {
    # Normal
    "GSM4909253": ("N-0092-total", "normal", "0092", "total"),
    "GSM4909254": ("N-0019-total", "normal", "0019", "total"),
    "GSM4909255": ("N-0280-epi", "normal", "0280", "epi"),
    "GSM4909256": ("N-0093-epi", "normal", "0093", "epi"),
    "GSM4909257": ("N-0093-total", "normal", "0093", "total"),
    "GSM4909258": ("N-1469-epi", "normal", "1469", "epi"),
    "GSM4909259": ("N-0408-epi", "normal", "0408", "epi"),
    "GSM4909260": ("N-1105-epi", "normal", "1105", "epi"),
    "GSM4909261": ("N-0230.17-total", "normal", "0230.17", "total"),
    "GSM4909262": ("N-0064-epi", "normal", "0064", "epi"),
    "GSM4909263": ("N-0064-total", "normal", "0064", "total"),
    "GSM4909264": ("N-0230.16-epi", "normal", "0230.16", "epi"),
    "GSM4909265": ("N-0233-total", "normal", "0233", "total"),
    "GSM4909266": ("N-0169-total", "normal", "0169", "total"),
    "GSM4909267": ("N-0123-epi", "normal", "0123", "epi"),
    "GSM4909268": ("N-0123-total", "normal", "0123", "total"),
    "GSM4909269": ("N-0342-epi", "normal", "0342", "epi"),
    "GSM4909270": ("N-0342-total", "normal", "0342", "total"),
    "GSM4909271": ("N-0288-total", "normal", "0288", "total"),
    "GSM4909272": ("N-0021-total", "normal", "0021", "total"),
    "GSM4909273": ("N-0275-epi", "normal", "0275", "epi"),
    "GSM4909274": ("N-0275-total", "normal", "0275", "total"),
    "GSM4909275": ("N-0372-epi", "normal", "0372", "epi"),
    "GSM4909276": ("N-0372-total", "normal", "0372", "total"),
    # BRCA1 preneoplastic
    "GSM4909277": ("B1-0894", "brca1", "0894", "total"),
    "GSM4909278": ("B1-0033", "brca1", "0033", "total"),
    "GSM4909279": ("B1-0023", "brca1", "0023", "total"),
    "GSM4909280": ("B1-0090", "brca1", "0090", "total"),
    # TNBC
    "GSM4909281": ("TN-0126", "tnbc", "0126", "total"),
    "GSM4909282": ("TN-0135", "tnbc", "0135", "total"),
    "GSM4909283": ("TN-0106", "tnbc", "0106", "total"),
    "GSM4909284": ("TN-0114-T2", "tnbc", "0114", "total"),
    # TNBC (BRCA1)
    "GSM4909285": ("TN-B1-4031", "tnbc_brca1", "4031", "total"),
    "GSM4909286": ("TN-B1-0131", "tnbc_brca1", "0131", "total"),
    "GSM4909287": ("TN-B1-0554", "tnbc_brca1", "0554", "total"),
    "GSM4909288": ("TN-B1-0177", "tnbc_brca1", "0177", "total"),
    # HER2+
    "GSM4909289": ("HER2-0308", "her2", "0308", "total"),
    "GSM4909290": ("HER2-0337", "her2", "0337", "total"),
    "GSM4909291": ("HER2-0031", "her2", "0031", "total"),
    "GSM4909292": ("HER2-0069", "her2", "0069", "total"),
    "GSM4909293": ("HER2-0161", "her2", "0161", "total"),
    "GSM4909294": ("HER2-0176", "her2", "0176", "total"),
    # ER+
    "GSM4909295": ("ER-0319", "er", "0319", "total"),
    "GSM4909296": ("ER-0001", "er", "0001", "total"),
    "GSM4909297": ("ER-0125", "er", "0125", "total"),
    "GSM4909298": ("ER-0360", "er", "0360", "total"),
    "GSM4909299": ("ER-0114-T3", "er", "0114", "total"),
    "GSM4909300": ("ER-0032", "er", "0032", "total"),
    "GSM4909301": ("ER-0042", "er", "0042", "total"),
    "GSM4909302": ("ER-0025", "er", "0025", "total"),
    "GSM4909303": ("ER-0151", "er", "0151", "total"),
    "GSM4909304": ("ER-0163", "er", "0163", "total"),
    "GSM4909305": ("ER-0029-7C", "er", "0029", "total"),
    "GSM4909306": ("ER-0029-9C", "er", "0029", "total"),
    "GSM4909307": ("ER-0040-T", "er", "0040", "T"),
    "GSM4909308": ("ER-0040-LN", "er", "0040", "LN"),
    "GSM4909309": ("ER-0043-T", "er", "0043", "T"),
    "GSM4909310": ("ER-0043-LN", "er", "0043", "LN"),
    "GSM4909311": ("ER-0056-T", "er", "0056", "T"),
    "GSM4909312": ("ER-0056-LN", "er", "0056", "LN"),
    "GSM4909313": ("ER-0064-T", "er", "0064", "T"),
    "GSM4909314": ("ER-0064-LN", "er", "0064", "LN"),
    "GSM4909315": ("ER-0167-T", "er", "0167", "T"),
    "GSM4909316": ("ER-0167-LN", "er", "0167", "LN"),
    "GSM4909317": ("ER-0173-T", "er", "0173", "T"),
    "GSM4909318": ("ER-0173-LN", "er", "0173", "LN"),
    # ER+ (male)
    "GSM4909319": ("mER-0178", "er_male", "0178", "total"),
    "GSM4909320": ("mER-0068-T", "er_male", "0068", "T"),
    "GSM4909321": ("mER-0068-LN", "er_male", "0068", "LN"),
}

# ---------------------------------------------------------------------------
# Coarse tissue groupings for analysis
# ---------------------------------------------------------------------------

NORMAL_TYPES = {"normal"}
PRENEOPLASTIC_TYPES = {"brca1"}
MALIGNANT_TYPES = {"tnbc", "tnbc_brca1", "her2", "er", "er_male"}

# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------

TISSUE_COLORS: dict[str, str] = {
    "normal": "#2ca02c",
    "brca1": "#ff7f0e",
    "tnbc": "#d62728",
    "tnbc_brca1": "#9467bd",
    "her2": "#e377c2",
    "er": "#1f77b4",
    "er_male": "#17becf",
}

MALIGNANCY_COLORS: dict[str, str] = {
    "normal": "#2ca02c",
    "preneoplastic": "#ff7f0e",
    "malignant": "#d62728",
}

# ---------------------------------------------------------------------------
# 380-gene Xenium breast panel (for cross-project mapping to
# segmentation-benchmark). Full list loaded at runtime from the
# segmentation-benchmark repo if available; these are the trajectory-relevant
# subset for quick reference.
# ---------------------------------------------------------------------------

XENIUM_TRAJECTORY_GENES = [
    "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "KRT15", "KRT17",
    "ESR1", "PGR", "GATA3", "FOXA1",
    "EPCAM", "MKI67", "TOP2A",
    "PTPRC", "PECAM1", "VWF",
    "COL1A1", "DCN",
]
