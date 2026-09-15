# Cellular Composition & Transcriptional Landscape: Primary vs Metastatic HCC

## Research Question
**How does the cellular composition and transcriptional landscape differ
between primary and metastatic hepatocellular carcinoma (HCC)?**

## Data
Re-analysis of public scRNA-seq data from:
> Lu et al. 2022, *A single-cell atlas of the multicellular ecosystem of
> primary and metastatic hepatocellular carcinoma*, Nature Communications
> 13:4594. https://doi.org/10.1038/s41467-022-32283-3

Processed data (GEO accession **GSE149614**) — 71,915 cells across 10
patients and 4 tissue types (NTL, PT, PVTT, MLN); 67,101 cells after QC.

## How to run
Everything is in one script: `scripts/full_analysis.R`. It is organized
into 15 numbered sections, each self-contained and commented with what it
does and what result it produces. Run top to bottom.

```r
# install once:
install.packages(c("Seurat", "tidyverse", "data.table", "R.utils",
                    "pheatmap", "patchwork", "BiocManager"))
BiocManager::install("GEOquery")

# then run:
source("scripts/full_analysis.R")
```

**Note:** Section 4 (building 21 per-sample Seurat objects) is the slowest
step (~20–40 minutes) and should not be interrupted once started.

## Results
See **[RESULTS_SUMMARY.md](RESULTS_SUMMARY.md)** for the full write-up
answering the research question, with explanations of every figure and
table.

- `results/` — CSV tables (composition stats, differential expression per
  cell type)
- `figures/` — PDF plots (UMAPs, composition bar chart, DE heatmap)

## Repository structure
```
├── README.md
├── RESULTS_SUMMARY.md      ← answers to the research question
├── scripts/
│   └── full_analysis.R     ← the complete, working pipeline
├── results/                ← output CSV tables
└── figures/                ← output plots
```

Raw/intermediate data files (`data/`) are not included in this repository
— they are large (multi-GB) and are re-downloadable via Section 1 of the
script directly from GEO (GSE149614).

## Citation
Lu, Y., Yang, A., Quan, C. et al. A single-cell atlas of the multicellular
ecosystem of primary and metastatic hepatocellular carcinoma. *Nat Commun*
13, 4594 (2022). https://doi.org/10.1038/s41467-022-32283-3
