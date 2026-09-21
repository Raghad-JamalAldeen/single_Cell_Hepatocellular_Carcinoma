# Primary vs. Metastatic HCC — Single-Cell RNA-seq Analysis

## Overview
This project re-analyzes public single-cell RNA-seq (scRNA-seq) data to
investigate how the cellular composition and transcriptional landscape of
hepatocellular carcinoma (HCC) differ between the primary tumor and
metastatic sites.

**Research question:** How does the cellular composition and
transcriptional landscape differ between primary and metastatic
hepatocellular carcinoma?

## Dataset
Data are from Lu et al. (2022), *A single-cell atlas of the multicellular
ecosystem of primary and metastatic hepatocellular carcinoma*, **Nature
Communications** 13:4594 (https://doi.org/10.1038/s41467-022-32283-3).
Processed counts and cell-level metadata were obtained from GEO under
accession **[GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)**.

| | |
|---|---|
| Patients | 10 |
| Tissue types | NTL (non-tumor liver), PT (primary tumor), PVTT (portal vein tumor thrombus), MLN (metastatic lymph node) |
| Total cells (raw) | 71,915 |
| Cells after QC | 67,101 |
| Genes | 25,712 |

## Analysis pipeline
The full workflow is implemented in `scripts/full_analysis.R`, organized
into 15 commented, self-contained sections meant to be run top to bottom:

| # | Section | What it does |
|---|---|---|
| 0 | Setup | Load required libraries |
| 1 | Download | Fetch GSE149614 supplementary files from GEO |
| 2 | Decompress | Unpack the ~3.5 GB counts matrix to disk |
| 3 | Identify samples | Parse the 21 sample codes from cell barcodes |
| 4 | Build objects | Create one Seurat object per sample (memory-safe, streamed) |
| 5 | Merge & annotate | Merge all samples; attach patient/tissue/cell-type metadata |
| 6 | Quality control | Filter by gene count (200–8000) and mitochondrial % (<10%) |
| 7 | Normalize | Log-normalization, variable feature selection, scaling |
| 8 | PCA | Dimensionality reduction |
| 9 | Integration | Batch correction across patients (RPCA method) |
| 10 | Clustering | UMAP embedding + graph-based clustering (29 clusters) |
| 11 | Annotation | Assign cell types by majority vote against published labels |
| 12 | Composition analysis | Cell-type proportions, PT vs. metastatic sites |
| 13 | Differential expression | Gene expression differences per cell type, PT vs. PVTT |
| 14–15 | Figures | Generate all final plots |

**Note:** Section 4 (building 21 per-sample Seurat objects) is the slowest
step (~20–40 minutes) and should not be interrupted once started.
RPCA integration (Section 9) was used instead of the CCA method from the
original paper for speed/memory reasons on limited hardware.

## How to run
```r
install.packages(c("Seurat", "tidyverse", "data.table", "R.utils",
                    "pheatmap", "patchwork", "ggrepel", "BiocManager"))
BiocManager::install("GEOquery")
source("scripts/full_analysis.R")
```

## Repository structure
```
├── README.md
├── RESULTS_SUMMARY.md      ← full write-up answering the research question
├── scripts/
│   └── full_analysis.R     ← complete, commented analysis pipeline
├── results/                ← output CSV tables
└── figures/                ← output plots (01–08, numbered in viewing order)
```
Raw/intermediate data files (multi-GB) are not included in this
repository — they are re-downloadable via Section 1 of the script
directly from GEO (GSE149614).

## Figures

| File | Shows |
|---|---|
| `01_umap_by_celltype_paper.png` | UMAP colored by the original authors' cell-type labels (validation) |
| `02_final_umap_celltype.png` | UMAP colored by our own cell-type annotation |
| `03_final_umap_tissue.png` | UMAP colored by tissue of origin |
| `04_final_umap_celltype_tissue.png` | Cell type and tissue UMAPs side by side |
| `05_composition_stacked_bar.png` | Cell-type proportions per sample, by tissue — **Part 1 of the research question** |
| `06_barchart_top_genes_hepatocyte.png` | Top 15 DE genes (hepatocytes, PT vs PVTT), ranked by fold change |
| `07_dotplot_top_genes_hepatocyte.png` | Same top 15 genes, showing expression level and % of cells expressing, per tissue |
| `08_volcano_hepatocyte_PT_vs_PVTT.png` | All hepatocyte genes, PT vs PVTT — **strongest evidence for Part 2** |

## Results tables

| File | Contents |
|---|---|
| `composition_PT_vs_Metastatic.csv` | Per-cell-type proportions and Wilcoxon test p-values, PT vs. metastatic |
| `DE_Hepatocyte_PT_vs_PVTT.csv` | Full differential expression results for hepatocytes |
| `DE_Myeloid_PT_vs_PVTT.csv` | Full differential expression results for myeloid cells |
| `DE_TNK_PT_vs_PVTT.csv` | Full differential expression results for T/NK cells |

## Key findings

**1. Cellular composition:** No statistically significant differences in
major cell-type proportions were found between PT and metastatic sites
(PVTT/MLN) — all p > 0.05, Wilcoxon rank-sum test. A directional trend was
observed: T/NK cell proportion decreases at metastatic sites (17.5% → 9.9%)
while hepatocyte (39.2% → 48.2%) and myeloid (27.3% → 32.2%) proportions
increase. This is consistent with the original study's observation that
PVTT/MLN more closely resemble PT than NTL in overall composition.

**2. Transcriptional landscape:** Clear, statistically significant
differences were found *within* matched cell types:
- **Hepatocytes:** normal liver identity genes (*ALB*, *AKR1C1*,
  *ALDH1A1*) are higher in PT; EMT/stemness markers (*EPCAM*, *KRT19*,
  *CLDN6*) are higher in PVTT — consistent with a "pro-metastatic
  hepatocyte" phenotype.
- **Myeloid cells:** the immunosuppressive gene *MIF* is elevated in
  PVTT, suggesting a more immunosuppressive microenvironment at the
  metastatic site.
- **T/NK cells:** top DE genes (*CD24*, *PRAME*) are not typical lymphoid
  markers, likely reflecting ambient RNA contamination in the
  densely-packed PVTT tissue — a noted limitation.

**Conclusion:** Metastatic progression in HCC appears to be driven
primarily by transcriptional reprogramming within existing cell
populations, rather than by large shifts in overall cell-type
composition — consistent with the conclusions of Lu et al. (2022).

See **[RESULTS_SUMMARY.md](RESULTS_SUMMARY.md)** for the complete
write-up with full interpretation and limitations.

## Citation
Lu, Y., Yang, A., Quan, C. et al. A single-cell atlas of the multicellular
ecosystem of primary and metastatic hepatocellular carcinoma. *Nat Commun*
13, 4594 (2022). https://doi.org/10.1038/s41467-022-32283-3
