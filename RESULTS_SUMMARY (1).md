# Results: Primary vs. Metastatic Hepatocellular Carcinoma

## Research Question
How does the cellular composition and transcriptional landscape differ
between primary and metastatic hepatocellular carcinoma (HCC)?

## Dataset
Single-cell RNA-seq data from Lu et al. (2022, *Nat Commun* 13:4594; GEO
accession GSE149614) were re-analyzed. The dataset comprises 67,101 cells
(post-QC) from 10 HCC patients across four tissue types: non-tumor liver
(NTL), primary tumor (PT), portal vein tumor thrombus (PVTT), and
metastatic lymph node (MLN).

## Figure Index

| File | Description |
|---|---|
| `01_umap_by_celltype_paper.png` | UMAP colored by the original authors' cell-type labels — validates our clustering against the published atlas |
| `02_final_umap_celltype.png` | UMAP colored by our own cell-type annotation |
| `03_final_umap_tissue.png` | UMAP colored by tissue of origin (NTL/PT/PVTT/MLN) |
| `04_final_umap_celltype_tissue.png` | Cell type and tissue UMAPs side by side |
| `05_composition_stacked_bar.png` | Cell-type proportions per sample, by tissue — **answers Part 1** |
| `06_barchart_top_genes_hepatocyte.png` | Top 15 DE genes, PT vs PVTT hepatocytes, ranked by fold change |
| `07_dotplot_top_genes_hepatocyte.png` | Same top 15 genes, showing expression level and percent-expressing cells per tissue |
| `08_volcano_hepatocyte_PT_vs_PVTT.png` | All hepatocyte genes, PT vs PVTT — **strongest visual evidence for Part 2** |

---

## Part 1: Cellular Composition

### Figure: `figures/05_composition_stacked_bar.png`
Stacked bar plot showing the proportion of each major cell type
(Hepatocyte, T/NK, Myeloid, B, Endothelial, Fibroblast) per sample,
faceted by tissue type.

### Table: `results/composition_PT_vs_Metastatic.csv`

| Cell type | Mean proportion, PT | Mean proportion, Metastatic | Direction | p-value |
|---|---|---|---|---|
| T/NK | 17.5% | 9.9% | Decreased in metastatic sites | 0.469 |
| Hepatocyte | 39.2% | 48.2% | Increased in metastatic sites | 0.692 |
| Myeloid | 27.3% | 32.2% | Increased in metastatic sites | 0.692 |
| Fibroblast | 4.1% | 2.8% | Decreased in metastatic sites | 0.811 |
| Endothelial | 6.1% | 1.2% | Decreased in metastatic sites | 0.864 |
| B | 6.4% | 5.7% | Comparable | 0.937 |

*Metastatic group = PVTT + MLN samples pooled. Statistical comparison by
two-sided Wilcoxon rank-sum test on per-sample proportions.*

### Interpretation
No statistically significant differences were observed (all p > 0.05),
consistent with the modest sample size (only 3 patients contributed PVTT
tissue, and 1 patient contributed MLN tissue), which limits statistical
power for compositional comparisons. However, the directional trends are
consistent with the original study: T/NK cell proportion trends lower in
metastatic sites, while malignant hepatocytes and myeloid cells trend
higher. This aligns with the original report that PVTT and MLN more
closely resemble PT than NTL in overall cell-type composition — indicating
that the dominant biological differences between primary and metastatic
HCC are unlikely to be driven by gross shifts in cell-type abundance, but
rather by transcriptional reprogramming within each cell type (see Part 2).

---

## Part 2: Transcriptional Landscape

Differential expression (DE) analysis was performed within each major cell
type, comparing PT against PVTT (Wilcoxon rank-sum test, min.pct = 0.1,
|log2FC| threshold = 0.25).

### 1. Hepatocytes (malignant epithelial cells)
**Table:** `results/DE_Hepatocyte_PT_vs_PVTT.csv`
**Figures:** `figures/06_barchart_top_genes_hepatocyte.png`,
`figures/07_dotplot_top_genes_hepatocyte.png`,
`figures/08_volcano_hepatocyte_PT_vs_PVTT.png`

- **Upregulated in PT:** *ALB*, *AKR1C1*, *AKR1C3*, *ALDH1A1*, *GC*,
  *CFHR1* — genes associated with normal hepatocyte metabolic function
  (notably *ALB*, encoding albumin).
- **Upregulated in PVTT:** *EPCAM*, *KRT19*, *CLDN6*, *FBLN1*, *TUBB2B*,
  *WFDC2* — genes associated with stemness and epithelial-mesenchymal
  transition (EMT), a process implicated in tumor invasion and metastatic
  spread.

**Interpretation:** Malignant hepatocytes at the metastatic site show
reduced expression of normal hepatocyte identity genes and increased
expression of stemness/EMT markers, consistent with the "pro-metastatic
hepatocyte" phenotype described in the original study.

### 2. Myeloid cells
**Table:** `results/DE_Myeloid_PT_vs_PVTT.csv`

- **Upregulated in PVTT:** ***MIF*** — an immunosuppressive ligand
  reported in the original study to signal through CD74 on
  tumor-associated macrophages — along with *S100A8*, *S100A10*, and
  *VIM*.
- **Upregulated in PT:** *ALB*, *APOC3*, *APOH*, *TTR* — likely reflecting
  ambient RNA contamination from neighboring hepatocytes rather than
  genuine myeloid expression.

**Interpretation:** Elevated *MIF* expression in myeloid cells at the
metastatic site supports a more immunosuppressive microenvironment in
regions of tumor spread.

### 3. T/NK cells
**Table:** `results/DE_TNK_PT_vs_PVTT.csv`

- **Upregulated in PVTT:** *CD24*, *PRAME*, *CRABP1*, *MDK* — genes not
  typically associated with lymphoid identity, more consistent with
  ambient RNA contamination from adjacent malignant cells in the
  densely tumor-infiltrated PVTT microenvironment. This result should be
  interpreted with caution and flags a limitation of the analysis.

---

## Summary: Answer to the Research Question

1. **Cellular composition:** No statistically significant differences in
   major cell-type proportions were detected between PT and metastatic
   sites (PVTT/MLN), though a trend toward reduced T/NK infiltration and
   increased hepatocyte/myeloid abundance was observed. Cell-type
   composition alone does not robustly distinguish primary from
   metastatic HCC in this dataset.

2. **Transcriptional landscape:** Clear, statistically significant
   transcriptional differences were identified within matched cell types.
   Malignant hepatocytes at metastatic sites downregulate normal liver
   identity genes and upregulate EMT/stemness markers. Myeloid cells at
   metastatic sites upregulate the immunosuppressive gene *MIF*.

**Conclusion:** Metastatic progression in HCC appears to be driven
primarily by transcriptional reprogramming within existing cell
populations — rather than by wholesale changes in the cellular composition
of the tumor microenvironment. This is consistent with the conclusions of
the original study (Lu et al., 2022).

---

## Limitations
- Small patient cohort (n = 10; only 3 patients contributed PVTT tissue
  and 1 contributed MLN tissue), limiting statistical power for
  compositional comparisons.
- Some DE results (notably in T/NK cells) may be confounded by ambient
  RNA contamination from adjacent malignant cells.
- RPCA integration was used in place of the CCA method used in the
  original study, due to memory constraints; this yielded 29 clusters
  versus the 53 reported originally, though major cell-type annotations
  (validated against the authors' own per-cell labels) closely matched
  published proportions.
