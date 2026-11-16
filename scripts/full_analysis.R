## =================================================================
## HCC scRNA-seq Analysis: Primary vs Metastatic Hepatocellular
## Carcinoma — Cellular Composition & Transcriptional Landscape
##
## Data source: Lu et al. 2022, Nat Commun 13:4594
## https://doi.org/10.1038/s41467-022-32283-3
## Processed data: GEO accession GSE149614
##
## Research question:
## How does the cellular composition and transcriptional landscape
## differ between primary and metastatic hepatocellular carcinoma?
##
## HOW TO USE THIS SCRIPT:
## Run section by section, in order, from top to bottom. Each
## section is self-contained and prints/saves a specific result
## (a table, a figure, or a checkpoint .rds file). Do not skip
## sections — later sections depend on objects created earlier.
##
## IMPORTANT — before running: edit the path marked
## "CHANGE THIS" below to match your own machine.
## =================================================================


## =================================================================
## SECTION 0 — Setup: libraries and library path
## =================================================================

## CHANGE THIS: set to wherever your R packages are installed.
## (Only needed if you installed packages to a custom location,
## e.g. a drive other than C: due to limited disk space.)
# .libPaths("D:/R-library")

library(GEOquery)
library(Seurat)
library(tidyverse)
library(data.table)
library(R.utils)
library(pheatmap)
library(patchwork)
library(ggrepel)


## =================================================================
## SECTION 1 — Download data from GEO
## Result: 3 files appear in data/
## =================================================================

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

getGEOSuppFiles(
  GEO           = "GSE149614",
  baseDir       = "data",
  makeDirectory = FALSE
)

list.files("data")
## Expect to see:
##   GSE149614_HCC.metadata.updated.txt.gz
##   GSE149614_HCC.scRNAseq.S71915.count.txt.gz
##   GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz


## =================================================================
## SECTION 2 — Decompress the counts file
## Result: data/unzipped/GSE149614_HCC.scRNAseq.S71915.count.txt
## =================================================================

dir.create("data/unzipped", showWarnings = FALSE)

counts_txt_path <- "data/unzipped/GSE149614_HCC.scRNAseq.S71915.count.txt"

if (!file.exists(counts_txt_path)) {
  gunzip(
    "data/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
    destname  = counts_txt_path,
    remove    = FALSE,
    overwrite = TRUE
  )
}

file.exists(counts_txt_path)   # should be TRUE


## =================================================================
## SECTION 3 — Identify samples from the counts file header
## Result: sample_prefixes = character vector of 21 sample codes
## =================================================================

header    <- fread(counts_txt_path, nrows = 0)
col_names <- colnames(header)
rm(header); gc()

## Sample codes are embedded in the cell barcode prefix, e.g.
## "HCC01T_AAACCTGAGGGCATGT" -> sample "HCC01T"
## Suffix meaning: T = primary Tumor (PT), N = Normal liver (NTL),
##                 P = PVTT, L = Lymph node metastasis (MLN)
sample_prefixes <- unique(sub("_.*", "", col_names[-1]))

length(sample_prefixes)   # should be 21
sample_prefixes


## =================================================================
## SECTION 4 — Build one Seurat object per sample (memory-safe)
## Result: data/seurat_list_21_samples.rds
## This section is the slowest step (~20-40 min). Do not interrupt
## it once started.
## =================================================================

seurat_list <- list()

for (s in sample_prefixes) {

  message("Processing sample: ", s)

  cols_s <- c("V1", col_names[grepl(paste0("^", s, "_"), col_names)])

  counts_s <- fread(counts_txt_path, select = cols_s)

  gene_names_s <- counts_s$V1
  mat_s <- as.matrix(counts_s[, -1])
  rownames(mat_s) <- gene_names_s
  rm(counts_s)

  seurat_list[[s]] <- CreateSeuratObject(
    counts       = mat_s,
    project      = s,
    min.cells    = 0,
    min.features = 0
  )

  rm(mat_s)
  gc()
}

message("All samples processed. Total objects: ", length(seurat_list))

saveRDS(seurat_list, "data/seurat_list_21_samples.rds")


## =================================================================
## SECTION 5 — Merge samples, fix cell names, attach metadata
## Result: data/seu_merged_raw.rds  (25,712 genes x 71,915 cells)
## =================================================================

seurat_list <- readRDS("data/seurat_list_21_samples.rds")

seu_merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

## Fix duplicated sample-name prefix introduced by add.cell.ids
## (original barcodes already contained the sample name)
new_names <- sub("^([^_]+)_\\1_", "\\1_", colnames(seu_merged))
seu_merged <- RenameCells(seu_merged, new.names = new_names)

## Attach the authors' cell-level metadata
cell_meta <- fread("data/GSE149614_HCC.metadata.updated.txt.gz")
cell_meta_ordered <- cell_meta[match(colnames(seu_merged), cell_meta$Cell), ]

stopifnot(sum(is.na(cell_meta_ordered$Cell)) == 0)  # sanity check

seu_merged$patient        <- cell_meta_ordered$patient
seu_merged$site           <- cell_meta_ordered$site
seu_merged$celltype_paper <- cell_meta_ordered$celltype
seu_merged$stage          <- cell_meta_ordered$stage
seu_merged$virus          <- cell_meta_ordered$virus

## Standardize tissue labels to match paper terminology
seu_merged$tissue <- recode(seu_merged$site,
  "Normal" = "NTL",
  "Tumor"  = "PT",
  "PVTT"   = "PVTT",
  "Lymph"  = "MLN"
)

table(seu_merged$tissue)          # NTL=28687 PT=34414 PVTT=5971 MLN=2843
table(seu_merged$celltype_paper)  # matches paper's Fig. 1 counts

saveRDS(seu_merged, "data/seu_merged_raw.rds")


## =================================================================
## SECTION 6 — Quality Control
## Result: data/seu_merged_after_qc.rds (67,101 cells retained)
## =================================================================

seu_merged <- readRDS("data/seu_merged_raw.rds")

seu_merged[["percent.mt"]] <- PercentageFeatureSet(seu_merged, pattern = "^MT-")

VlnPlot(seu_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0)

## Same thresholds used by the original paper
seu_merged <- subset(
  seu_merged,
  subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10
)

dim(seu_merged)   # 25712 67101

saveRDS(seu_merged, "data/seu_merged_after_qc.rds")


## =================================================================
## SECTION 7 — Normalization, feature selection, scaling
## Result: data/seu_merged_after_scaling.rds
## =================================================================

seu_merged <- readRDS("data/seu_merged_after_qc.rds")

## Split by patient so each patient can be normalized/integrated
## as its own batch
seu_merged[["RNA"]] <- split(seu_merged[["RNA"]], f = seu_merged$patient)

seu_merged <- NormalizeData(seu_merged)
seu_merged <- FindVariableFeatures(seu_merged)

## Scale only the ~2000 variable genes (scaling all 25,712 genes
## can exhaust RAM on machines with 16 GB)
seu_merged <- ScaleData(seu_merged, features = VariableFeatures(seu_merged))

saveRDS(seu_merged, "data/seu_merged_after_scaling.rds")


## =================================================================
## SECTION 8 — PCA
## Result: data/seu_merged_after_pca.rds
## =================================================================

seu_merged <- readRDS("data/seu_merged_after_scaling.rds")

seu_merged <- RunPCA(seu_merged)

saveRDS(seu_merged, "data/seu_merged_after_pca.rds")


## =================================================================
## SECTION 9 — Integration across patients (batch correction)
## Result: data/seu_merged_after_integration.rds
##
## NOTE: CCAIntegration can hang/run very slowly on large datasets
## with limited RAM. RPCAIntegration is much faster and gave a
## clean result in this analysis.
## =================================================================

seu_merged <- readRDS("data/seu_merged_after_pca.rds")

seu_merged <- IntegrateLayers(
  object         = seu_merged,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca"
)

saveRDS(seu_merged, "data/seu_merged_after_integration.rds")


## =================================================================
## SECTION 10 — UMAP + Clustering
## Result: data/seu_merged_after_clustering.rds
## =================================================================

seu_merged <- readRDS("data/seu_merged_after_integration.rds")

seu_merged <- RunUMAP(seu_merged, reduction = "integrated.rpca", dims = 1:30)
seu_merged <- FindNeighbors(seu_merged, reduction = "integrated.rpca", dims = 1:30)
seu_merged <- FindClusters(seu_merged, resolution = 1.0)

## 29 clusters found

saveRDS(seu_merged, "data/seu_merged_after_clustering.rds")


## =================================================================
## SECTION 11 — Manual marker review and published-label annotation
##
## Manual cluster labels were assigned from marker expression and
## compared with published per-cell labels for the reviewed run.
## Published labels are retained for downstream analyses.
## Lineage annotation alone does not establish malignancy.
## =================================================================

seu_merged <- readRDS("data/seu_merged_after_clustering.rds")

DefaultAssay(seu_merged) <- "RNA"

seu_merged <- JoinLayers(seu_merged)

Idents(seu_merged) <- "seurat_clusters"

dir.create(
  "results/annotation",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)

# Preserve the original authors' per-cell labels.
if (!"celltype_paper" %in% names(seu_merged[[]])) {

  stopifnot(
    "celltype" %in% names(seu_merged[[]])
  )

  seu_merged$celltype_paper <- seu_merged$celltype
}

stopifnot(
  !anyNA(seu_merged$celltype_paper)
)


## 11.1 — Identify markers without using published identities
## Use up to 200 cells per cluster for this descriptive screen.

set.seed(123)

cells_by_cluster <- split(
  colnames(seu_merged),
  as.character(seu_merged$seurat_clusters)
)

marker_cells <- unlist(
  lapply(
    cells_by_cluster,
    function(x) {
      sample(x, min(length(x), 200))
    }
  ),
  use.names = FALSE
)

marker_counts <- as(
  seu_merged[["RNA"]]$counts[
    ,
    marker_cells,
    drop = FALSE
  ],
  "dgCMatrix"
)

marker_obj <- CreateSeuratObject(
  counts = marker_counts,
  min.cells = 0,
  min.features = 0
)

marker_obj$cluster <- as.character(
  seu_merged$seurat_clusters[
    match(
      colnames(marker_obj),
      colnames(seu_merged)
    )
  ]
)

marker_obj <- NormalizeData(marker_obj)

Idents(marker_obj) <- "cluster"

annotation_markers <- FindAllMarkers(
  marker_obj,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.5,
  max.cells.per.ident = 200,
  random.seed = 123,
  densify = FALSE
)

if (nrow(annotation_markers) == 0) {
  stop(
    "No markers returned. Check the marker-test warnings."
  )
}

write.csv(
  annotation_markers,
  "results/annotation/all_cluster_markers.csv",
  row.names = FALSE
)

top20 <- do.call(
  rbind,
  lapply(
    split(
      annotation_markers,
      annotation_markers$cluster
    ),
    function(x) {
      head(
        x[order(-x$avg_log2FC), ],
        20
      )
    }
  )
)

write.csv(
  top20,
  "results/annotation/top20_markers_per_cluster.csv",
  row.names = FALSE
)

rm(marker_counts, marker_obj)
gc()


## 11.2 — Review canonical-marker expression across all cells

annotation_genes <- c(
  "ALB", "SERPINA1", "HNF4A",
  "EPCAM", "KRT19",
  "CD3D", "CD3E", "NKG7", "GNLY",
  "LYZ", "CD14", "CD68", "CD163",
  "CD1C", "CLEC9A", "CLEC4C",
  "MS4A1", "CD79A", "JCHAIN", "MZB1",
  "VWF", "PECAM1",
  "COL1A1", "COL1A2", "DCN",
  "RGS5", "PDGFRB", "TPSAB1"
)

annotation_genes <- intersect(
  annotation_genes,
  rownames(seu_merged)
)

if (length(annotation_genes) == 0) {
  stop("No canonical markers found. Check gene identifiers.")
}

if (!"data" %in% Layers(seu_merged[["RNA"]])) {
  seu_merged <- NormalizeData(seu_merged)
}

p_annotation <- DotPlot(
  seu_merged,
  features = annotation_genes,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

ggsave(
  "figures/annotation_markers_by_cluster.png",
  plot = p_annotation,
  width = 14,
  height = 10,
  dpi = 150
)


## 11.3 — Record our marker-based manual decisions
## These decisions apply to the reviewed 36-cluster run only.
## Plasma cells are grouped under B for broad-label comparison.

manual_review <- read.table(
  header = TRUE,
  sep = "|",
  colClasses = "character",
  text = "cluster|manual_label|supporting_markers
0|T/NK|TRGC1;TRGC2;TRDC;CD8A
1|Hepatocyte|CYP2D6;BHMT;FGL1;AHSG
2|T/NK|CD8A;CD8B;GZMK;CD3E
3|Myeloid|C1QA;C1QC;FOLR2;CD209
4|T/NK|IL7R;CD40LG;CCR7;CD3G
5|T/NK|KLRF1;NKG7;GNLY;TRDC
6|Hepatocyte|ARG1;FABP1;APOA2;OTC
7|Myeloid|CD5L;MARCO;TIMD4;CD163
8|Myeloid|SPP1;GPNMB;TREM2;MMP12
9|Epithelial|EPCAM;KRT19;CHGA;CHGB
10|Myeloid|FCN1;S100A8;S100A9;LST1
11|Hepatocyte|HNF4A;FABP1;APOA2;CDC45
12|Endothelial|VWF;PLVAP;CD34;ROBO4
13|B|JCHAIN;MZB1;TNFRSF17;DERL3
14|T/NK|FOXP3;IL2RA;CTLA4;TRAC
15|Epithelial|HNF4A;FABP1;EPCAM;KRT19
16|Myeloid|CD1C;CLEC10A;FCER1A;CD1E
17|B|MS4A1;CD19;CD79A;BANK1
18|Fibroblast|COL1A1;COL1A2;COL3A1;DCN
19|Hepatocyte|HNF4A;APOA2;CCNB1;PLK1
20|Hepatocyte|SLCO1B3;SLC22A1;CYP2C9;HAMP
21|T/NK|CD3D;CD3E;TRAC;MKI67
22|Endothelial|CLEC4G;CLEC4M;STAB2;OIT3
23|Mural|MYH11;CNN1;TAGLN;ACTA2
24|Endothelial|GJA5;SEMA3G;BMX;CLDN5
25|Myeloid|FCGR1A;TREM2;F13A1;MKI67
26|Myeloid|CLEC9A;XCR1;BATF3;FLT3
27|T/NK|CD8A;CD8B;CD3D;LAG3
28|B|MZB1;TNFRSF17;DERL3;RRM2
29|T/NK|CD6;TRAT1;IL7R;IFNG
30|Epithelial|KRT7;CFTR;MMP7;PKHD1
31|Unresolved|PECAM1;CDH5;RGS5;PDGFRB
32|Myeloid|LILRA4;CLEC4C;SPIB;GZMB
33|Myeloid|TPSAB1;TPSB2;CPA3;KIT
34|B|MZB1;JCHAIN;XBP1;CD79A
35|Myeloid|CD5L;MARCO;TIMD4;CD163"
)

manual_review$needs_review <-
  manual_review$cluster %in%
  as.character(c(0, 5, 9, 11, 15, 19, 25, 28, 31))


## 11.4 — Check exact cell-to-cluster membership
## Do not transfer manual labels merely because cluster IDs match.

cell_names <- colnames(seu_merged)

cluster_ids <- as.character(
  seu_merged$seurat_clusters
)

stopifnot(
  !anyNA(cluster_ids),
  anyDuplicated(cell_names) == 0
)

ord <- order(
  cell_names,
  method = "radix"
)

clustering_text <- paste(
  paste(
    cell_names[ord],
    cluster_ids[ord],
    sep = "\t"
  ),
  collapse = "\n"
)

fingerprint_file <- tempfile("annotation-clustering-")

writeBin(
  charToRaw(enc2utf8(clustering_text)),
  fingerprint_file
)

current_fingerprint <- unname(
  tools::md5sum(fingerprint_file)
)

unlink(fingerprint_file)

reviewed_fingerprint <- "51d7ab76ca3e18deb4df04f28008c50d"


## 11.5 — Compare manual and published labels when applicable

seu_merged$celltype_manual <- NA_character_
seu_merged$annotation_needs_review <- NA

annotation_comparison <- data.frame(
  Manual = character(),
  Published = character(),
  Freq = integer()
)

comparison_percent <- data.frame()

annotation_agreement <- data.frame(
  status = "Manual comparison not performed: different clustering",
  cells_compared = 0L,
  percent_agreement = NA_real_
)

if (current_fingerprint == reviewed_fingerprint) {

  idx <- match(
    cluster_ids,
    manual_review$cluster
  )

  stopifnot(!anyNA(idx))

  seu_merged$celltype_manual <-
    manual_review$manual_label[idx]

  seu_merged$annotation_needs_review <-
    manual_review$needs_review[idx]

  comparison <- table(
    Manual = seu_merged$celltype_manual,
    Published = seu_merged$celltype_paper
  )

  annotation_comparison <- as.data.frame(
    comparison
  )

  comparison_percent <- round(
    100 * prop.table(comparison, margin = 1),
    1
  )

  annotation_agreement <- data.frame(
    status = "Compared marker-based manual labels with published labels",
    cells_compared = ncol(seu_merged),
    percent_agreement = 100 * mean(
      seu_merged$celltype_manual ==
        seu_merged$celltype_paper
    )
  )

} else {

  warning(
    "Different clustering: manual labels were NOT transferred. ",
    "Published labels will still be used downstream."
  )
}

manual_review$reviewed_clustering_fingerprint <-
  reviewed_fingerprint

write.csv(
  manual_review,
  "results/annotation/reviewed_manual_decisions.csv",
  row.names = FALSE
)

write.csv(
  annotation_comparison,
  "results/annotation/manual_vs_published_counts.csv",
  row.names = FALSE
)

write.csv(
  comparison_percent,
  "results/annotation/manual_vs_published_percent.csv"
)

write.csv(
  annotation_agreement,
  "results/annotation/agreement_all_cells.csv",
  row.names = FALSE
)

print(annotation_agreement)


## 11.6 — Retain published labels for downstream analyses
## Manual labels remain separate for transparency.

seu_merged$major_celltype <- seu_merged$celltype_paper

table(seu_merged$major_celltype)

saveRDS(
  seu_merged,
  "data/seu_final_annotated.rds"
)


## =================================================================
## SECTION 12 — RESEARCH QUESTION, PART 1: Cellular composition
## PT vs metastatic sites (PVTT, MLN)
## Result: results/composition_PT_vs_Metastatic.csv
##         figures/05_composition_stacked_bar.png
## =================================================================

seu_merged <- readRDS("data/seu_final_annotated.rds")

## Per-sample cell-type proportions
prop_df <- seu_merged@meta.data %>%
  group_by(patient, tissue, major_celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient, tissue) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

## Visual: stacked bar of composition per sample, faceted by tissue
p_bar <- prop_df %>%
  ggplot(aes(x = patient, y = prop, fill = major_celltype)) +
  geom_col() +
  facet_grid(~ tissue, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(y = "Proportion of cells", x = NULL, fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figures/05_composition_stacked_bar.png", p_bar, width = 12, height = 6, dpi = 150)

## Statistical comparison: Primary (PT) vs Metastatic (PVTT+MLN)
comp_pt_vs_met <- prop_df %>%
  filter(tissue %in% c("PT", "PVTT", "MLN")) %>%
  mutate(group = ifelse(tissue == "PT", "Primary", "Metastatic")) %>%
  group_by(major_celltype) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(prop ~ group)$p.value,
      error = function(e) NA
    ),
    mean_primary    = mean(prop[group == "Primary"], na.rm = TRUE),
    mean_metastatic = mean(prop[group == "Metastatic"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p_value)

print(comp_pt_vs_met, width = Inf)

write_csv(comp_pt_vs_met, "results/composition_PT_vs_Metastatic.csv")


## =================================================================
## SECTION 13 — RESEARCH QUESTION, PART 2: Transcriptional
## landscape — differential expression within matched cell types,
## Primary (PT) vs Metastatic (PVTT)
## Result: results/DE_<CellType>_PT_vs_PVTT.csv  (x3)
## =================================================================

Idents(seu_merged) <- seu_merged$tissue

## ---- Hepatocytes (the malignant cells themselves) ----
hep_cells <- subset(seu_merged, subset = major_celltype == "Hepatocyte" &
                                          tissue %in% c("PT", "PVTT"))

de_hep_pt_vs_pvtt <- FindMarkers(
  hep_cells, ident.1 = "PT", ident.2 = "PVTT",
  test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25
)

write.csv(de_hep_pt_vs_pvtt, "results/DE_Hepatocyte_PT_vs_PVTT.csv")

## ---- Myeloid cells (macrophages etc.) ----
myeloid_cells <- subset(seu_merged, subset = major_celltype == "Myeloid" &
                                              tissue %in% c("PT", "PVTT"))

de_myeloid_pt_vs_pvtt <- FindMarkers(
  myeloid_cells, ident.1 = "PT", ident.2 = "PVTT",
  test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25
)

write.csv(de_myeloid_pt_vs_pvtt, "results/DE_Myeloid_PT_vs_PVTT.csv")

## ---- T/NK cells ----
tnk_cells <- subset(seu_merged, subset = major_celltype == "T/NK" &
                                          tissue %in% c("PT", "PVTT"))

de_tnk_pt_vs_pvtt <- FindMarkers(
  tnk_cells, ident.1 = "PT", ident.2 = "PVTT",
  test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25
)

write.csv(de_tnk_pt_vs_pvtt, "results/DE_TNK_PT_vs_PVTT.csv")


## =================================================================
## SECTION 14 — Final figures
## Result: figures/01_umap_by_celltype_paper.png
##         figures/02_final_umap_celltype.png
##         figures/03_final_umap_tissue.png
##         figures/04_final_umap_celltype_tissue.png
## =================================================================

## UMAP validation: colored by the authors' original cell-type labels
p0 <- DimPlot(seu_merged, reduction = "umap", group.by = "celltype_paper")
ggsave("figures/01_umap_by_celltype_paper.png", p0, width = 7, height = 6, dpi = 150)

## UMAP: our own cell-type annotation
p1 <- DimPlot(seu_merged, reduction = "umap", group.by = "major_celltype",
              label = TRUE) + NoLegend() + ggtitle("Cell types")
ggsave("figures/02_final_umap_celltype.png", p1, width = 7, height = 6, dpi = 150)

## UMAP: tissue of origin
p2 <- DimPlot(seu_merged, reduction = "umap", group.by = "tissue") +
  ggtitle("Tissue")
ggsave("figures/03_final_umap_tissue.png", p2, width = 7, height = 6, dpi = 150)

## Combined UMAP: cell type + tissue side by side
ggsave("figures/04_final_umap_celltype_tissue.png", p1 + p2, width = 12, height = 5, dpi = 150)


## =================================================================
## SECTION 15 — Top DE genes: bar chart, dot plot, volcano plot
## (Hepatocytes, PT vs PVTT)
## Result: figures/06_barchart_top_genes_hepatocyte.png
##         figures/07_dotplot_top_genes_hepatocyte.png
##         figures/08_volcano_hepatocyte_PT_vs_PVTT.png
## =================================================================

top_genes <- rownames(de_hep_pt_vs_pvtt[order(de_hep_pt_vs_pvtt$p_val_adj), ])[1:15]

hep_only <- subset(seu_merged, subset = major_celltype == "Hepatocyte" &
                                         tissue %in% c("PT", "PVTT"))

## ---- Bar chart: quick visual summary of the top 15 genes ----
top_genes_df <- de_hep_pt_vs_pvtt[top_genes, ] %>%
  rownames_to_column("gene") %>%
  mutate(gene = fct_reorder(gene, avg_log2FC))

p_bar_genes <- ggplot(top_genes_df, aes(x = gene, y = avg_log2FC, fill = avg_log2FC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "Higher in PT", "FALSE" = "Higher in PVTT")) +
  theme_minimal(base_size = 13) +
  labs(title = "Top 15 Differentially Expressed Genes",
       subtitle = "Hepatocytes: Primary Tumor vs PVTT",
       x = NULL, y = expression(log[2]~"Fold Change"),
       fill = NULL) +
  theme(legend.position = "top")

ggsave("figures/06_barchart_top_genes_hepatocyte.png", p_bar_genes, width = 8, height = 6, dpi = 150)

## ---- Dot plot: expression level + percent expressing, per tissue ----
p_dotplot <- DotPlot(hep_only, features = top_genes, group.by = "tissue") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 9)) +
  labs(title = "Top 15 DE Genes: Expression by Tissue", x = NULL, y = NULL)

ggsave("figures/07_dotplot_top_genes_hepatocyte.png", p_dotplot, width = 9, height = 4, dpi = 150)

## ---- Volcano plot: all hepatocyte genes, PT vs PVTT ----
volcano_df <- de_hep_pt_vs_pvtt %>%
  rownames_to_column("gene") %>%
  mutate(
    significant = p_val_adj < 0.05 & abs(avg_log2FC) > 1,
    label       = ifelse(significant & rank(p_val_adj) <= 15, gene, "")
  )

p_volcano <- ggplot(volcano_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_text_repel(aes(label = label), size = 3.2, max.overlaps = 20,
                   color = "black", fontface = "italic") +
  scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(
    title    = "Differential Expression: Hepatocytes, Primary Tumor vs PVTT",
    subtitle = "Genes labeled: top 15 by adjusted p-value",
    x        = expression(log[2]~"Fold Change (PT vs PVTT)"),
    y        = expression(-log[10]~"adjusted p-value")
  ) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(color = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave("figures/08_volcano_hepatocyte_PT_vs_PVTT.png", p_volcano, width = 8, height = 6, dpi = 150)


## =================================================================
## END OF ANALYSIS
## See results/ for all CSV tables and figures/ for all plots
## (numbered 01-08 in logical viewing order).
## See RESULTS_SUMMARY.md for the write-up answering the research
## question based on these outputs.
## =================================================================
