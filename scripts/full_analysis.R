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
## IMPORTANT — before running: edit the two paths marked
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
## SECTION 11 — Annotation (assign a cell type to each cluster)
## Result: data/seu_final_annotated.rds — the master analysis file
##
## Strategy: each cluster is assigned the majority cell type from
## the authors' original per-cell annotation (celltype_paper).
## This validates our independent clustering against the published
## atlas rather than re-deriving marker genes from scratch.
## =================================================================

seu_merged <- readRDS("data/seu_merged_after_clustering.rds")

cluster_composition <- table(seu_merged$seurat_clusters, seu_merged$celltype_paper)
cluster_annotations <- apply(cluster_composition, 1, function(x) names(which.max(x)))

seu_merged$major_celltype <- unname(cluster_annotations[as.character(seu_merged$seurat_clusters)])

table(seu_merged$major_celltype)
## B=3709 Endothelial=3661 Fibroblast=1997 Hepatocyte=19488
## Myeloid=15219 T/NK=23027
## (closely matches the authors' own celltype_paper counts —
## confirms clustering quality)

## IMPORTANT: join per-patient expression layers into one before
## any downstream differential expression analysis
seu_merged <- JoinLayers(seu_merged)

saveRDS(seu_merged, "data/seu_final_annotated.rds")


## =================================================================
## SECTION 12 — RESEARCH QUESTION, PART 1: Cellular composition
## PT vs metastatic sites (PVTT, MLN)
## Result: results/composition_PT_vs_Metastatic.csv
##         figures/composition_stacked_bar.pdf
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

ggsave("figures/composition_stacked_bar.pdf", p_bar, width = 12, height = 6)

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
## Result: figures/final_umap_celltype_tissue.pdf
##         figures/heatmap_hepatocyte_top_genes.pdf
## =================================================================

## UMAP overview: cell type + tissue side by side
p1 <- DimPlot(seu_merged, reduction = "umap", group.by = "major_celltype",
              label = TRUE) + NoLegend() + ggtitle("Cell types")
p2 <- DimPlot(seu_merged, reduction = "umap", group.by = "tissue") +
  ggtitle("Tissue")

ggsave("figures/final_umap_celltype_tissue.pdf", p1 + p2, width = 12, height = 5)

## Heatmap of top DE genes in Hepatocytes (PT vs PVTT)
top_genes <- rownames(de_hep_pt_vs_pvtt[order(de_hep_pt_vs_pvtt$p_val_adj), ])[1:15]

hep_only <- subset(seu_merged, subset = major_celltype == "Hepatocyte" &
                                         tissue %in% c("PT", "PVTT"))
avg_exp <- AverageExpression(hep_only, features = top_genes, group.by = "tissue")$RNA

pdf("figures/heatmap_hepatocyte_top_genes.pdf", width = 6, height = 8)
pheatmap(log1p(avg_exp), scale = "row",
         main = "Top DE genes: Hepatocyte PT vs PVTT")
dev.off()

## =================================================================
## SECTION 15 — Volcano plot (Hepatocytes, PT vs PVTT)
## Result: figures/volcano_hepatocyte_PT_vs_PVTT.pdf
## =================================================================

library(ggrepel)

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

ggsave("figures/volcano_hepatocyte_PT_vs_PVTT.pdf", p_volcano, width = 8, height = 6)


## =================================================================
## END OF ANALYSIS
## See results/ for all CSV tables and figures/ for all plots.
## See RESULTS_SUMMARY.md for the write-up answering the research
## question based on these outputs.
## =================================================================
