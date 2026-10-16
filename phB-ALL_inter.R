# ============================================================
# Single-cell RNA-seq analysis pipeline for PH/CML blast crisis samples
# ============================================================
# This script performs a comprehensive single-cell analysis, including:
#   1. Data loading and Seurat object creation per sample
#   2. Quality control (QC) and filtering
#   3. Normalization, variable feature selection
#   4. Scaling and PCA
#   5. Multi-sample integration using Harmony
#   6. Clustering, UMAP visualization
#   7. Marker gene identification
#   8. Cell type prediction
#   9. CytoTRACE2 analysis
#   10. Visualizing lineage marker expression
#   11. Manual cell type annotation
#   12. Differential expression between different disease
#   13. Cell cycle scoring
#   14. Signature scoring with AddModuleScore
# ============================================================

# --- set work directory ---
me=system("whoami", intern = T)
setwd(paste0("/public/", me, "/projects/CMLblast/CMLharmony/"))


# Load required libraries
library(data.table)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(Matrix)
library(patchwork)
library(RColorBrewer)
library(gplots)
library(magrittr)

# Set working directory
setwd(paste0("/public/", me, "/projects/CMLblast/PHall/"))

# ============================================================
# Step 1: Create Seurat objects from 10X data
# ============================================================
# Read 10X data for each sample and create Seurat objects.
# Barcodes are modified to include a sample-specific suffix to avoid collisions.
# min.cells = 10: keep genes detected in at least 10 cells
# min.features = 200: keep cells with at least 200 genes detected

me=system("whoami", intern = T)

## case1: BLBC: B-Lymphoblastic Blast Crisis;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case1_ZL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-1")
P1_BLBC <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10,
                           min.features = 200,
                           project = "P1_BLBC")

## case2: BLBC: B-Lymphoblastic Blast Crisis;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case2_SL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-2")
P2_BLBC <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10,
                           min.features = 200,
                           project = "P2_BLBC")

## case3: BLBC: B-Lymphoblastic Blast Crisis;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case3_CL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-3")
P3_BLBC <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10,
                           min.features = 200,
                           project = "P3_BLBC")



# PH_ALL1: PH+ B-ALL patient 1 (initial diagnosis)
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/PH001/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x), 1, 16), "-4")
PH_ALL1 <- CreateSeuratObject(counts = sce.10x,
                            min.cells = 10,
                            min.features = 200,
                            project = "PH_ALL1")

# PH_ALL2: PH+ B-ALL patient 2 (initial diagnosis)
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/PH002/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x), 1, 16), "-5")
PH_ALL2 <- CreateSeuratObject(counts = sce.10x,
                            min.cells = 10,
                            min.features = 200,
                            project = "PH_ALL2")


rm(sce.10x)



# ============================================================
# Step 2: Quality control (QC) and filtering
# ============================================================

# Compute mitochondrial, ribosomal, and hemoglobin gene percentages for each sample
# Helper function to add QC metrics
add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB_m <- match(HB.genes, rownames(obj@assays$RNA))
  HB.genes <- rownames(obj@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  obj[["percent.HB"]] <- PercentageFeatureSet(obj, features = HB.genes)
  return(obj)
}

P1_BLBC <- add_qc_metrics(P1_BLBC)
P2_BLBC <- add_qc_metrics(P2_BLBC)
P3_BLBC <- add_qc_metrics(P3_BLBC)
PH_ALL1 <- add_qc_metrics(PH_ALL1)
PH_ALL2 <- add_qc_metrics(PH_ALL2)


# Merge all samples for QC visualization
PHall <- merge(P1_BLBC, y = c(P2_BLBC, P3_BLBC, PH_ALL1, PH_ALL2))

# Plot QC metrics before filtering
p <- VlnPlot(PHall, group.by = "orig.ident",
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.HB"),
             log = TRUE, pt.size = 0.01)
ggsave(plot = p, file = "PHall_QC_before.pdf", width = 12, height = 12)

rm(PHall)

# Filter cells based on QC thresholds
# nFeature_RNA > 200 & < 6000: remove low-quality cells and doublets
# percent.mt < 10 (or <20 for P2_BLBC): remove dying cells
# percent.ribo < 54 (or <44.7, <52, <42 for some samples): remove cells with extreme ribosomal content
# percent.HB < 20: remove erythrocyte-contaminated cells (not applied here, but could be)

P1_BLBC <- subset(P1_BLBC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
P2_BLBC <- subset(P2_BLBC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
P3_BLBC <- subset(P3_BLBC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
PH_ALL1 <- subset(PH_ALL1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & percent.ribo < 54)
PH_ALL2 <- subset(PH_ALL2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & percent.ribo < 44.7)

# Merge all samples again for downstream analysis
PHall <- merge(P1_BLBC, y = c(P2_BLBC, P3_BLBC, PH_ALL1, PH_ALL2))

# Plot QC metrics after filtering
p <- VlnPlot(PHall, group.by = "orig.ident",
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.HB"),
             log = TRUE, pt.size = 0.01)
ggsave(plot = p, file = "PHall_QC_after.pdf", width = 12, height = 12)

# Create a subset excluding healthy donors for some analyses
PH <- merge(P1_BLBC, y = c(P2_BLBC, P3_BLBC, PH_ALL1, PH_ALL2))

# ============================================================
# Step 3: Normalization and variable feature selection
# ============================================================
# LogNormalize: divide by total counts per cell, multiply by scale factor (10000), then log1p
PH <- NormalizeData(object = PH, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable genes (HVGs) using variance stabilizing transformation (vst)
PH <- FindVariableFeatures(object = PH, selection.method = "vst", nfeatures = 3000)

# Identify top 10 variable genes
top10 <- head(VariableFeatures(PH), 10)
top10

# ============================================================
# Step 4: Scaling and PCA
# ============================================================
# Scale data: mean = 0, variance = 1 for each gene
# vars.to.regress = "percent.mt": regress out mitochondrial percentage
# Note: this step is memory-intensive; may be skipped if resources are limited

all.genes <- rownames(PH)
vars.to.regress = c("percent.mt")
PH <- ScaleData(PH, features = all.genes, vars.to.regress = vars.to.regress)

# Run PCA on scaled data
PH <- RunPCA(PH, features = VariableFeatures(object = PH), npcs = 50, verbose = FALSE)

# Examine PCA results
print(PH[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(PH, dims = 1:2, reduction = "pca")
DimPlot(PH, reduction = "pca") + NoLegend()

# Heatmap of PCA
p <- DimHeatmap(PH, dims = 1:15, cells = 500, balanced = TRUE)
ggsave(plot = p, file = "PH.DimHeatmap.pdf")

# Elbow plot to determine number of PCs
p <- ElbowPlot(PH, ndims = 50)
ggsave(plot = p, file = "PH.ElbowPlot.pdf")

# ============================================================
# Step 5: Multi-sample integration (batch correction)
# ============================================================
# Use Harmony integration to correct batch effects between samples
# orig.reduction = "pca": use PCA as input
# new.reduction = "harmony": name of the integrated reduction
PH <- IntegrateLayers(
  object = PH, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# Check available reductions
names(PH@reductions)

# Join layers for downstream analysis
PH <- JoinLayers(PH)

# Add sample and type metadata
# Sample_group: assign sample names based on cell numbers
Sample_group <- c(rep("P1_BLBC", 7220), rep("P2_BLBC", 8897), rep("P3_BLBC", 11013),
                  rep("PH_ALL1", 2294), rep("PH_ALL2", 3216))
Sample_group <- factor(Sample_group, levels = c("PH_ALL1", "PH_ALL2", "P1_BLBC", "P2_BLBC", "P3_BLBC"))
PH@meta.data$Sample <- Sample_group

# Type_group: assign disease types (CML-BLBC vs B-ALL)
Type_group <- c(rep("CML-BLBC", 7220), rep("CML-BLBC", 8897), rep("CML-BLBC", 11013),
                rep("B-ALL", 2294), rep("B-ALL", 3216))
Type_group <- factor(Type_group, levels = c("B-ALL", "CML-BLBC"))
PH@meta.data$Type <- Type_group

# Save metadata
write.csv(PH@meta.data, file = "PH_metaData.csv", row.names = TRUE)

# ============================================================
# Step 6: Cell clustering
# ============================================================
# Find neighbors using Harmony reduction
PH <- FindNeighbors(PH, reduction = "harmony", dims = 1:20)

# Cluster cells at multiple resolutions
PH <- FindClusters(PH, resolution = 1.2, cluster.name = "harmony_clusters12")
PH <- FindClusters(PH, resolution = 1.1, cluster.name = "harmony_clusters11")
PH <- FindClusters(PH, resolution = 1, cluster.name = "harmony_clusters10")
PH <- FindClusters(PH, resolution = 0.9, cluster.name = "harmony_clusters9")
PH <- FindClusters(PH, resolution = 0.8, cluster.name = "harmony_clusters8")
PH <- FindClusters(PH, resolution = 0.7, cluster.name = "harmony_clusters7")
PH <- FindClusters(PH, resolution = 0.6, cluster.name = "harmony_clusters6")
PH <- FindClusters(PH, resolution = 0.5, cluster.name = "harmony_clusters5")
PH <- FindClusters(PH, resolution = 0.4, cluster.name = "harmony_clusters4")
PH <- FindClusters(PH, resolution = 0.3, cluster.name = "harmony_clusters3")
PH <- FindClusters(PH, resolution = 0.2, cluster.name = "harmony_clusters2")
PH <- FindClusters(PH, resolution = 0.1, cluster.name = "harmony_clusters1")

# Visualize clustering tree
library(clustree)
p <- clustree(PH, prefix = "harmony_clusters")
ggsave(plot = p, file = "PH.clustree_harmony.pdf", height = 28, width = 14)

# Run UMAP on Harmony reduction
PH <- RunUMAP(PH, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

# Save UMAP plots by Sample and Type
pdf("UMAP_plot_integrated.harmony_by_Sample.pdf", height = 5, width = 6)
p1 <- DimPlot(PH, pt.size = 0.5, reduction = 'umap.harmony',
              cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Sample))),
              group.by = "Sample")
print(p1)
dev.off()

pdf("UMAP_plot_integrated.harmony_by_Samplesplit.pdf", height = 6, width = 18)
p1 <- DimPlot(PH, split.by = "Sample", pt.size = 0.5, reduction = 'umap.harmony',
              cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Sample))),
              group.by = "Sample")
print(p1)
dev.off()

pdf("UMAP_plot_integrated.harmony_by_Type.pdf", height = 5, width = 6)
p1 <- DimPlot(PH, pt.size = 0.5, reduction = 'umap.harmony',
              cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Type))),
              group.by = "Type")
print(p1)
dev.off()

pdf("UMAP_plot_integrated.harmony_by_splitType.pdf", height = 4, width = 10)
p1 <- DimPlot(PH, split.by = "Type", pt.size = 0.5, reduction = 'umap.harmony',
              cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Type))),
              group.by = "Type")
print(p1)
dev.off()


# Bar plots showing cluster distribution by Sample
p <- ggplot(PH@meta.data, aes(x = harmony_clusters12))
p <- p + geom_bar(aes(fill = Sample), position = "fill")
ggsave("cluster.harmony_SampledisperAllCell_1.2.pdf", device = "pdf")

# ... (similar bar plots for other resolutions)



# ============================================================
# Step 7: Identify cell markers between clusters
# ============================================================
# Set active identity to cluster
Idents(PH) <- "harmony_clusters8"

# Rename cluster labels to start from 1
cluster_add1_idents <- as.integer(PH@active.ident)
cluster_add1_idents <- factor(cluster_add1_idents, levels = 1:length(unique(cluster_add1_idents)))
names(cluster_add1_idents) <- rownames(PH@meta.data)
PH@active.ident <- cluster_add1_idents
PH@meta.data$cluster <- cluster_add1_idents

# Plot UMAP colored by cluster
pdf("UMAP_plot_by_harmonyCluster.pdf", height = 8, width = 10)
p1 <- DimPlot(PH, cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$cluster))),
              pt.size = 0.5, label = TRUE, reduction = "umap.harmony",
              group.by = "cluster", label.size = 6)
print(p1)
dev.off()


pdf("UMAP_plot_by_harmonyCluster_splitSample.pdf",width = 24)
p1 <- DimPlot(PH, cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$cluster))),
              pt.size = 0.5, label = TRUE, reduction = "umap.harmony",
              group.by = "cluster",label.size = 4, split.by = "Sample")
print(p1)
dev.off()


pdf("UMAP_plot_by_harmonyCluster_splitType.pdf",width = 16)
p1 <- DimPlot(PH, cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$cluster))),
              pt.size = 0.5, label = TRUE, reduction = "umap.harmony",
              group.by = "cluster",label.size = 4, split.by = "Type")
print(p1)
dev.off()


# Find all markers for each cluster using Wilcoxon test
DefaultAssay(PH) <- "RNA"
var_genes <- VariableFeatures(PH)
library(presto)
diff_wilcox <- FindAllMarkers(PH, features = var_genes, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1)

# For each cluster, plot top 12 markers as FeaturePlot and VlnPlot
for (a_cluster in levels(PH@active.ident)) {
  mks <- diff_wilcox[diff_wilcox$cluster == a_cluster, ]
  mks1 <- mks[order(mks$avg_log2FC, decreasing = TRUE), "gene"][1:12]
  # FeaturePlot
  otfile <- paste("harmonyCluster", a_cluster, "_wilcox_markers_umap.pdf", sep = "")
  pdf(file = otfile, width = 16, height = 10)
  p <- FeaturePlot(PH, features = mks1, cols = c("grey", "red"), reduction = "umap.harmony")
  print(p)
  dev.off()
  # VlnPlot
  outfile <- paste("harmonyCluster", a_cluster, "_VlnPlot_wilcox_markers.pdf", sep = "")
  pdf(file = outfile, width = 8, height = 20)
  p1 <- VlnPlot(PH, features = mks1,
                cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$cluster))),
                pt.size = 0, ncol = 1)
  print(p1)
  dev.off()
}


# Write differential expression results to files
for (a_cluster in levels(PH@active.ident)) {
  subdiff <- diff_wilcox[diff_wilcox$cluster == a_cluster, ]
  otfile <- paste("harmonyCluster", a_cluster, "_diff_wilcox.xls", sep = "")
  write.table(subdiff, sep = "\t", file = otfile, row.names = FALSE, quote = FALSE)
  otfile <- paste("harmonyCluster", a_cluster, "_diff_wilcox_significant.xls", sep = "")
  write.table(subdiff[subdiff[, 1] <= 0.05 & abs(subdiff[, 2]) >= 0.585, ],
              sep = "\t", file = otfile, row.names = FALSE, quote = FALSE)
}



# ============================================================
# Step 8: Cell type prediction
# ============================================================
# scHCL algorithm for automatic cell type annotation
library(scHCL)

# Compute average expression per cluster
PH.AVG <- AverageExpression(object = PH, assays = "RNA", layer = "counts",
                            group.by = "cluster", return.seurat = TRUE)
cluster_mean_counts <- GetAssayData(PH.AVG[["RNA"]], layer = "counts")
cluster_mean_counts_matrix <- as.matrix(cluster_mean_counts)

# Run scHCL
ct_cluster <- scHCL(scdata = cluster_mean_counts_matrix, numbers_plot = 10)
write.csv(ct_cluster$scHCL_probility, file = "PH.AVG_scHCL_ByCluster.csv")



# ============================================================
# Step 9: CytoTRACE2 analysis
# ============================================================
library(CytoTRACE2)

# Run CytoTRACE2 on the Seurat object
cytotrace2_result <- cytotrace2(PH, is_seurat = TRUE, slot_type = "counts", species = 'human')

# Inspect results
head(cytotrace2_result, 2)

# Visualize CytoTRACE2 results
annotation <- data.frame(phenotype = cytotrace2_result@meta.data$cluster) %>%
  set_rownames(., colnames(cytotrace2_result))
plots <- plotData(cytotrace2_result = cytotrace2_result, annotation = annotation, is_seurat = TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
p <- (p1 + p2 + p3 + p4) + plot_layout(ncol = 2)
ggsave(plot = p, "CytoTRACE2_UMAP_plotData.pdf")

# UMAP by CytoTRACE2_Potency split by Sample
p <- DimPlot(cytotrace2_result, reduction = "umap.harmony", split.by = "Sample",
             cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(cytotrace2_result@meta.data$CytoTRACE2_Potency))),
             label = FALSE, group.by = "CytoTRACE2_Potency") + theme(aspect.ratio = 1)
ggsave(plot = p, "CytoTRACE2_UMAP_Potency_splitSample.pdf", width = 20, height = 5)

# FeaturePlot for CytoTRACE2_Relative
p <- FeaturePlot(cytotrace2_result, features = "CytoTRACE2_Relative",
                 cols = viridis::viridis(100), order = TRUE, reduction = "umap.harmony")
ggsave(plot = p, "CytoTRACE2_UMAP_Relative.pdf", width = 6, height = 5)

# Violin plot for CytoTRACE2_Score by cluster
p <- VlnPlot(cytotrace2_result, group.by = "cluster", features = "CytoTRACE2_Score",
             pt.size = 0, ncol = 1, layer = "data",
             cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(cytotrace2_result@meta.data$cluster))))
ggsave(plot = p, file = "vln_CytoTRACE2_score.pdf", width = 6, height = 2)

# Bar plots for CytoTRACE2_Potency distribution
p <- ggplot(cytotrace2_result@meta.data, aes(x = cluster))
p <- p + geom_bar(aes(fill = CytoTRACE2_Potency), position = "fill")
ggsave("cluster.harmony_PotencydisperAllCell.pdf", device = "pdf")

p <- ggplot(cytotrace2_result@meta.data, aes(x = Sample))
p <- p + geom_bar(aes(fill = CytoTRACE2_Potency), position = "fill")
ggsave("sample.harmony_PotencydisperAllCell.pdf", device = "pdf")



# ============================================================
# Step 10: Visualizing lineage marker expression
# ============================================================

Idents(PH) <- "cluster"

# Define marker gene lists
Cycling_marker <- c("CD34", "MCM7", "TYMS", "PCNA", "GINS2", "MKI67", "TOP2A", "ASPM", "PCLAF", "CCNB1")
HSPC_marker <- c("CD34", "FLT3", "SOX4", "AVP", "SPINK2", "LAPTM4B", "CAVIN2", "NPR3", "PRSS57",
                 "AIF1", "CRHBP", "HTR1F", "MYCT1", "MSRB3", "ANGPT1", "FAM30A", "GATA2")
ErythrocyticCell_marker <- c("HBA1", "HBA2", "HBB", "HBG2", "AHSP", "CA1", "GYPA", "HBG1",
                             "SLC4A1", "APOC1", "KLF1", "HBD", "FAM178B", "RHCE", "BLVRB")
progenitor_marker <- c("CD34", "PRSS57", "DNTT", "SMIM3", "MT1X", "SOCS2", "HPS4", "MLXIP",
                       "VPREB1", "TSC22D", "CEBPA", "APOE")
GMP_marker <- c("MPO", "ELANE", "CTSG", "AZU1", "LYST", "LYZ", "CEBPD", "MNDA")
MEP_marker <- c("GATA1", "KLF1", "APOE", "AHSP", "CA1")
MyloidCell_marker <- c("CD68", "FCGR3B", "CD163", "ITGAM", "CD1C", "MPO", "ITGB2", "CD244",
                       "TNFRSF8", "CD33", "ANPEP")
Monocyte_marker <- c("S100A8", "MNDA", "CD68", "FCN1", "CD14", "S100A12", "FCGR3A", "CSF1R",
                     "RNASE2", "UGCG")
ClaMOcell_marker <- c("S100A8", "S100A12", "CD14", "LYZ", "CD68", "RETN", "CFD", "SERPINA1",
                      "LGALS2", "CEBPB", "VCAN")
NonClaMOcell_marker <- c("FCGR3A", "AIF1", "MS4A7", "LST1", "FTH1", "C1QA", "MAFB", "CDKN1C",
                         "HES4", "HMOX1", "CSF1R", "LYPD2")
pDC_marker <- c("CD1C", "FCER1A", "IL3RA", "UGCG", "LILRA4", "CLEC4C", "SMPD3", "SERPINF1",
                "MS4A7", "RNASE2", "PPBP", "PF4", "TUBB1", "GP9")
MacroCHECK_marker <- c("CALR", "CD47", "SIRPA", "CD24", "SIGLEC10", "VSIG4", "VCAM1", "SLAMF7")
Tcell_marker <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD4", "CD2", "CD5", "CD7", "CD8B", "CD27", "CD70")
naiveT_marker <- c("CD3D", "CD3E", "CD3G", "CD4", "CD27", "LEF1", "IL7R", "IL32")
effectT_marker <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD70", "CD27", "CCL4", "CCL5",
                    "GZMK", "GZMA", "CST7", "NKG7")
NKcell_marker <- c("GNLY", "NKG7", "CTSW", "FGFBP2", "KLRD1", "KLRB1", "PRF1", "FCER1G",
                   "KLRC1", "KLRF1", "GZMB", "TRDC", "CD247")
NKT_marker <- c("CD3D", "TRBC1", "CD8A", "CCR7", "GPR183", "GZMA", "GZMB", "GZMH", "KLRB1",
                "GNLY", "SPON2", "CCL3", "CCL5", "HSPA6")
NKTCHECK_marker <- c("PRF1", "GZMA", "GZMB", "GZMH", "CST7", "KLRB1", "KLRD1", "KLRG1",
                     "KLRC1", "KLRC2", "PDCD1", "HAVCR2", "CTLA4", "LAG3", "TIGIT")
TCHECK_marker <- c("CD274", "PDCD1", "CTLA4", "CD80", "CD86", "CD28", "LAG3", "HAVCR2", "TIGIT",
                   "TNFRSF9", "TOX", "TOX2", "CMTM4", "SATB1", "VSIG4", "TGIF1", "VCAM1")
NKCHECK_marker <- c("KLRB1", "CLEC2D", "KLRG1", "KLRC2", "LILRB1", "CD96", "PVR", "IL18R1",
                    "IL18RAP", "KLRD1")
THCHECK_marker <- c("NCR3", "IL7R", "CD247", "CCR7", "KIR2DL3", "IL17F", "IL23R", "RORC",
                    "CD96", "FOXP3")
Cyto_features <- c('GNLY', 'NKG7', 'PRF1', 'GZMB', 'GZMA', 'GZMH', 'IFNG', 'CST7', 'TNFSF10')
B_marker <- c("CD24", "EBF1", "MME", "VPREB1", "PAX5", "CD19", "CD79A", "MS4A1", "BANK1",
              "CD27", "CD70", "MZB1", "IGLL5", "JCHAIN")
naiveB_marker <- c("IGKC", "IGHD", "IGHM", "IGLL1", "CD83", "VPREB3", "LHPP", "NSMCE", "ADAM28")
MemoryBCell_marker <- c("CD79A", "MS4A1", "IFI44L", "CD27", "CD70", "GPR183", "EBI3", "IGLL1",
                        "IGHG1", "CD81")
Plasmablasts_marker <- c("IGHG2", "IGHG3", "IGHG4", "IGLC3", "IGLV1", "JCHAIN", "CD38",
                         "ZBP1", "XBP1", "MZB1", "HSPA5", "MANF", "CALR")
pro_BT_marker <- c("CD34", "DNTT", "SOCS2", "HPS4", "VPREB1", "CD24", "MME", "CD79A", "MS4A1",
                   "CD3D", "MPO", "ANPEP")
AllCell_marker <- c("CD34", "FLT3", "SOX4", "GATA1", "KLF1", "HBB", "PF4", "MPO", "LYZ",
                    "VCAN", "S100A8", "CD68", "CD163", "MAFB", "CD1C", "LILRA4", "CD79A",
                    "MS4A1", "PAX5", "CD27", "CD70", "MZB1", "LEF1", "CD3E", "CD8A", "CCL5",
                    "KLRB1", "CD247", "PTPRC")
AllCell_marker1 <- c("CRHBP", "AVP", "CD34", "FLT3", "KIT", "GATA2", "GATA1", "ITGA2B", "KLF1",
                     "TFRC", "CD36", "HBD", "HBA1", "CA1", "GYPA", "ALAS2", "VWF", "GP9", "PF4",
                     "CSF3R", "MPO", "ELANE", "LYZ", "IRF8", "CD14", "MME", "VPREB1", "CD79A",
                     "IGHM", "IGHD", "EBF1", "VPREB3", "CD24", "PAX5", "E2F2", "CD19", "MS4A1",
                     "MKI67", "PLK1")
AllCell_marker2 <- c("AVP", "CD34", "FLT3", "GATA2", "GATA1", "KLF1", "CD36", "HBD", "CA1",
                     "GP9", "MPO", "LYZ", "MME", "CD79A", "CD24", "EBF1", "PAX5", "CD19",
                     "MS4A1", "MKI67", "CD3D")

# Function to generate VlnPlot and FeaturePlot for a marker set
plot_markers <- function(obj, markers, prefix, width_vln = 8, height_vln = 24, width_feat = 16, height_feat = 10) {
  p <- VlnPlot(obj, group.by = "cluster", features = markers, pt.size = 0, ncol = 1, layer = "data",
               cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(obj@meta.data$cluster))))
  ggsave(plot = p, file = paste0(prefix, "_VlnPlot.pdf"), width = width_vln, height = height_vln)
  p <- FeaturePlot(obj, features = markers, cols = c("grey", "red"), reduction = "umap.harmony")
  ggsave(plot = p, file = paste0(prefix, "_FeaturePlot.pdf"), width = width_feat, height = height_feat)
}

# Apply to each marker set
plot_markers(PH, Cycling_marker, "vlnPH_resolution0.9_Cycling_marker_Data")
plot_markers(PH, HSPC_marker, "vlnPH_resolution0.9_HSPC_marker_Data")
plot_markers(PH, ErythrocyticCell_marker, "vlnPH_resolution0.9_ErythrocyticCell_marker_Data")
plot_markers(PH, progenitor_marker, "vlnPH_resolution0.9_progenitor_marker_Data")
plot_markers(PH, GMP_marker, "vlnPH_resolution0.9_GMP_marker_Data")
plot_markers(PH, MEP_marker, "vlnPH_resolution0.9_MEP_marker_Data")
plot_markers(PH, MyloidCell_marker, "vlnPH_resolution0.9_MyloidCell_marker_Data")
plot_markers(PH, Monocyte_marker, "vlnPH_resolution0.9_Monocyte_marker_Data")
plot_markers(PH, ClaMOcell_marker, "vlnPH_resolution0.9_ClaMOcell_marker_Data")
plot_markers(PH, NonClaMOcell_marker, "vlnPH_resolution0.9_NonClaMOcell_marker_Data")
plot_markers(PH, pDC_marker, "vlnPH_resolution0.9_pDC_marker_Data")
plot_markers(PH, MacroCHECK_marker, "vlnPH_resolution0.9_MacroCHECK_marker_Data")
plot_markers(PH, Tcell_marker, "vlnPH_resolution0.9_Tcell_marker_Data")
plot_markers(PH, naiveT_marker, "vlnPH_resolution0.9_naiveT_marker_Data")
plot_markers(PH, effectT_marker, "vlnPH_resolution0.9_effectT_marker_Data")
plot_markers(PH, NKcell_marker, "vlnPH_resolution0.9_NKcell_marker_Data")
plot_markers(PH, NKT_marker, "vlnPH_resolution0.9_NKT_marker_Data")
plot_markers(PH, NKTCHECK_marker, "vlnPH_resolution0.9_NKTCHECK_marker_Data")
plot_markers(PH, TCHECK_marker, "vlnPH_resolution0.9_TCHECK_marker_Data")
plot_markers(PH, NKCHECK_marker, "vlnPH_resolution0.9_NKCHECK_marker_Data")
plot_markers(PH, THCHECK_marker, "vlnPH_resolution0.9_THCHECK_marker_Data")
plot_markers(PH, Cyto_features, "vlnPH_resolution0.9_Cyto_marker_Data")
plot_markers(PH, B_marker, "vlnPH_resolution0.9_B_marker_Data")
plot_markers(PH, naiveB_marker, "vlnPH_resolution0.9_naiveB_marker_Data")
plot_markers(PH, MemoryBCell_marker, "vlnPH_resolution0.9_MemoryBCell_marker_Data")
plot_markers(PH, Plasmablasts_marker, "vlnPH_resolution0.9_Plasmablasts_marker_Data")
plot_markers(PH, pro_BT_marker, "vlnPH_resolution0.9_pro_BT_marker_Data")
plot_markers(PH, AllCell_marker, "dotplotPH_resolution0.9_AllCell_marker_Data")
plot_markers(PH, AllCell_marker1, "dotplotPH_resolution0.9_AllCell_marker1_Data")
plot_markers(PH, AllCell_marker2, "dotplotPH_resolution0.9_AllCell_marker2_Data")



# ============================================================
# Step 11: Manual cell type annotation
# ============================================================

# Map cluster IDs to major cell types
identify_celltypes <- c("B", "B", "B", "B", "B", "B", "B", "B", "B", "T_NK", "T_NK",
                        "B", "B", "HSPC", "ProMono", "B", "B", "EEP")
names(identify_celltypes) <- 1:length(identify_celltypes)
clusters <- PH@meta.data$cluster
celltypes <- identify_celltypes[clusters]
celltypes <- factor(celltypes, level = c("HSPC", "EEP", "ProMono", "B", "T_NK"))
PH@meta.data$majorCelltype <- celltypes

# More detailed cell types
identify_celltypes <- c("proB", "preB", "preB", "preB", "preB", "proB", "preB", "preB",
                        "proB", "T_NK", "T_NK", "preB", "preB", "HSPC", "ProMono", "preB",
                        "B", "EEP")
names(identify_celltypes) <- 1:length(identify_celltypes)
clusters <- PH@meta.data$cluster
celltypes <- identify_celltypes[clusters]
PH@meta.data$Celltype <- celltypes
PH@meta.data$Celltype <- factor(PH@meta.data$Celltype,
                                level = c("HSPC", "EEP", "ProMono", "proB", "preB", "B", "T_NK"))

# Plot UMAP by Celltype
pdf("UMAP_plot_resolution0.9_by_Celltype.pdf", width = 8, height = 6)
p1 <- DimPlot(PH, cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Celltype))),
              pt.size = 0.5, reduction = "umap.harmony", group.by = "Celltype")
print(p1)
dev.off()

# Plot UMAP by Celltype split by Type
pdf("UMAP_plot_resolution0.9_by_Celltype_splitType.pdf", width=10, height=6)
p1 <- DimPlot(PH, cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Celltype))),
              pt.size = 0.5, reduction = "umap.harmony", group.by = "Celltype", split.by = "Type")
print(p1)
dev.off()

# Plot UMAP by Celltype split by Sample
pdf("UMAP_plot_resolution0.9_by_Celltype_splitSample.pdf", width=28, height=5)
p1 <- DimPlot(PH, cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Celltype))),
              pt.size = 0.5, reduction = "umap.harmony", group.by = "Celltype", split.by = "Sample")
print(p1)
dev.off()



# Set active identity to Celltype
celltypes <- PH@meta.data$Celltype
names(celltypes) <- rownames(PH@meta.data)
PH@active.ident <- celltypes

# DotPlot for selected markers by Celltype
p <- DotPlot(PH, features = AllCell_marker3, cols = c("grey", "#498EA4"))
ggsave(plot = p, file = "dotplotPH_resolution0.9_AllCell_marker3_Data_byCelltype.pdf", width = 8, height = 3)



# ============================================================
# Step 12: Differential expression between CML-BLBC and Ph+ B-ALL
# ============================================================
DefaultAssay(PH) <- "RNA"

PH@meta.data$Type.majorCelltype <- paste(PH@meta.data$Type, PH@meta.data$majorCelltype, sep = "_")
Idents(PH) <- "Type.majorCelltype"
BCmarkers <- FindMarkers(PH, ident.1 = "CML-BLBC_B", ident.2 = "B-ALL_B", verbose = FALSE, test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)
write.csv(BCmarkers, file="BCmarkers_wilcox.csv", row.names=TRUE)



# ============================================================
# Step 13: Cell cycle scoring
# ============================================================
# Use CellCycleScoring to assign S and G2M scores
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(PH))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(PH))

PH <- CellCycleScoring(PH, g2m.features = g2m_genes, s.features = s_genes)

# Plot cell cycle scores
p1 <- VlnPlot(PH, "S.Score", pt.size = 0, group.by = "Phase")
p2 <- VlnPlot(PH, "G2M.Score", pt.size = 0, group.by = "Phase")
p1 + p2
ggsave("CellCycle_VlnPlot.pdf", width = 4, height = 4)

# Bar plots of cell cycle phase distribution
p <- ggplot(PH@meta.data, aes(x = cluster))
p <- p + geom_bar(aes(fill = factor(Phase)), position = "fill")
ggsave("cluster.harmony_Phase_resolution0.9_dodge.pdf", device = "pdf")

# Ridge plot for cell cycle markers
pdf("CellCycle_RidgePlot.pdf", width = 16, height = 10)
p <- RidgePlot(PH, cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(factor(PH@meta.data$Phase)))),
               group.by = "Phase", features = c("MKI67", "PCNA", "TOP2A", "HMGB2", "CENPF", "UBE2C",
                                                "TYMS", "NUSAP1", "CDC20", "TPX2", "GINS2", "UHRF1",
                                                "HELLS", "NASP"))
print(p)
dev.off()

# UMAP plots of cell cycle phase
pdf("UMAP_plot_by_Phase.pdf",width = 5, height = 4)
p <- DimPlot(PH, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(factor(PH@meta.data$Phase)))), 
  pt.size = 0.5, label = FALSE, reduction = "umap.harmony", group.by = "Phase", label.size = 4)
print(p)
dev.off()




# ============================================================
# Step 14: Signature scoring with AddModuleScore
# ============================================================
all.genes <- rownames(PH)


# Load gene signatures from files
# --- Apoptosis signature ---
Apoptosis = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/Apoptosis.txt"),sep='\t',header=T,stringsAsFactors =F)
# Keep only genes that are present in the Seurat object
Apoptosis <- Apoptosis[Apoptosis$gene %in% all.genes,]
# Convert to a list (AddModuleScore expects a list of gene sets)
Apoptosis <- as.list(Apoptosis)


# --- P53 pathway signature ---
P53 = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/P53_pathway.txt"),sep='\t',header=T,stringsAsFactors =F)
P53 <- P53[P53$gene %in% all.genes,]
P53 <- as.list(P53$gene)


# --- Proliferation signature ---
Proliferation = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/Proliferation.txt"),sep='\t',header=T,stringsAsFactors =F)
Proliferation <- Proliferation[Proliferation$gene %in% all.genes,]
Proliferation <- as.list(Proliferation$gene)


# --- Tumor inflammation signature ---
Tumor_inflammation = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/Tumor_inflammation.txt"),sep='\t',header=T,stringsAsFactors =F)
Tumor_inflammation <- Tumor_inflammation[Tumor_inflammation$gene %in% all.genes,]
Tumor_inflammation <- as.list(Tumor_inflammation)


# --- IL6 signature ---
IL6 = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/IL6.txt"),sep='\t',header=T,stringsAsFactors =F)
IL6 <- IL6[IL6$gene %in% all.genes,]
IL6 <- as.list(IL6)


# --- IL12 signature ---
IL12 = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/IL12.txt"),sep='\t',header=T,stringsAsFactors =F)
IL12 <- IL12[IL12$gene %in% all.genes,]
IL12 <- as.list(IL12$gene)


# --- BCR (B-cell receptor) signature ---
BCR = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/BCR.txt"),,sep='\t',header=T,stringsAsFactors =F)
BCR <- BCR[BCR$gene %in% all.genes,]
BCR <- as.list(BCR$gene)


# Add module scores
PH <- AddModuleScore(object = PH, features = Apoptosis, nbin = 25, ctrl = 20, name = 'Apoptosis')
PH <- AddModuleScore(object = PH, features = P53, nbin = 25, ctrl = 20, name = 'P53')
PH <- AddModuleScore(object = PH, features = Proliferation, nbin = 25, ctrl = 20, name = 'Proliferation')
PH <- AddModuleScore(object = PH, features = IL12, nbin = 25, ctrl = 20, name = 'IL12')
PH <- AddModuleScore(object = PH, features = BCR, nbin = 25, ctrl = 20, name = 'BCR')
PH <- AddModuleScore(object = PH, features = Tumor_inflammation, nbin = 25, ctrl = 20, name = 'Tumor_inflammation')
PH <- AddModuleScore(object = PH, features = IL6, nbin = 25, ctrl = 20, name = 'IL6_production')


# Aggregate multiple score columns into single scores
# --- BCR score ---
metaBCR <- PH@meta.data %>% select(matches("BCR"))
BCRscore <- apply(metaBCR, 1, mean)
PH@meta.data$BCRscores <- BCRscore

# --- IL12 score ---
metaIL2 <- PH@meta.data %>% select(matches("IL2"))
IL2Score <- apply(metaIL2, 1, mean)
PH@meta.data$IL2Scores <- IL2Score

# --- IL6 score ---
metaIL6 <- CML@meta.data %>% select(matches("IL6"))
IL6Score <- apply(metaIL6, 1, mean)
CML@meta.data$IL6Scores <- IL6Score

# --- Tumor inflammation score ---
metaTumor_inflammation <- CML@meta.data %>% select(matches("Tumor_inflammation"))
Tumor_inflammationScore <- apply(metaTumor_inflammation, 1, mean)
CML@meta.data$Tumor_inflammationScores <- Tumor_inflammationScore

# --- Proliferation score ---
metaProliferation <- CML@meta.data %>% select(matches("Proliferation"))
ProliferationScore <- apply(metaProliferation, 1, mean)
CML@meta.data$ProliferationScores <- ProliferationScore

# --- P53 score ---
metaP53 <- CML@meta.data %>% select(matches("P53"))
P53Score <- apply(metaP53, 1, mean)
CML@meta.data$P53Scores <- P53Score

# --- Apoptosis score ---
metaApoptosis <- CML@meta.data %>% select(matches("Apoptosis"))
ApoptosisScore <- apply(metaApoptosis, 1, mean)
CML@meta.data$ApoptosisScores <- ApoptosisScore


# Collect all signature scores
type <- PH@meta.data$Type
names(type) <- rownames(PH@meta.data)
PH@active.ident <- type

signatureScore1 <- colnames(PH@meta.data)[grep("Score", colnames(PH@meta.data))]
signatureScore2 <- colnames(PH@meta.data)[grep("score", colnames(PH@meta.data))]
signatureScores <- c(signatureScore1, signatureScore2)

# Violin plots for signature scores by Type
pdf("signatureScores_VlnPlot_byType_new.pdf", width = 8, height = 8)
p <- VlnPlot(PH, features = signatureScores, group.by = "Type",
             cols = colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(PH@meta.data$Type))),
             pt.size = 0, ncol = 5)
print(p)
dev.off()


# Boxplots for selected scores
mdata <- PH@meta.data
s1 <- mdata[mdata$Type == "PH", ]

pdf("boxplot_Tumor_inflammationScores_byType.pdf", width = 3, height = 5)
boxplot(Tumor_inflammationScores ~ Type, data = mdata, outline = FALSE, col = NA,
        border = c("Violet", "DodgerBlue"))
dev.off()

