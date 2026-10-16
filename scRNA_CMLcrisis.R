# ============================================================
# CML Blast Crisis scRNA-seq Analysis Pipeline
# ============================================================
# This script performs a comprehensive single-cell analysis, including:
#   1. Data loading and Seurat object creation per sample
#   2. Quality control (QC) and filtering
#   3. Normalization and variable feature selection
#   4. Scaling and PCA
#   5. Multi-sample integration using Harmony
#   6. Assign sample metadata
#   7. Clustering and UMAP visualization
#   8. Cell type annotation
#   9. Differential expression between different disease status
#   10. Gene signature scoring
#   11. Pseudotime trajectory analysis
#   12. Add BCR::ABL1 fusion to Seurat metadata
# ============================================================
# Author: Xuehong Zhang
# Description: Full pipeline for CP-BC CML single-cell analysis


# --- Load required packages ---
library(Seurat)                       # Core single-cell analysis
library(SeuratObject)            # Seurat object infrastructure
library(ggplot2)                     # Visualization
library(DoubleFinder)           # Doublet detection
library(dplyr)                         # Data manipulation
library(Matrix)                       # Sparse matrix operations
library(patchwork)                # Multi-panel figure assembly
library(RColorBrewer)          # Color palettes for heatmaps and grouped visualizations
library(gplots)                       # Heatmap and visualization
library(monocle2)                 # Pseudotime trajectory analysis
library(magrittr)                    # Pipe operators for readable workflows
library(scHCL)                      # Reference-based cell type annotation



#### Step 1: Create Seurat object 

######### Read in single-cell data ########

### Assume we have a Seurat object named sce.all. By default, the sample origin information for each cell is stored in the metadata slot orig.ident (sce.all@meta.data$orig.ident).
### However, the contents of orig.ident are numeric values starting from 1. If there are N samples, there will be N numbers representing these N samples. In downstream plotting and analysis, if sample information is only numeric, one must map it back to the original sample names to check, which is not intuitive.
### Therefore, we need to replace orig.ident in the Seurat object with the actual sample names.

# For 10X single-cell sequencing samples, the barcode sequence length is 16 bp.
### After reading 10X single-cell sequencing samples, meta.data$orig.ident already contains sample names.

# ---  Create Seurat object with minimal gene/cell thresholds ---
# min.cells = 10: keep genes detected in at least 10 cells
# min.features = 200: keep cells with at least 200 detected genes

me=system("whoami", intern = T)
## case1: CP: Chronic-Phase;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case1_ZBL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-1")
P1_CP <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10, 
                           min.features = 200,  
                           project = "P1_CP")

## case1: BLBC: B-Lymphoblastic Blast Crisis;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case1_ZL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-2")
P1_BLBC <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10,
                           min.features = 200,
                           project = "P1_CPBC")

## case2: CP: Chronic-Phase;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case2_SBL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-3")
P2_CP <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10,
                           min.features = 200,
                           project = "P2_CP")

## case2: BLBC: B-Lymphoblastic Blast Crisis;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case2_SL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-4")
P2_BLBC <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10,
                           min.features = 200,
                           project = "P2_CPBC")

## case3: BLBC: B-Lymphoblastic Blast Crisis;
sce.10x <- Read10X(data.dir = paste0("/public/", me, "/projects/CMLblast/cellranger/case3_CL/"))
colnames(sce.10x) <- paste0(substr(colnames(sce.10x),1,16),"-5")
P3_BLBC <- CreateSeuratObject(counts = sce.10x,
                           min.cells = 10,
                           min.features = 200,
                           project = "P3_BLBC")


rm(sce.10x)




# ============================================================
# Step 2: Quality control (QC) and filtering
# ============================================================

# ---- 2.1 Calculate mitochondrial gene percentage ---------
# Human mitochondrial genes start with "MT-" (mouse: "mt-")
# Low-quality/dying cells often show high mitochondrial expression

P1_CP[["percent.mt"]] <- PercentageFeatureSet(P1_CP, pattern = "^MT-")
P1_BLBC[["percent.mt"]] <- PercentageFeatureSet(P1_BLBC, pattern = "^MT-")
P2_CP[["percent.mt"]] <- PercentageFeatureSet(P2_CP, pattern = "^MT-")
P2_BLBC[["percent.mt"]] <- PercentageFeatureSet(P2_BLBC, pattern = "^MT-")
P3_BLBC[["percent.mt"]] <- PercentageFeatureSet(P3_BLBC, pattern = "^MT-")



# ============================================================
# ---- 2.2 Ribosomal gene percentage
# ============================================================
# Ribosomal protein genes are typically prefixed with "RPS" or "RPL"
# in human (e.g., RPS6, RPL13) and "Rps" or "Rpl" in mouse.
# High ribosomal content may indicate actively translating cells,
# but extremely high values can also reflect low-quality or stressed cells.

P1_CP[["percent.ribo"]] <- PercentageFeatureSet(P1_CP, pattern = "^RP[SL]")
P1_BLBC[["percent.ribo"]] <- PercentageFeatureSet(P1_BLBC, pattern = "^RP[SL]")
P2_CP[["percent.ribo"]] <- PercentageFeatureSet(P2_CP, pattern = "^RP[SL]")
P2_BLBC[["percent.ribo"]] <- PercentageFeatureSet(P2_BLBC, pattern = "^RP[SL]")
P3_BLBC[["percent.ribo"]] <- PercentageFeatureSet(P3_BLBC, pattern = "^RP[SL]")



# ============================================================
# ---- 2.3 Erythrocyte (hemoglobin) gene percentage
# ============================================================
# Erythrocyte contamination is common in blood or tissue dissociations.
# Human hemoglobin genes include HBA1, HBA2, HBB, HBD, HBE1,
# HBG1, HBG2, HBM, HBQ1, and HBZ.

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")

# 2.3.1 Find the indices of hemoglobin genes in the Seurat object's RNA assay.
#    `match()` returns the position of each gene in `rownames(P1_CP@assays$RNA)`.
#    Genes not present in the dataset will be assigned NA.
HB_m <- match(HB.genes, rownames(P1_CP@assays$RNA)) 

# 2.3.2 Retrieve the actual gene names that exist in the dataset using the indices.
#    This replaces the original `HB.genes` vector with only the genes found.
HB.genes <- rownames(P1_CP@assays$RNA)[HB_m] 

# 2.3.3 Remove NA values (i.e., genes that were not found in the dataset).
#    The final `HB.genes` contains only valid hemoglobin gene symbols.
HB.genes <- HB.genes[!is.na(HB.genes)] 

# 2.3.4 Calculate the percentage of counts mapping to the filtered hemoglobin genes
#    for each cell, and store the result in the Seurat object's metadata
#    under the column name `percent.HB`.
P1_CP[["percent.HB"]]<-PercentageFeatureSet(P1_CP, features=HB.genes) 

# P1_BLBC
HB_m <- match(HB.genes, rownames(P1_BLBC@assays$RNA)) 
HB.genes <- rownames(P1_BLBC@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
P1_BLBC[["percent.HB"]]<-PercentageFeatureSet(P1_BLBC, features=HB.genes) 

# P2_CP
HB_m <- match(HB.genes, rownames(P2_CP@assays$RNA)) 
HB.genes <- rownames(P2_CP@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
P2_CP[["percent.HB"]]<-PercentageFeatureSet(P2_CP, features=HB.genes) 

# P2_BLBC
HB_m <- match(HB.genes, rownames(P2_BLBC@assays$RNA)) 
HB.genes <- rownames(P2_BLBC@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
P2_BLBC[["percent.HB"]]<-PercentageFeatureSet(P2_BLBC, features=HB.genes) 

# P3_BLBC
HB_m <- match(HB.genes, rownames(P3_BLBC@assays$RNA)) 
HB.genes <- rownames(P3_BLBC@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
P3_BLBC[["percent.HB"]]<-PercentageFeatureSet(P3_BLBC, features=HB.genes) 



# ============================================================
# ---- 2.4 QC metric visualization
# ============================================================

# nFeature_RNA: number of genes detected per cell.
# nCount_RNA: total number of RNA molecules per cell.
# Generally, high values may indicate high cellular expression levels,
# but excessively high values can also be a sign of doublets or low-quality cells.

# Violin plots for all 5 QC metrics

CML1 <- merge(P1_CP, y=c(P1_BLBC, P2_CP, P2_BLBC, P3_BLBC))

slotNames(CML1)
p <- VlnPlot(CML1, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.HB"),log =T,pt.size =0.01)
ggsave(plot=p, file="CML_QC_before.pdf",width =12, height =12)

rm(CML1)


# ============================================================
#  ---- 2.5 Filter cells based on QC thresholds 
# ============================================================

# In practice, we usually only need to filter on percent.mt, i.e., to remove dead cells.
# How do we perform data filtering? We can do this by examining violin plots.
# In the plots, the distributions of nFeature_RNA, nCount_RNA, and percent.mt are crucial.
# We need to remove those discrete/outlier data points,
# because such abnormal data may affect the accuracy of downstream analyses.


# Thresholds should be adjusted according to the dataset.
#   nFeature_RNA > 200  - remove empty droplets / low-quality cells
#   nFeature_RNA < 6000 - remove doublets
#   percent.mt   < 10    - remove cells with high mitochondrial contamination
#   percent.ribo < 40   - optional; remove cells with extreme ribosomal content
#   percent.hb   < 20    - optional; remove erythrocyte-contaminated cells



dim(P1_CP)
# [1] 16228 10139
P1_CP <- subset(P1_CP, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000  & 
                   percent.mt < 10)
# Check number of cells after filtering
dim(P1_CP)
# [1] 16228  8707
rm(P1_CP1)

dim(P1_BLBC)
# [1] 16134  8938
P1_BLBC <- subset(P1_BLBC, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000  & 
                   percent.mt < 10)
dim(P1_BLBC)
# [1] 16134  7220
rm(P1_BLBC1)


dim(P2_CP)
# [1] 16285  9258
P2_CP <- subset(P2_CP, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000  & 
                   percent.mt < 20)
dim(P2_CP)
# [1] 16285  6270
rm(P2_CP1)

dim(P2_BLBC)
#[1] 16154 13230
P2_BLBC <- subset(P2_BLBC, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000  & 
                   percent.mt < 20)
dim(P2_BLBC)
#[1] 16154  8897
rm(P2_BLBC1)

dim(P3_BLBC)
#[1] 15188 11943
P3_BLBC <- subset(P3_BLBC, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000  & 
                   percent.mt < 10)
dim(P3_BLBC)
#[1] 15188 11013
rm(P3_BLBC1)


# ============================================================
# ---- 2.6 QC metric visualization after filtering
# ============================================================

# Violin plots for all 5 QC metrics

CML <- merge(P1_CP, y=c(P1_BLBC, P2_CP, P2_BLBC, P3_BLBC))

slotNames(CML)
p <- VlnPlot(CML, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.HB"),log =T,pt.size =0.01)
ggsave(plot=p, file="CML_QC_before.pdf",width =12, height =12)




# ============================================================
# Step 3: Normalization and Detect High Variable Genes analysis
# ============================================================

# ============================================================
# ---- 3.1 Normalize data
# ============================================================
# Normalize data to remove sequencing depth differences
# LogNormalize: divide counts by total counts per cell, multiply by scale.factor, then log1p
# Default method: LogNormalize
# scale.factor = 10000; Scale each cell to 10,000 counts

CML <- NormalizeData(object = CML, normalization.method = "LogNormalize", scale.factor = 10000)
#An object of class Seurat 
#18969 features across 42107 samples within 1 assay 
#Active assay: RNA (18969 features, 2000 variable features)
 #10 layers present: counts.P1_CP, counts.P1_BLBC, counts.P2_CP, counts.P2_BLBC, counts.P3_BLBC, data.P1_CP, data.P1_BLBC, data.P2_CP, data.P2_BLBC, data.P3_BLBC


# ============================================================
# ---- 3.2 Identify highly variable genes (HVGs)
# ============================================================
# These genes show high cell-to-cell variation and drive downstream clustering
# Default method: vst; Variance stabilizing transformation
# nfeatures = 2000; Select top 2,000 HVGs
CML  <- FindVariableFeatures(object=CML, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
 top10 <- head(VariableFeatures(CML), 10)
 top10
#   [1] "MPO"    "LGALS1" "AVP"    "ANXA1"  "LYZ"    "AZU1"   "JCHAIN" "UBE2C"  "IGKC"  
#[10] "MKI67



# ============================================================
# Step 4: Scaling and PCA
# ============================================================

# ---- 4.1 Scale data: mean = 0, variance = 1 for each gene (required for PCA)
####  we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
####  Shifts the expression of each gene, so that the mean expression across cells is 0
####  Scales the expression of each gene, so that the variance across cells is 1
####  This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
####  The results of this are stored in pbmc[["RNA"]]$scale.data
####  By default, only variable features are scaled.
####  You can specify the features argument to scale additional features


all.genes <- rownames(CML)
CML <- ScaleData(CML, features = all.genes)


# ---- 4.2 Principal Component Analysis (PCA) on HVGs
## Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to 
## choose a different subset (if you do want to use a custom subset of features, make sure you pass these to ScaleData first).
## For the first principal components, Seurat outputs a list of genes with the most positive and negative loadings, representing modules of genes that exhibit either correlation (or anti-correlation) 
## across single-cells in the dataset.

CML <- RunPCA(CML, features = VariableFeatures(object = CML), npcs = 50, verbose = FALSE)
slotNames(CML)
###  [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
 ###[6] "neighbors"    "reductions"   "images"       "project.name" "misc"        
###[11] "version"      "commands"     "tools"

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
# Examine and visualize PCA results a few different ways

# Inspect PCA results
print(CML[["pca"]], dims = 1:20, nfeatures = 10)

# Visualize PCA
VizDimLoadings(CML, dims = 1:2, reduction = "pca")
DimPlot(CML, reduction = "pca") + NoLegend()

p <- DimHeatmap(CML, dims = 1:15, cells = 500, balanced = TRUE)
ggsave(plot=p, file="CML.DimHeatmap.pdf")


# Determine number of PCs to retain (elbow plot)
# Typically choose the number at the elbow (e.g., 15-30)
p <- ElbowPlot(CML, ndims = 50)
ggsave(plot=p, file="CML.ElbowPlot.pdf")

## 


# ============================================================
# Step 5: Multi-sample integration using Harmony
# ============================================================
# Seurat v5 provides IntegrateLayers() as a unified interface for
# dataset integration. It corrects batch effects by finding shared
# biological structure across samples, so that cells of the same type
# co-cluster regardless of their sample origin.
# Choose one of the following methods:
##  •   Anchor-based CCA integration (method=CCAIntegration)
##  •   Anchor-based RPCA integration (method=RPCAIntegration)
##  •   Harmony (method=HarmonyIntegration)
##  •  JointPCAIntegration



# --- Harmony Integration ---
# Fast and widely used; requires the harmony package.
CML <- IntegrateLayers(
  object = CML, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


names(CML@reductions)
[1] "pca"     "harmony"

CML <- JoinLayers(CML)


table(CML@meta.data$orig.ident)
####   P1_CP  P1_BLBC   P2_CP      P2_BLBC      P3_BLBC 
###      8707            7220          6270         8897            11013 


str(CML@meta.data$orig.ident)
 ## chr [1:42107] "P1_CP" "P1_CP" "P1_CP" "P1_CP" "P1_CP" "P1_CP" ...



# ============================================================
# Step 6: Assign sample metadata
# ============================================================

##  ---- 6.1  Define different samples and add the sample information to the meta.data of the Seurat object.
########## 
Sample_group <- c(rep("P1_CP",8707),rep("P1_BLBC",7220),rep("P2_CP",6270),rep("P2_BLBC",8897),rep("P3_BLBC",11013))
str(Sample_group)
 ##  chr [1:42107] "P1_CP" "P1_CP" "P1_CP" "P1_CP" "P1_CP" "P1_CP" ...

sample_group_uniq <- unique(Sample_group)
Sample_group <- factor(Sample_group,levels=c("P1_CP","P2_CP","P1_BLBC","P2_BLBC","P3_BLBC"))
levels(Sample_group)
 ##  [1] "P1_CP" "P2_CP" "P1_BLBC"  "P2_BLBC"  "P3_BLBC"
CML@meta.data$Sample <- Sample_group
head(CML@meta.data)



##  ---- 6.2  Define different disease status and add the type information to the meta.data of the Seurat object.
########## 
Type_group <- c(rep("Chronic phase",8707),rep("Blastic crisis",7220),rep("Chronic phase",6270),rep("Blastic crisis",8897),rep("Blastic crisis",11013))
str(Type_group)
 ## chr [1:42107] "Chronic phase" "Chronic phase" "Chronic phase" "Chronic phase" ...

type_group_uniq <- unique(Type_group)
Type_group <- factor(Type_group,levels=c("Chronic phase","Blastic crisis"))
levels(Type_group)
## [1] "Chronic phase"  "Blastic crisis"
CML@meta.data$Type <- Type_group
head(CML@meta.data)

write.csv(CML@meta.data, file="CML_metaData.csv", row.names=TRUE)




# ============================================================
# Step 7: Cell Clustering and UMAP visualization
# ============================================================
# Construct KNN graph based on Euclidean distance in PCA space

#### ---7.1  find neighbors
###Use first 20 PCs (adjust based on ElbowPlot)
CML <- FindNeighbors(CML, reduction = "harmony", dims = 1:20)
##  Computing nearest neighbor graph
##  Computing SNN


#### ---7.2  find clusters
# Cluster cells using Louvain/Leiden algorithm
# Resolution: higher values yield more clusters
CML <- FindClusters(object = CML, resolution = 1.2, cluster.name = "harmony_clusters12")## 24 clusters
CML <- FindClusters(object = CML, resolution = 1.1, cluster.name = "harmony_clusters11")## 23 clusters
CML <- FindClusters(object = CML, resolution = 1, cluster.name = "harmony_clusters10")## 22 clusters
CML <- FindClusters(object = CML, resolution = 0.9, cluster.name = "harmony_clusters9")## 21 clusters
CML <- FindClusters(object = CML, resolution = 0.8, cluster.name = "harmony_clusters8")## 21 clusters
CML <- FindClusters(object = CML, resolution = 0.7, cluster.name = "harmony_clusters7")## 20 clusters
CML <- FindClusters(object = CML, resolution = 0.6, cluster.name = "harmony_clusters6")## 17 clusters
CML <- FindClusters(object = CML, resolution = 0.5, cluster.name = "harmony_clusters5")## 14 clusters
CML <- FindClusters(object = CML, resolution = 0.4, cluster.name = "harmony_clusters4")## 14 clusters
CML <- FindClusters(object = CML, resolution = 0.3, cluster.name = "harmony_clusters3")## 10 clusters
CML <- FindClusters(object = CML, resolution = 0.2, cluster.name = "harmony_clusters2")## 8 clusters
CML <- FindClusters(object = CML, resolution = 0.1, cluster.name = "harmony_clusters1")## 5 clusters


#### ---7.3  UMAP for visualization
CML <- RunUMAP(CML, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")


# UMAP plot
pdf("UMAP_plot_integrated.harmony_by_Sample.pdf",height = 5, width = 6)
p1 <- DimPlot(object = CML,  pt.size = 0.5, reduction = 'umap.harmony', cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Sample))),group.by = "Sample")
print(p1)
dev.off()

pdf("UMAP_plot_integrated.harmony_by_Samplesplit.pdf",height = 6, width = 18)
p1 <- DimPlot(object = CML,  split.by = "Sample", pt.size = 0.5, reduction = 'umap.harmony', cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Sample))),group.by = "Sample")
print(p1)
dev.off()

pdf("UMAP_plot_integrated.harmony_by_Type.pdf",height = 5, width = 6)
p1 <- DimPlot(object = CML,  pt.size = 0.5, reduction = 'umap.harmony', cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Type))),group.by = "Type")
print(p1)
dev.off()

pdf("UMAP_plot_integrated.harmony_by_splitType.pdf",height = 4, width = 10)
p1 <- DimPlot(object = CML,  split.by = "Type",pt.size = 0.5, reduction = 'umap.harmony', cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Type))),group.by = "Type")
print(p1)
dev.off()

pdf("UMAP_plot_integrated.harmony_by_Typesplit.pdf",height = 6, width = 18)
p1 <- DimPlot(object = CML,  split.by = "Sample", pt.size = 0.5, reduction = 'umap.harmony', cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Type))),group.by = "Type")
print(p1)
dev.off()



pdf("all_umap.harmony.pdf")
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters1))), label = TRUE, group.by = "harmony_clusters1")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters2))), label = TRUE, group.by = "harmony_clusters2")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters3))), label = TRUE, group.by = "harmony_clusters3")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters4))), label = TRUE, group.by = "harmony_clusters4")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters5))), label = TRUE, group.by = "harmony_clusters5")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters6))), label = TRUE, group.by = "harmony_clusters6")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters7))), label = TRUE, group.by = "harmony_clusters7")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters8))), label = TRUE, group.by = "harmony_clusters8")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters9))), label = TRUE, group.by = "harmony_clusters9")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters10))), label = TRUE, group.by = "harmony_clusters10")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters11))), label = TRUE, group.by = "harmony_clusters11")+theme(aspect.ratio=1)
DimPlot(CML, reduction = "umap.harmony", pt.size = 0.5, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$harmony_clusters12))), label = TRUE, group.by = "harmony_clusters12")+theme(aspect.ratio=1)
dev.off()



## ---7.4  Selected the resolution to detemine the final clusters

## Draw a cluster tree to help us detemine the resolution 
library(clustree)
p <- clustree(CML,prefix = "harmony_clusters")
ggsave(plot=p, file="CML.clustree_harmony.pdf",height = 28, width = 14)


##### Set "harmony_clusters4" as the default cluster
Idents(CML) <-  "harmony_clusters9"

#At this point, we can observe that the classification numbers of the clusters start from 0, and the cell categories (of factor data type)
table(CML@active.ident)
##      0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
## 5858 4118 3838 3180 3057 2911 2805 1990 1846 1678 1568 1533 1488 1464 1244  900  846 
##   17   18   19   20 
##  698  490  428  167 


##### Next, we add a "Cluster" attribute to meta.data, so that the classification numbers of the clusters start from 1.
cluster_add1_idents <- as.integer(CML@active.ident)
cluster_add1_idents <- factor(cluster_add1_idents,levels=1:length(unique(cluster_add1_idents)))
names(cluster_add1_idents) <- rownames(CML@meta.data)

CML@active.ident <- cluster_add1_idents
CML@meta.data$cluster <- cluster_add1_idents
head(CML@meta.data)
table(CML@meta.data$cluster)
##     1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
## 5858 4118 3838 3180 3057 2911 2805 1990 1846 1678 1568 1533 1488 1464 1244  900  846 
 ##  18   19   20   21 
 ## 698  490  428  167 


pdf("UMAP_plot_by_harmonyCluster.pdf",height = 8, width = 10)
p1 <- DimPlot(object = CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))), pt.size = 0.5, label = TRUE, reduction = "umap.harmony", 
    group.by = "cluster", label.size = 6)
print(p1)
dev.off()

pdf("UMAP_plot_by_harmonyCluster_splitSample.pdf",width = 24)
p1 <- DimPlot(CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))), pt.size = 0.5, label = TRUE, reduction = "umap.harmony", 
    group.by = "cluster", label.size = 4, split.by = "Sample")
print(p1)
dev.off()

pdf("UMAP_plot_by_harmonyCluster_splitType.pdf",width = 16)
p1 <- DimPlot(CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))), pt.size = 0.5, label = TRUE, reduction = "umap.harmony", 
    group.by = "cluster", label.size = 4, split.by = "Type")
print(p1)
dev.off()






# ============================================================
# Step 8: Cell type annotation
# ============================================================
# ----- 8.1 Find marker genes
# ------------------------------------------------
# Find marker genes for each cluster
##  min.pct: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. 
#### Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
####                   p_val  avg_log2FC pct.1 pct.2     p_val_adj cluster    gene


DefaultAssay(object = CML) <- "RNA"
var_genes <- VariableFeatures(CML)
library(presto)

#### wilcox method
diff_wilcox <- FindAllMarkers(object = CML,  features = var_genes, only.pos = TRUE, logfc.threshold = 0.25, test.use="wilcox", min.pct = 0.1)


####  Visualizing the expression of marker genes in each cluster
for(a_cluster in levels(CML@active.ident)){
   mks <- diff_wilcox[diff_wilcox$cluster==a_cluster,]
   mks1 <- mks[order(mks$avg_log2FC,decreasing = T),"gene"][1:12]
   otfile <- paste("harmonyCluster",a_cluster,"_wilcox_markers_umap.pdf", sep="")
   pdf(file = otfile, width = 16, height = 10)
   p <- FeaturePlot(CML, features = mks1, cols = c("grey", "red"), reduction = "umap.harmony")
   print(p)
   dev.off()
   outfile <- paste("harmonyCluster",a_cluster,"_VlnPlot_wilcox_markers.pdf", sep="")
   pdf(file = outfile, width = 8, height = 20)
   p1 <- VlnPlot(CML, features = mks1, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))), pt.size = 0, ncol = 1)
   print(p1)
   dev.off()
}


####  Generating the table of the differential expressed marker genes for each cluster
for(a_cluster in levels(CML@active.ident))
{
     subdiff <- diff_wilcox[diff_wilcox$cluster == a_cluster,]
     otfile <- paste("harmonyCluster",a_cluster,"_diff_wilcox.xls", sep="")
     write.table(subdiff, sep="\t", file=otfile, row.names=FALSE, quote=FALSE)
  #    n <- dim(subdiff_info)[2]
     otfile <- paste("harmonyCluster",a_cluster, "_diff_wilcox_significant.xls", sep="")
     write.table(subdiff[subdiff[,1] <= 0.05 & abs(subdiff[,2]) >= 0.585,], sep="\t", file=otfile, row.names=FALSE, quote=FALSE)
}



# ----- 8.2 Visualizing the expression of well konwon lineage-specific marker genes
# ------------------------------------------------
Idents(CML) <- "cluster"
Cycling_marker=c("CD34", "MCM7","TYMS", "PCNA","GINS2", "MKI67", "TOP2A", "ASPM", "PCLAF", "CCNB1")

p <- VlnPlot(CML, group.by = "cluster", features =Cycling_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_Cycling_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =Cycling_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_Cycling_marker_Data.pdf",width = 16, height = 10)

HSPC_marker=c("CD34","FLT3","SOX4","AVP", "SPINK2","LAPTM4B","CAVIN2","NPR3","PRSS57","AIF1","CRHBP", "HTR1F", "MYCT1", "MSRB3",  "ANGPT1", "FAM30A", "GATA2")

p <- VlnPlot(CML, group.by = "cluster", features =HSPC_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_HSPC_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =HSPC_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_HSPC_marker_Data.pdf",width = 16, height = 10)

ErythrocyticCell_marker=c("HBA1", "HBA2","HBB", "HBG2", "AHSP", "CA1", "GYPA", "HBG1", "SLC4A1","APOC1","KLF1", "HBD", "FAM178B", "RHCE", "BLVRB")

p <- VlnPlot(CML, group.by = "cluster", features =ErythrocyticCell_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_ErythrocyticCell_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =ErythrocyticCell_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_ErythrocyticCell_marker_Data.pdf",width = 16, height = 10)

progenitor_marker=c("CD34", "PRSS57","DNTT", "SMIM3","MT1X", "SOCS2", "HPS4", "MLXIP", "VPREB1", "TSC22D","CEBPA","APOE")

p <- VlnPlot(CML, group.by = "cluster", features =progenitor_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_progenitor_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =progenitor_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_progenitor_marker_Data.pdf",width = 16, height = 10)

GMP_marker=c("MPO", "ELANE","CTSG", "AZU1","LYST", "LYZ", "CEBPD", "MNDA")

p <- VlnPlot(CML, group.by = "cluster", features =GMP_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_GMP_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =progenitor_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_GMP_marker_Data.pdf",width = 16, height = 10)

MEP_marker=c("GATA1", "KLF1","APOE", "AHSP", "CA1")

p <- VlnPlot(CML, group.by = "cluster", features =MEP_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_MEP_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =progenitor_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_MEP_marker_Data.pdf",width = 16, height = 10)


MyloidCell_marker=c("CD68","FCGR3B", "CD163","ITGAM", "CD1C", "MPO",  "ITGB2",  "CD244", "TNFRSF8", "CD33", "ANPEP")

p <- VlnPlot(CML, group.by = "cluster", features =MyloidCell_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_MyloidCell_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =MyloidCell_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_MyloidCell_marker_Data.pdf",width = 16, height = 10)

Monocyte_marker=c("S100A8", "MNDA", "CD68","FCN1","CD14", "S100A12", "FCGR3A", "CSF1R", "RNASE2", "UGCG")

p <- VlnPlot(CML, group.by = "cluster", features =Monocyte_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_Monocyte_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =Monocyte_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_Monocyte_marker_Data.pdf",width = 16, height = 10)

ClaMOcell_marker=c("S100A8", "S100A12","CD14", "LYZ", "CD68", "RETN", "CFD", "SERPINA1", "LGALS2", "CEBPB", "VCAN")

p <- VlnPlot(CML, group.by = "cluster", features =ClaMOcell_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_ClaMOcell_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =ClaMOcell_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_ClaMOcell_marker_Data.pdf",width = 16, height = 10)

NonClaMOcell_marker=c("FCGR3A", "AIF1","MS4A7", "LST1", "FTH1", "C1QA","MAFB", "CDKN1C", "HES4", "HMOX1", "CSF1R", "LYPD2")

p <- VlnPlot(CML, group.by = "cluster", features =NonClaMOcell_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_NonClaMOcell_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =NonClaMOcell_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_NonClaMOcell_marker_Data.pdf",width = 16, height = 10)

pDC_marker=c("CD1C","FCER1A","IL3RA","UGCG","LILRA4","CLEC4C", "SMPD3", "SERPINF1","MS4A7", "RNASE2","PPBP", "PF4","TUBB1", "GP9")

p <- VlnPlot(CML, group.by = "cluster", features =pDC_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_pDC_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =pDC_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_pDC_marker_Data.pdf",width = 16, height = 10)

MacroCHECK_marker=c("CALR","CD47", "SIRPA", "CD24","SIGLEC10","VSIG4","VCAM1","SLAMF7")

p <- VlnPlot(CML, group.by = "cluster", features =MacroCHECK_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_MacroCHECK_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =MacroCHECK_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_MacroCHECK_marker_Data.pdf",width = 16, height = 10)


Tcell_marker=c("CD3D", "CD3E","CD3G", "CD8A", "CD4", "CD2", "CD5", "CD7", "CD8B", "CD27", "CD70")

p <- VlnPlot(CML, group.by = "cluster", features =Tcell_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_Tcell_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =Tcell_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_Tcell_marker_Data.pdf",width = 16, height = 10)

naiveT_marker=c("CD3D", "CD3E","CD3G",  "CD4","CD27",  "LEF1", "IL7R", "IL32")

p <- VlnPlot(CML, group.by = "cluster", features =naiveT_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_naiveT_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =naiveT_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_naiveT_marker_Data.pdf",width = 16, height = 10)

effectT_marker=c("CD3D", "CD3E","CD3G",  "CD8A","CD8B", "CD70", "CD27",  "CCL4", "CCL5", "GZMK", "GZMA", "CST7", "NKG7")

p <- VlnPlot(CML, group.by = "cluster", features =effectT_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_effectT_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =effectT_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_effectT_marker_Data.pdf",width = 16, height = 10)

NKcell_marker=c("GNLY", "NKG7" ,"CTSW", "FGFBP2","KLRD1", "KLRB1", "PRF1","FCER1G", "KLRC1", "KLRF1", "GZMB", "TRDC", "CD247")

p <- VlnPlot(CML, group.by = "cluster", features =NKcell_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_NKcell_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =NKcell_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_NKcell_marker_Data.pdf",width = 16, height = 10)

NKT_marker=c("CD3D", "TRBC1", "CD8A","CCR7","GPR183", "GZMA", "GZMB","GZMH","KLRB1", "GNLY", "SPON2", "CCL3","CCL5","HSPA6")

p <- VlnPlot(CML, group.by = "cluster", features =NKT_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_NKT_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =NKT_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_NKT_marker_Data.pdf",width = 16, height = 10)


NKTCHECK_marker=c( "PRF1","GZMA", "GZMB","GZMH","CST7", "KLRB1", "KLRD1","KLRG1", "KLRC1","KLRC2", "PDCD1", "HAVCR2","CTLA4", "LAG3","TIGIT")

p <- VlnPlot(CML, group.by = "cluster", features =NKTCHECK_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_NKTCHECK_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =NKTCHECK_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_NKTCHECK_marker_Data.pdf",width = 16, height = 10)

TCHECK_marker=c("CD274", "PDCD1", "CTLA4","CD80","CD86","CD28", "LAG3", "HAVCR2","TIGIT","TNFRSF9","TOX","TOX2","CMTM4","SATB1","VSIG4","TGIF1","VCAM1")

p <- VlnPlot(CML, group.by = "cluster", features =TCHECK_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_TCHECK_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =TCHECK_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_TCHECK_marker_Data.pdf",width = 16, height = 10)

NKCHECK_marker=c("KLRB1", "CLEC2D", "KLRG1","KLRC2","LILRB1","CD96","PVR","IL18R1","IL18RAP","KLRD1")

p <- VlnPlot(CML, group.by = "cluster", features =NKCHECK_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_NKCHECK_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =NKCHECK_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_NKCHECK_marker_Data.pdf",width = 16, height = 10)

THCHECK_marker=c("NCR3", "IL7R", "CD247","CCR7","KIR2DL3","IL17F", "IL23R","RORC","CD96","FOXP3")

p <- VlnPlot(CML, group.by = "cluster", features =THCHECK_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_THCHECK_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =THCHECK_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_THCHECK_marker_Data.pdf",width = 16, height = 10)


B_marker=c("CD24", "EBF1","MME","VPREB1","PAX5","CD19", "CD79A","MS4A1","BANK1","CD27","CD70","MZB1","IGLL5","JCHAIN")

p <- VlnPlot(CML, group.by = "cluster", features =B_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_B_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =B_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_B_marker_Data.pdf",width = 16, height = 10)

naiveB_marker=c("IGKC", "IGHD","IGHM","IGLL1","CD83", "VPREB3", "LHPP", "NSMCE", "ADAM28")

p <- VlnPlot(CML, group.by = "cluster", features =naiveB_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_v_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =naiveB_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_naiveB_marker_Data.pdf",width = 16, height = 10)

MemoryBCell_marker=c("CD79A","MS4A1","IFI44L","CD27","CD70","GPR183", "EBI3","IGLL1","IGHG1","CD81")

p <- VlnPlot(CML, group.by = "cluster", features =MemoryBCell_marker,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_MemoryBCell_marker_Data.pdf",width =8, height =24)
p <- FeaturePlot(CML, features =MemoryBCell_marker, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_MemoryBCell_marker_Data.pdf",width = 16, height = 10)


pro_BT_marker=c("CD34", "DNTT","SOCS2", "HPS4","VPREB1", "CD24","MME","CD79A","MS4A1","CD3D","MPO","ANPEP")

p <- VlnPlot(CML, group.by = "cluster", split.by = "Type", features =pro_BT_marker,  pt.size = 0, ncol = 1,  layer = "data", 
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Type))))
ggsave(plot=p, file="vlnCML_resolution0.9_pro_BT_marker_Data.pdf",width =8, height =24)



AllCell_marker=c("CD34", "FLT3","SOX4", "GATA1", "KLF1",  "HBB", "PF4", "MPO","LYZ", "VCAN","S100A8", "CD68", "CD163","MAFB","CD1C","LILRA4","CD79A","MS4A1", "PAX5","CD27","CD70", 
   "MZB1", "LEF1",  "CD3E", "CD8A", "CCL5","KLRB1", "CD247", "PTPRC")

p <- DotPlot(CML, features =AllCell_marker, cols=c("blue", "red"))
ggsave(plot=p, file="dotplotCML_resolution0.9_AllCell_marker_Data.pdf",width =20, height =12)
p <- VlnPlot(CML, group.by = "cluster", features =AllCell_marker,  pt.size = 0, ncol = 2,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_AllCell_marker_Data.pdf",width =10, height =24)
p <- FeaturePlot(CML, features =AllCell_marker, cols = c("grey", "red"), order = T, reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_AllCell_marker_Data.pdf",width = 16, height = 10)

CML <- readRDS("/public/zxh/projects/CMLblast/figure/CMLharmony/CML_phase.rds")
Idents(CML) <- "cluster"

AllCell_marker1=c("CRHBP","AVP", "CD34", "FLT3","KIT", "GATA2", "GATA1","ITGA2B",  "KLF1", "TFRC", "CD36", "HBD", "HBA1", "CA1","GYPA","ALAS2","VWF","GP9","PF4",
    "CSF3R", "MPO", "ELANE", "LYZ","IRF8", "CD14", "MME", "VPREB1", "CD79A","IGHM", "IGHD", "EBF1", "VPREB3", "CD24","PAX5", "E2F2", "CD19","MS4A1", "MKI67",  "PLK1")

p <- DotPlot(CML, features =AllCell_marker1, cols=c("grey", "red"))
ggsave(plot=p, file="dotplotCML_resolution0.9_AllCell_marker1_Data.pdf",width =24, height =10)
p <- VlnPlot(CML, group.by = "cluster", features =AllCell_marker1,  pt.size = 0, ncol = 2,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_AllCell_marker1_Data.pdf",width =10, height =36)

AllCell_marker2=c("AVP", "CD34", "FLT3","GATA2", "GATA1","KLF1","CD36", "HBD","CA1","GP9", "MPO", "LYZ", "MME", "CD79A","CD24","EBF1","PAX5","CD19","MS4A1", "MKI67","CD3D")

p <- DotPlot(CML, features =AllCell_marker2, cols=c("grey", "red"))
ggsave(plot=p, file="dotplotCML_resolution0.9_AllCell_marker2_Data.pdf",width =14, height =6)
p <- VlnPlot(CML, group.by = "cluster", features =AllCell_marker2,  pt.size = 0, ncol = 2,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_AllCell_marker2_Data.pdf",width =10, height =20)





# ----- 8.3 Inferring the relative differentiation states to determine the Hematopoietic Stem/Progenitor Cells (HSPC)
# ------------------------------------------------

# Load CytoTRACE2 package
library(CytoTRACE2)

# Support direct reading of Seurat object files
cytotrace2_result<-cytotrace2(CML, is_seurat=TRUE, slot_type="counts", species='human')

# Check the return value
names(cytotrace2_result)
# [1] "RNA"          "RNA_nn"       "RNA_snn"      "pca"          "harmony"     
# [6] "tsne.harmony" "umap.harmony"

head(cytotrace2_result,2)
 #                              barcode   EBNA.2    LMP.1 CytoTRACE2_Score
# AAACCCACAGTCTGGC-1 AAACCCACAGTCTGGC-1 Negative Negative        0.0317993
# AAACCCAGTGCTATTG-1 AAACCCAGTGCTATTG-1 Negative Negative        0.1070740
 #                   CytoTRACE2_Potency CytoTRACE2_Relative preKNN_CytoTRACE2_Score
# AAACCCACAGTCTGGC-1     Differentiated          0.03795011              0.01760164
# AAACCCAGTGCTATTG-1     Differentiated          0.13004431              0.10211434
#                    preKNN_CytoTRACE2_Potency
# AAACCCACAGTCTGGC-1            Differentiated
# AAACCCAGTGCTATTG-1            Differentiated


# Note: CytoTRACE2_Relative is the specific numerical result of the score;
#       CytoTRACE2_Potency is the six-category differentiation potency classification assigned by CytoTRACE2.
#
# CytoTRACE2_Score represents the differentiation potential score for each cell
# (range 0 to 1) assigned by the software; higher scores indicate greater
# differentiation potential and a lower degree of differentiation.
#
# CytoTRACE2_Potency is a differentiation potency classification based on the
# differentiation potential score:
#   differentiated - differentiated cells;
#   multipotent    - can differentiate into multiple types of cells;
#   oligopotent    - can differentiate into a few types of cells;
#   unipotent      - can differentiate into one type of cell.
#
# CytoTRACE2_Relative is the standardized (Z-Score) differentiation potential score,
# ranging from [0,1], which is convenient for setting gradient colors in plots.



# Result visualization; CytoTRACE2 provides multiple visualization outputs.
annotation <- data.frame(phenotype=cytotrace2_result@meta.data$cluster) %>% set_rownames(., colnames(cytotrace2_result))
plots <- plotData(cytotrace2_result=cytotrace2_result, annotation=annotation, is_seurat=TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
p <- (p1+p2+p3+p4) + plot_layout(ncol = 2)
ggsave(plot=p, "CytoTRACE2_UMAP_plotData.pdf")


p <- DimPlot(cytotrace2_result, reduction = "umap.harmony", split.by = "Sample", cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(cytotrace2_result@meta.data$CytoTRACE2_Potency))), label = FALSE, group.by = "CytoTRACE2_Potency")+theme(aspect.ratio=1)
ggsave(plot=p, "CytoTRACE2_UMAP_Potency_splitSample.pdf",width =20, height =5)

p <- FeaturePlot(cytotrace2_result, features = "CytoTRACE2_Relative", cols = viridis::viridis(100),  order = TRUE, reduction = "umap.harmony")
ggsave(plot=p, "CytoTRACE2_UMAP_Relative.pdf",width =6, height =5)

p <- VlnPlot(cytotrace2_result, group.by = "cluster", features ="CytoTRACE2_Score",  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(cytotrace2_result@meta.data$cluster))))
ggsave(plot=p, file="vln_CytoTRACE2_score.pdf",width =6, height =2)

p <- ggplot(cytotrace2_result@meta.data, aes(x=cluster))
p <- p + geom_bar(aes(fill=CytoTRACE2_Potency), position="fill")
ggsave("cluster.harmony_PotencydisperAllCell.pdf",device = "pdf")

p <- ggplot(cytotrace2_result@meta.data, aes(x=Sample))
p <- p + geom_bar(aes(fill=CytoTRACE2_Potency), position="fill")
ggsave("sample.harmony_PotencydisperAllCell.pdf",device = "pdf")

############




# ----- 8.4 Cell-type inference using scHCL
# ------------------------------------------------

#### The return of scHCL() is a list which contains 4 parts.
#### cors_matrix: Pearson correlation coefficient matrix of each cell and cell type.
#### top_cors: equals to numbers_plot
#### scHCL: the most relevant cell type for each query cell
#### scHCL_probility: the top n relevant cell types for each query cell

# scHCL has two parameters , single cell expression matrix(scdata) and 
# the number of most similar cell types
library(scHCL)

########### cell type annotation by clusters
avg_exp_cluster <- sapply(1:21,function(x){
  data <- CML@assays$RNA@layers$counts[,CML@meta.data$cluster == x]
  data1 <- apply(data, 2, function(i){(i/sum(i))*100000})
  data2 <- apply(data1, 1, mean)
  return(data2)
})

CML.AVG <- AverageExpression(object = CML, assays= "RNA", layer = "counts", group.by = "cluster", return.seurat = TRUE)
cluster_mean_counts <- GetAssayData(CML.AVG[["RNA"]], layer = "counts")
cluster_mean_counts_matrix <- as.matrix(cluster_mean_counts)

ct_cluster <- scHCL(scdata = cluster_mean_counts_matrix, numbers_plot = 10)
head(ct_cluster$scHCL_probility)
ct_cluster$scHCL_probility$Cell <- as.character(ct_cluster$scHCL_probility$Cell )
ct_cluster$scHCL_probility$`Cell type` <- as.character(ct_cluster$scHCL_probility$`Cell type`)
md_cluster <- CML.AVG@meta.data
md_cluster$Cell <- rownames(CML.AVG@meta.data)
CML.AVG@meta.data$`Cell type` <- md_cluster$`Cell type`
CML.AVG@meta.data$Score <- md_cluster$Score
write.csv(ct_cluster$scHCL_probility, file = "CML.AVG_scHCL_ByCluster.csv")






# ----- 8.5 Mapping cluster ids to Cell types by manual
# ------------------------------------------------
identify_celltypes <- c("proB","preB","preB","MEP","preB","EEP","ImmatureB","EEP","preB","proB","preB","GMP","BMEP","HSPC", "preB", "EEP","preB","ProMono","proB","EEP","T")
names(identify_celltypes) <- 1:length(identify_celltypes)
clusters <- CML@meta.data$cluster
celltypes <- identify_celltypes[clusters]
CML@meta.data$Celltype <- celltypes
head(CML@meta.data$Celltype)

table(CML@meta.data$Celltype)
##          ImmatureB    BMEP     EEP     GMP    HSPC     MEP    preB    proB    ProMono       T 
##                   2805    1488    6229      1533    1464       3180   16517    8026     698           167 



## Set the level of cell types
CML@meta.data$Celltype <- factor(CML@meta.data$Celltype, level=c("HSPC","BMEP","MEP","EEP","GMP","ProMono","proB","preB","ImmatureB","T"))


## UMAP visulization based on cell types
pdf("UMAP_plot_resolution0.9_by_Celltype.pdf", width=7, height=6)
p1 <- DimPlot(object = CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Celltype))), pt.size = 0.5, reduction = "umap.harmony", group.by = "Celltype")
print(p1)
dev.off()

pdf("UMAP_plot_resolution0.9_by_Celltype_splitType.pdf", width=10, height=6)
p1 <- DimPlot(object = CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Celltype))), pt.size = 0.5, reduction = "umap.harmony", group.by = "Celltype", split.by = "Type")
print(p1)
dev.off()

pdf("UMAP_plot_resolution0.9_by_Celltype_splitSample.pdf", width=28, height=5)
p1 <- DimPlot(object = CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Celltype))), pt.size = 0.5, reduction = "umap.harmony",  group.by = "Celltype", split.by = "Sample")
print(p1)
dev.off()


## Draw dotplots and violin plots of representative marker genes for each cell type
celltypes <- CML@meta.data$Celltype
names(celltypes) <- rownames(CML@meta.data)
CML@active.ident <- celltypes

AllCell_marker3=c("AVP", "CD34", "FLT3","GATA2","KLF1","CA1","GP9", "MPO","LYZ", "MME", "CD79A","CD24","EBF1","CD19","MS4A1", "CD3D","MKI67","PTPRC")

p <- DotPlot(CML, features =AllCell_marker3, cols=c("grey", "#498EA4"))
ggsave(plot=p, file="dotplotCML_resolution0.9_AllCell_marker3_Data_byCelltype.pdf",width =8, height =3)
p <- VlnPlot(CML, group.by = "Celltype", features =AllCell_marker3,  pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_AllCell_marker3_byCelltype.pdf",width = 4, height =34)

AllCell_marker4=c("AVP", "KLF1","CA1","GP9", "MPO","LYZ", "MME", "CD24","MS4A1", "CD3D","MKI67","PTPRC")
p <- FeaturePlot(CML, features =AllCell_marker4, cols = c("grey", "red"), reduction = "umap.harmony")
ggsave(plot=p, file="FeaturePlot_resolution0.9_AllCell_marker4_Data.pdf",width = 16, height = 10)


AllCell_marker5=c("CD34","GATA1","KLF1", "CD79A","CD24","DNTT","CD19","VPREB1", "VPREB3","RAG1","RAG2","CD83")
ggsave(plot=p, file="dotplotCML_resolution0.9_AllCell_marker3_Data_byCelltype.pdf",width =8, height =3)
p <- VlnPlot(CML, group.by = "Celltype", features =AllCell_marker5, split.by = "Type", pt.size = 0, ncol = 1,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$cluster))))
ggsave(plot=p, file="vlnCML_resolution0.9_AllCell_marker5_byCelltype.pdf",width = 6, height =24)

upCluster_marker=c("ITGA6","C1QTNF4","ECM1","MDM2","CDKN2A","MSH6","VPREB3","JCHAIN","LTB","CD37")
p <- VlnPlot(sceHSPC, group.by = "cluster",features =upCluster_marker,  split.by = "Type", pt.size = 0, ncol = 5,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceHSPC@meta.data$Type))))
ggsave(plot=p, file="vlnsceHSPC_upCluster_Data.pdf",width =10, height =3.5)




# ----- 8.6 Draw dotplots and violin plots of representative marker genes for myeloid clusters and lymphoid clusters respectively
# ------------------------------------------------


##  For all 9 myeloid clusters  ###
Idents(CML) <- "cluster"
sceM <-  subset(CML, idents = c(4, 6, 8, 12,13, 14, 16, 18, 20))

sceM@meta.data$cluster <-  factor(sceM@meta.data$cluster, levels = c("14","12","18","13","4","16", "6","20","8"))
heteroM_marker=c("CD34","AVP","CRHBP","SPINK2","RAB32","AZU1","ELANE","SELL","CEBPD","S100A9","FCER1G","TYROBP","GATA1","GATA2","PBX1","RAB27B",
    "FCER1A","ITGA2B","APOE","PLK1","HBA1","CA2")
p <- VlnPlot(sceM, group.by = "cluster",features =heteroM_marker,  pt.size = 0, ncol = 2,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceM@meta.data$cluster))))
ggsave(plot=p, file="vlnSCEm_heteroM_marker_Data.pdf",width =6, height =16)
p <- DotPlot(sceM, features =heteroM_marker, group.by = "cluster", cols=c("grey", "#498EA4"))
ggsave(plot=p, file="dotplotSCEm_heteroM_marker_Data.pdf",width =8, height =3)

sceM@meta.data$Celltype <- factor(CML@meta.data$Celltype, level=c("HSPC","BMEP","MEP","EEP","GMP","ProMono","proB","preB","B","T"))

p <- ggplot(sceM@meta.data, aes(x=Celltype))
p <- p + geom_bar(aes(fill= Sample),  position="fill")
ggsave("celltype_sample_resolution0.9_dodge.pdf",device = "pdf")


##  For all 11 lymphoid clusters  ###
sceB <-  subset(CML, idents = c(1, 2, 3, 5, 7, 9, 10, 11, 15, 17, 19))

sceB@meta.data$cluster <-  factor(sceB@meta.data$cluster, levels = c("1","10","19","11","2","3", "5","15","17","9","7"))
heteroB_marker=c("CD34","DNTT","CD24","ECM1","ITGA6","CTHRC1","JCHAIN","ADGRF1","GRASP","TK1","S100A16","ID2",
    "VPREB1","TCL1A","MME","VPREB3","EBF1","PAX5","IGLL1","IGHD","IFI27","CD70")
p <- VlnPlot(sceB, group.by = "cluster",features =heteroB_marker,  pt.size = 0, ncol = 2,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceB@meta.data$cluster))))
ggsave(plot=p, file="vlnSCEb_heteroB_marker_Data.pdf",width =7, height =17)
p <- DotPlot(sceB, features =heteroB_marker, group.by = "cluster", cols=c("grey", "#498EA4"))
ggsave(plot=p, file="dotplotSCEb_heteroB_marker_Data.pdf",width =8, height =4)


p <- ggplot(sceB@meta.data, aes(x=Celltype))
p <- p + geom_bar(aes(fill= Sample),  position="fill")
ggsave("celltype_sampleB_resolution0.9_dodge.pdf",device = "pdf")


# Draw violin plots of representative marker genes for HSPC
Idents(CML) <- "cluster"
sceHSPC <-  subset(CML, idents = c(14, 11))

upHSPC_marker=c("CD24","S100A10","LMNA","GBP2","ANXA2")
p <- VlnPlot(sceHSPC, group.by = "cluster",features =upHSPC_marker,  split.by = "Type", pt.size = 0, ncol = 5,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceHSPC@meta.data$Type))))
ggsave(plot=p, file="vlnsceHSPC_UPgenes_Data.pdf",width =10, height =2)

upCluster_marker=c("ITGA6","C1QTNF4","ECM1","MDM2","CDKN2A","MSH6","VPREB3","JCHAIN","LTB","CD37")
p <- VlnPlot(sceHSPC, group.by = "cluster",features =upCluster_marker,  split.by = "Type", pt.size = 0, ncol = 5,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceHSPC@meta.data$Type))))
ggsave(plot=p, file="vlnsceHSPC_upCluster_Data.pdf",width =10, height =3.5)





# ============================================================
# Step: 9 Differential expression analysis
# ============================================================
# Compare two groups 

"Chronic phase","Blastic crisis"
Idents(CML) <- "Type"
CMLmarkers <- FindMarkers(CML, ident.1 = "Blastic crisis", ident.2 = "Chronic phase", verbose = FALSE, test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)
write.csv(CMLmarkers, file="CMLmarkers_byType_wilcox.csv", row.names=TRUE)

"P1_CP","P2_CP","P1_BLBC","P2_BLBC","P3_BLBC"
Idents(CML) <- "Sample"
case1BLmarkers <- FindMarkers(CML, ident.1 = "P1_BLBC", ident.2 = "P1_CP", verbose = FALSE, test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)
write.csv(case1BLmarkers, file="case1BLBCmarkers_wilcox.csv", row.names=TRUE)
case2BLmarkers <- FindMarkers(CML, ident.1 = "P2_BLBC", ident.2 = "P2_CP", verbose = FALSE, test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)
write.csv(case2BLmarkers, file="case2BLBCmarkers_wilcox.csv", row.names=TRUE)


CML@meta.data$Type.Celltype <- paste(CML@meta.data$Type, CML@meta.data$Celltype, sep = "_")
Idents(CML) <- "Type.Celltype"
HSPC_BCmarkers <- FindMarkers(CML, ident.1 = "Blastic crisis_HSPC", ident.2 = "Chronic phase_HSPC", verbose = FALSE, test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)
write.csv(HSPC_BCmarkers, file="HSPC_BCmarkers_wilcox.csv", row.names=TRUE)

CML@meta.data$Type.cluster <- paste(CML@meta.data$Type, CML@meta.data$cluster, sep = "_")
Idents(CML) <- "Type.cluster"
cluster11_BCmarkers <- FindMarkers(CML, ident.1 = "Blastic crisis_11", ident.2 = "Chronic phase_11", verbose = FALSE, test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)
write.csv(cluster11_BCmarkers, file="cluster11_BCmarkers_wilcox.csv", row.names=TRUE)




# ============================================================
# Step: 10 Gene signature scoring
# ============================================================

# ----- 10.1 CellCycleScoring
# ------------------------------------------------

 ## Obtain the marker genes in G2M phase 
g2m_genes <- cc.genes$g2m.genes
# Extract the G2M phase marker genes from the matrix
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(CML)) 
#Obtain the marker genes in S phase
s_genes <- cc.genes$s.genes 
# Extract the S phase marker genes from the matrix    
s_genes <- CaseMatch(search=s_genes, match=rownames(CML)) 
# Cell cycle scoring using the CellCycleScoring function
CML <- CellCycleScoring(CML, g2m.features=g2m_genes, s.features=s_genes)


# Visualizing the distribution of markers throughout the entire cell cycle
phase <- factor(CML@meta.data$Phase)
names(phase) <- rownames(CML@meta.data)
CML@active.ident <- phase

p <- ggplot(CML@meta.data, aes(x=cluster))
p <- p + geom_bar(aes(fill= factor(Phase)),  position="fill")
ggsave("cluster.harmony_Phase_resolution0.9_dodge.pdf",device = "pdf")

p <- ggplot(CML@meta.data, aes(x=Type))
p <- p + geom_bar(aes(fill= factor(Phase)),  position="fill")
ggsave("Type_Phase_resolution0.9_dodge.pdf",device = "pdf")

pdf("CellCycle_RidgePlot.pdf", width = 16, height = 10)
p <- RidgePlot(CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(factor(CML@meta.data$Phase)))), group.by = "Phase", features = c("MKI67","PCNA", "TOP2A", "HMGB2","CENPF","UBE2C","TYMS","NUSAP1","CDC20","TPX2","GINS2","UHRF1","HELLS","NASP"))
print(p)
dev.off()

pdf("UMAP_plot_by_Phase.pdf",width = 5, height = 4)
p1 <- DimPlot(CML, cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(factor(CML@meta.data$Phase)))), pt.size = 0.5, label = FALSE, reduction = "umap.harmony", group.by = "Phase", label.size = 4)
print(p1)
dev.off()

p <- VlnPlot(CML, features = "G2M.Score", pt.size = 0, group.by = "cluster", split.by = "Phase")
ggsave(plot=p, "CellCycle_VlnPlot_splitPhase.pdf",width = 12, height = 4)



# Draw violin plots of cell cycle genes for myeloid and lymphoid clusters
p <- VlnPlot(sceM, group.by = "cluster",features =s_genes,  pt.size = 0, ncol = 4,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceM@meta.data$cluster))))
ggsave(plot=p, file="vlnSCEm_s_genes_Data.pdf",width =12, height =16)
p <- VlnPlot(sceM, group.by = "cluster",features =g2m_genes,  pt.size = 0, ncol = 4,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceM@meta.data$cluster))))
ggsave(plot=p, file="vlnSCEm_g2m_genes_Data.pdf",width =12, height =20)

p <- VlnPlot(sceB, group.by = "cluster",features =s_genes,  pt.size = 0, ncol = 4,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceB@meta.data$cluster))))
ggsave(plot=p, file="vlnSCEb_s_genes_Data.pdf",width =14, height =18)
p <- VlnPlot(sceB, group.by = "cluster",features =g2m_genes,  pt.size = 0, ncol = 4,  layer = "data",
    cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(sceB@meta.data$cluster))))
ggsave(plot=p, file="vlnSCEb_g2m_genes_Data.pdf",width =14, height =22)



### Draw heatmap plots of cell cycle genes for myeloid and lymphoid clusters
library(pheatmap)
library(viridis)
sgenes <- c("NASP","SLBP", "RPA2","TYMS","PCNA","FEN1","GINS2","HELLS","HMGB2","TUBB4B","TMPO", "CENPF","CKS1B", "MKI67", "BIRC5","TOP2A", "CDC20", "CENPE","AURKB")
sceM.AVG <- AverageExpression(object = sceM, assays= "RNA", layer = "counts", group.by = "cluster", return.seurat = TRUE)
cluster_mean_counts <- GetAssayData(sceM.AVG[["RNA"]], layer = "counts")
clusterTop10_mean_counts <- as.matrix(cluster_mean_counts[sgenes,])
p <- pheatmap(clusterTop10_mean_counts, color=colorRampPalette(c("navy", "white", "firebrick3"))(50), border_color="white", scale="row", cluster_rows=F, cluster_cols=F, fontsize_row = 6)
ggsave(plot=p, file="heatmap_sceM_sgenes.pdf",width =4, height =4)


sceB.AVG <- AverageExpression(object = sceB, assays= "RNA", layer = "counts", group.by = "cluster", return.seurat = TRUE)
cluster_mean_counts <- GetAssayData(sceB.AVG[["RNA"]], layer = "counts")
clusterTop10_mean_counts <- as.matrix(cluster_mean_counts[sgenes,])
p <- pheatmap(clusterTop10_mean_counts, color=colorRampPalette(c("navy", "white", "firebrick3"))(50), border_color="white", scale="row", cluster_rows=F, cluster_cols=F, fontsize_row = 6)
ggsave(plot=p, file="heatmap_sceB_sgenes.pdf",width =4, height =4)





# ----- 10.2 Gene signature scoring related to apoptosis, p53 signaling, proliferation, B-cell receptor (BCR) signaling, and IL6 and IL12 production
# ------------------------------------------------


#####  Gene signature scoring using the AddModuleScore  function ########

####  ------ Load the signature genes  ---------    #####

# Get all gene names from the Seurat object `CML`
all.genes <- rownames(CML)

# --- Long-term HSC (LT-HSC) signature ---
# Read the gene list from a tab-separated file
LT_HSC = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/LT-HSC.txt"),sep='\t',header=T,stringsAsFactors =F)
# Keep only genes that are present in the Seurat object
LT_HSC <- LT_HSC[LT_HSC$gene %in% all.genes,]
# Convert to a list (AddModuleScore expects a list of gene sets)
LT_HSC <- as.list(LT_HSC)


# --- Short-term HSC (ST-HSC) signature ---
ST_HSC = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/ST-HSC.txt"),sep='\t',header=T,stringsAsFactors =F)
ST_HSC <- ST_HSC[ST_HSC$gene %in% all.genes,]
ST_HSC <- as.list(ST_HSC)


# --- Quiescent HSC signature ---
Quiescent_HSC = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/Quiescent-HSC.txt"),sep='\t',header=T,stringsAsFactors =F)
Quiescent_HSC <- Quiescent_HSC[Quiescent_HSC$gene %in% all.genes,]
Quiescent_HSC <- as.list(Quiescent_HSC)


# --- Apoptosis signature ---
Apoptosis = read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/Apoptosis.txt"),sep='\t',header=T,stringsAsFactors =F)
Apoptosis <- Apoptosis[Apoptosis$gene %in% all.genes,]
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



# ============================================================
# ------ Calculate the signature scores ------
# ============================================================
# AddModuleScore computes a module score for each cell by comparing
# the average expression of the signature genes to a set of control genes.
# Parameters:
#   object  : Seurat object
#   features: list of gene sets (each element is a character vector of genes)
#   nbin    : number of bins for controlling for average expression
#   ctrl    : number of control genes per bin
#   name    : prefix for the metadata column name; the actual column will be
#             name followed by a number (e.g., "Apoptosis1")

####  ------ Caculate the signature score  ---------    #####
CML <- AddModuleScore(object = CML, features = Apoptosis, nbin = 25, ctrl = 20, name = 'Apoptosis')
CML <- AddModuleScore(object = CML, features = P53, nbin = 25, ctrl = 20, name = 'P53')
CML <- AddModuleScore(object = CML, features = Proliferation, nbin = 25, ctrl = 20, name = 'Proliferation')

CML <- AddModuleScore(object = CML, features = IL12, nbin = 25, ctrl = 20, name = 'IL12')
CML <- AddModuleScore(object = CML, features = BCR, nbin = 25, ctrl = 20, name = 'BCR')
CML <- AddModuleScore(object = CML, features = Tumor_inflammation, nbin = 25, ctrl = 20, name = 'Tumor_inflammation')
CML <- AddModuleScore(object = CML, features = IL6, nbin = 25, ctrl = 20, name = 'IL6_production')

CML <- AddModuleScore(object = CML, features = LT_HSC, nbin = 25, ctrl = 20, name = 'LT_HSC')
CML <- AddModuleScore(object = CML, features = ST_HSC, nbin = 25, ctrl = 20, name = 'ST_HSC')
CML <- AddModuleScore(object = CML, features = Quiescent_HSC, nbin = 25, ctrl = 20, name = 'Quiescent_HSC')



# ============================================================
# ------ Aggregate multiple score columns into a single score ------
# ============================================================
# For each signature, the AddModuleScore function may create multiple
# columns if the gene set is split? Actually it creates one column per
# call, named with the given name plus a suffix (e.g., "BCR1").
# The following code selects columns matching the signature name and
# calculates the row-wise mean to obtain a single score per cell.

# --- BCR score ---
metaBCR <- CML@meta.data %>% select(matches("BCR"))
BCRscore <- apply(metaBCR, 1, mean)
CML@meta.data$BCRscores <- BCRscore

# --- IL12 score ---
metaIL12 <- CML@meta.data %>% select(matches("IL12"))
IL12Score <- apply(metaIL12, 1, mean)
CML@meta.data$IL12Scores <- IL12Score

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

# --- Quiescent HSC score ---
metaQuiescent_HSC <- CML@meta.data %>% select(matches("Quiescent_HSC"))
Quiescent_HSCscore <- apply(metaQuiescent_HSC, 1, mean)
CML@meta.data$Quiescent_HSCscores <- Quiescent_HSCscore

# --- ST HSC score ---
metaST_HSC <- CML@meta.data %>% select(matches("ST_HSC"))
ST_HSCscore <- apply(metaST_HSC, 1, mean)
CML@meta.data$ST_HSCscores <- ST_HSCscore

# --- LT HSC score ---
metaLT_HSC <- CML@meta.data %>% select(matches("LT_HSC"))
LT_HSCscore <- apply(metaLT_HSC, 1, mean)
CML@meta.data$LT_HSCscores <- LT_HSCscore


# ============================================================
# ------ Visualize signature scores by cell type ------
# ============================================================

# Set the active identity to the "Type" column (cell type annotation)
type <- CML@meta.data$Type
names(type) <- rownames(CML@meta.data)
CML@active.ident <- type

# Collect all columns that contain "Score" or "score" in their names
signatureScore1 = colnames(CML@meta.data)[grep("Score", colnames(CML@meta.data))]
signatureScore2 = colnames(CML@meta.data)[grep("score", colnames(CML@meta.data))]
signatureScores <- c(signatureScore1, signatureScore2)

# Create violin plots for all signature scores, grouped by cell type
pdf("signatureScores_VlnPlot_byType_new.pdf", width = 6, height = 8)
p <- VlnPlot(CML, features = signatureScores, group.by = "Type", cols=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(CML@meta.data$Type))), pt.size = 0, ncol = 5)
print(p)
dev.off()




# Extract all metadata columns containing "Score" or "score"
signatureScore <- as.matrix(CML@meta.data %>% select(matches(c("Score","score"))))

# Transpose the matrix so that rows are signatures and columns are cells
meta_sigt <- t(signatureScore)

# Create a data frame with signature names as a column
score <- data.frame(row.names(meta_sigt))
colnames(score) <- "score"
meta_sig <- cbind(score,meta_sigt)

# Melt the data into long format: Signature, barcode, score
mdf <- reshape2::melt(meta_sig, id.vars="score")
colnames(mdf)=c("Signature","barcode","score")

# Extract cell barcodes and selected metadata (cluster, Type, Sample)
barcode <- row.names(CML@meta.data)
features <- CML@meta.data %>% select(matches(c("cluster","Type","Sample")))
clin <- cbind(barcode,features)

# Merge signature scores with clinical/metadata information
mdf_clin <- merge(mdf,clin,by="barcode",all.x=T)

# Define a two-color palette for the plot
s_color <- c("DodgerBlue","Salmon")


# Load plotting packages
library(tidyverse)
library(tidyplots)


# ============================================================
# Boxplot: Signature scores across all signatures, colored by Type
# ============================================================

p <- mdf_clin |>
tidyplot(x = Signature, y = score, color = Type) |>
add_boxplot(show_outliers = F, dodge_width = 0.8, box_width = 0.6) |>
adjust_y_axis_title("Signature Scores") |>
 adjust_size(width = 160,height = 60) |>
  adjust_colors(new_colors = s_color) |>
  add_test_asterisks(hide_info = TRUE, method = "wilcox_test") |>
  adjust_x_axis(rotate_labels = 90)

pdf("boxplot_signature_byType.pdf", width = 12, height = 4)
p
dev.off()


# ============================================================
# Boxplot: Signature scores by Celltype, split by Signature
# ============================================================

p <- mdf_clin |>
tidyplot(x = Celltype, y = score, color = Type, width = 160,height = 100) |>
add_boxplot(show_outliers = F, dodge_width = 0.8, box_width = 0.6) |>
adjust_y_axis_title("Signature Scores") |>
  adjust_colors(new_colors = s_color) |>
  add_test_asterisks(hide_info = TRUE, method = "wilcox_test") |>
split_plot(by = Signature, ncol = 5)

pdf("boxplot_signature_Celltype_byType.pdf", width =36, height =16)
p
dev.off()


# ============================================================
# Plot signature scores for a specific cell type or cluster
# ============================================================

# --- Only HSPC cells ---
p <- mdf_clin |>
filter(Celltype %in% c("HSPC")) |> 
tidyplot(x = Signature, y = score, color = Type, width = 160,height = 60) |>
add_boxplot(show_outliers = F, dodge_width = 0.8, box_width = 0.6) |>
adjust_y_axis_title("Signature Scores") |>
 adjust_x_axis(rotate_labels = 90) |>
  adjust_colors(new_colors = s_color) |>
  add_test_asterisks(hide_info = TRUE, method = "wilcox_test")

pdf("boxplot_signature_HSPC_byType.pdf", width = 12, height = 4)
p
dev.off()


# --- Only cluster 11 ---
p <- mdf_clin |>
filter(cluster %in% c("11")) |> 
tidyplot(x = Signature, y = score, color = Type, width = 160,height = 60) |>
add_boxplot(show_outliers = F, dodge_width = 0.8, box_width = 0.6) |>
adjust_y_axis_title("Signature Scores") |>
 adjust_x_axis(rotate_labels = 90) |>
  adjust_colors(new_colors = s_color) |>
  add_test_asterisks(hide_info = TRUE, method = "wilcox_test")

pdf("boxplot_signature_cluster11_byType.pdf", width = 12, height = 4)
p
dev.off()



# ============================================================
# Plot selected signatures for HSPC cells
# ============================================================

# --- HSC quiescence / LT-HSC / ST-HSC signatures ---

p <- mdf_clin |>
dplyr::filter(Celltype %in% c("HSPC")) |>
dplyr::filter(Signature %in% c("Quiescent_HSCscores","LT_HSCscores","ST_HSCscores")) |>
tidyplot(x = Celltype, y = score, color = Type, width = 160,height = 100) |>
add_boxplot(show_outliers = F, dodge_width = 0.8, box_width = 0.6) |>
adjust_y_axis_title("Signature Scores") |>
  adjust_colors(new_colors = s_color) |>
  add_test_asterisks(hide_info = TRUE, method = "wilcox_test") |>
 split_plot(by = Signature, ncol = 3)

pdf("boxplot_signature_HSCscore_byType.pdf", width =8, height =4)
p
dev.off()


# ============================================================
# Proliferation signatures for HSPC cells
# ============================================================

 p <- mdf_clin |>
dplyr::filter(Celltype %in% c("HSPC")) |>
dplyr::filter(Signature %in% c("S.Score","G2M.Score","ProliferationScores")) |>
tidyplot(x = Celltype, y = score, color = Type, width = 160,height = 100) |>
add_boxplot(show_outliers = F, dodge_width = 0.8, box_width = 0.6) |>
adjust_y_axis_title("Signature Scores") |>
  adjust_colors(new_colors = s_color) |>
  add_test_asterisks(hide_info = TRUE, method = "wilcox_test") |>
 split_plot(by = Signature, ncol = 3)

pdf("boxplot_signature_HSCproliferation_byType.pdf", width =8, height =4)
p
dev.off()


#######






# ============================================================
# Step 11: Pseudotime trajectory analysis using Monocle 2
# ============================================================
# This script performs pseudotime trajectory inference on a subset
# of cells using the Monocle 2 package. It includes data extraction,
# quality control, dimension reduction, trajectory construction, and
# visualization of the results.

library(monocle)

# ============================================================
# Extract the subset of cells for trajectory analysis
# ============================================================
# Select cells belonging to the specified samples

cellbarcodes <- rownames(CML@meta.data)[CML@meta.data$Sample %in% c("P1_CP","P2_CP","P1_BLBC","P2_BLBC")]
Paired <- subset(CML,cells = cellbarcodes)

# ============================================================
# Extract the expression matrix from the subset
# ============================================================
# Use raw counts from the RNA assay

monocl_data <- GetAssayData(Paired, assay = "RNA", slot = "counts")
monocl_data <- as.matrix(monocl_data)

# Remove genes with zero expression across all cells
monocl_data <- monocl_data[apply(monocl_data,1,sum) > 0,]


# ============================================================
# Prepare feature annotation data
# ============================================================
# Create a data frame with gene IDs and short names
feature_ann <- data.frame(gene_id=rownames(monocl_data),gene_short_name=rownames(monocl_data))
rownames(feature_ann) <- rownames(monocl_data)
feature_data <- new("AnnotatedDataFrame", data = feature_ann)

# ============================================================
# Prepare phenotype annotation data
# ============================================================
# Extract cell metadata
sample_ann <- Paired@meta.data

# Create a cluster label with prefix "Cluster"
cluster <- as.character(sample_ann$Cluster)
clusterN  <- factor(sapply(cluster,function(x) paste("Cluster",x,sep = "",collapse = NULL)))

# Build a simplified phenotype data frame
sample_ann <- data.frame(Cellbarcode=rownames(sample_ann),Celltype=sample_ann$Celltype,Cluster=sample_ann$cluster,Sample=sample_ann$Sample)
rownames(sample_ann) <- sample_ann$Cellbarcode
pheno_data <- new("AnnotatedDataFrame", data = sample_ann)



# ============================================================
# Construct the CellDataSet object
# ============================================================
# Important: choose the correct expressionFamily based on data type.
# - For FPKM values: use tobit()
# - For UMIs or transcript counts: use negbinomial.size()
# Here we assume the data are UMI counts and use negbinomial.size().
#
# CellDataSet structure:
#   exprs       : numeric matrix (rows = genes, columns = cells)
#   phenoData   : AnnotatedDataFrame (rows = cells, columns = cell attributes)
#   featureData : AnnotatedDataFrame (rows = genes, columns = gene attributes)
#
# NOTE: Do NOT normalize UMI data before creating the CellDataSet.
#       Monocle will perform all necessary normalization internally.

monocl_object <- newCellDataSet(monocl_data, phenoData = pheno_data,featureData = feature_data,expressionFamily=negbinomial.size())


# Rename for convenience
HSMM <- monocl_object


# ============================================================
# Normalization and dispersion estimation
# ============================================================
# Estimate size factors and dispersions for downstream analysis
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM, cores=10)


# Inspect the phenotype data (size factors added)；
head(pData(HSMM))

# Inspect the dispersion table
head(dispersionTable(HSMM))


# ============================================================
# Quality control: filter low-quality cells and genes
# ============================================================
# Detect genes expressed in at least min_expr counts

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

# Filter genes expressed in at least 10 cells
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
HSMM <- HSMM[expressed_genes,]


# Select genes for clustering based on mean expression
#Genes that aren’t highly expressed enough will not be used for clustering, since they may not give meaningful signal and would only add noise.
disp_table <- dispersionTable(HSMM)
clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, clustering_genes$gene_id)


# ============================================================
# Dimensionality reduction and clustering (optional)
# ============================================================
# Reduce dimensions using UMAP and cluster cells with Louvain

HSMM <- reduceDimension(HSMM, num_dim = 40, reduction_method = 'umap')
HSMM <- clusterCells(HSMM, method = "louvain")

# ============================================================
# 9. Select genes for trajectory inference
# ============================================================
# Option A: Use differential expression analysis based on existing clusters
diff_genes <- differentialGeneTest(HSMM, fullModelFormulaStr = "~ Cluster",
                                   cores = 10)


# Select top 3000 differentially expressed genes
ordering_genes <- row.names(subset(diff_genes, qval < 1e-3))[order(diff_genes$qval)][1:3000]
HSMM <- setOrderingFilter(HSMM, ordering_genes)


# ============================================================
# Trajectory inference with DDRTree
# ============================================================
# Reduce dimensions with DDRTree and order cells in pseudotime
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

plot_cell_trajectory(cds, color_by = "cell_type", cell_size = 1) +
  scale_color_brewer(name = "cell type", type = "qual", palette = "Set2")




# ============================================================
# Build trajectory (alternative approach)
# ============================================================
# Use variable genes identified by Seurat's FindVariableFeatures

vst_plot_data <- FindVariableFeatures(object=Paired@assays$RNA@layers$data,selection.method = "vst")
vst_plot_data <- vst_plot_data[vst_plot_data$vst.mean >= 0.1 & vst_plot_data$vst.mean <= 8 & vst_plot_data$vst.variance >= 1,]
vst_plot_data <- vst_plot_data[order(vst_plot_data$vst.variance,decreasing=TRUE),]
vst_plot_data <- vst_plot_data[1:1000,]

# Mark genes for ordering
HSMM <- setOrderingFilter(HSMM, rownames(vst_plot_data))

# Infer trajectory using DDRTree
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)


# ============================================================
# Visualization of trajectory results
# ============================================================

# --- Trajectory colored by State ---
pdf('pseudotime_trajectory_State2_Btype.pdf', height=5)
pp <- plot_cell_trajectory(HSMM,color_by = "State", show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_manual(values=brewer.pal(8,'Set2'))
print(pp)
dev.off()

# --- Trajectory colored by State, faceted by State ---
pdf('pseudotime_trajectory_StateSingle.pdf', width=12, height=8)
pp <- plot_cell_trajectory(HSMM,color_by = "State", show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_manual(values=brewer.pal(8,'Set2')) +facet_wrap(~State,nrow=2,ncol=4)
print(pp)
dev.off()

# --- Trajectory colored by Cluster ---
pdf('pseudotime_trajectory_cluster.pdf', height=5)
pp <- plot_cell_trajectory(HSMM,color_by = "Cluster",show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_manual(values=brewer.pal(12,'Paired'))
print(pp)
dev.off()

# --- Trajectory colored by Cluster, faceted by Cluster ---
pdf('pseudotime_trajectory_clusterSingle.pdf', width=12, height=10)
pp <- plot_cell_trajectory(HSMM,color_by = "Cluster",show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_manual(values=brewer.pal(12,'Paired'))+facet_wrap(~Cluster,nrow=3,ncol=4)
print(pp)
dev.off()

# --- Trajectory colored by Cluster, faceted by Sample ---
pdf('pseudotime_trajectory_clusterSample.pdf', width=8, height=4)
pp <- plot_cell_trajectory(HSMM,color_by = "Cluster",show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_manual(values=brewer.pal(12,'Paired'))+facet_wrap(~Sample)
print(pp)
dev.off()


# --- Trajectory colored by Sample ---
color <- c("blue","red")
pdf('pseudotime_trajectory_sample.pdf', height=5)
pp <- plot_cell_trajectory(HSMM,color_by = "Sample",show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_manual(values=color)
print(pp)
dev.off()

# --- Trajectory colored by Pseudotime, faceted by type ---
pdf('pseudotime_trajectory_pseudotime_Btype.pdf', height=5)
pp <- plot_cell_trajectory(HSMM,color_by = "Pseudotime",show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_gradient(high="red", low="lightgrey")+facet_wrap(~Type)
print(pp)
dev.off()

# --- Trajectory colored by Pseudotime, faceted by Sample ---
pdf('pseudotime_trajectory_pseudotimeSample.pdf', height=5)
pp <- plot_cell_trajectory(HSMM,color_by = "Pseudotime",show_branch_points=FALSE,show_backbone = T,show_tree = T,cell_link_size=1,cell_size=1)+scale_color_gradient(high="blue", low="lightgrey")+facet_wrap(~Sample)
print(pp)
dev.off()





#######

# ============================================================
# Step 12: Add BCR::ABL1 fusion to Seurat metadata
# ============================================================

# Read the fusion information table
# The file contains per-cell fusion status (e.g., "Fusion" column) and
# other metadata. `barcode` is expected to be one of the columns.

FusionInfo <- read.table(paste0("/public/", me, "projects/CMLblast/figure/CMLharmony/FusionInfo.txt"),header=TRUE, sep='\t',stringsAsFactors =TRUE)

# Set cell barcodes as row names so that AddMetaData can match them
# to the cells in the Seurat object.
rownames(FusionInfo) <- FusionInfo$barcode

# Add all columns from FusionInfo to the Seurat object's metadata.
# This includes the "Fusion" column indicating BCR-ABL1 fusion status.
CML <- AddMetaData(CML, metadata = FusionInfo)

# Inspect the updated metadata
head(CML@meta.data)

# Ensure the default assay is set to RNA for downstream plotting
DefaultAssay(CML) <- "RNA"


# ============================================================
# Visualize fusion status on the UMAP embedding
# ============================================================

# --- UMAP colored by fusion status ---
# `cols` defines colors for the two groups (e.g., Fusion-negative vs positive).
# `order = TRUE` ensures that cells with the fusion (red) are plotted on top.
# `reduction = "umap.harmony"` uses the Harmony-integrated UMAP.

p <- DimPlot(CML, cols=c("lightgrey", "red"), pt.size = 1.5,  order = TRUE, reduction = "umap.harmony", group.by = "Fusion")
ggsave(plot=p, file="DimPlot_CML_resolution0.9_Fusion_status_Data.pdf",width = 5, height = 4)

p <- DimPlot(CML, cols=c("lightgrey", "red"), pt.size = 1.5,  order = TRUE, reduction = "umap.harmony", group.by = "Fusion", split.by = "Sample")
ggsave(plot=p, file="DimPlot_CML_resolution0.9_Fusion_status_Data_splitSample.pdf",width = 20, height = 4)


