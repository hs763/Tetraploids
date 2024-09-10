cd /cephfs2/hannas/tetraploids/analysis
singularity shell --bind /cephfs2/hannas:/mnt /public/singularity/containers/nextflow/singlecell_latest.sif
cd /mnt/tetraploids/analysis/preprocessing
R

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(harmony)
library(future)

#reading in the data
path2data <- "/mnt/tetraploids/expdata/data_samptab/combined/all-sample/DGE_unfiltered/"
counts <- readMM(paste0(path2data, "count_matrix.mtx"))
metadata <- read.csv("/mnt/tetraploids/analysis/preprocessing/metadata_TBC.csv")
gene_info <- read.csv("/mnt/tetraploids/analysis/preprocessing/gene_info.csv")
barcodes <- read.csv("/mnt/tetraploids/analysis/preprocessing/barcodes.csv")

# gene_info <- all_genes
# duplicated_names <- gene_info$gene_name[duplicated(gene_info$gene_name) | duplicated(gene_info$gene_name, fromLast = TRUE)]
# gene_info$unique_name <- ifelse(gene_info$gene_name %in% duplicated_names, gene_info$gene_id, gene_info$gene_name)
# write.csv(gene_info, "/mnt/tetraploids/analysis/preprocessing/gene_info.csv", row.names = FALSE)

counts <- as(counts, "CsparseMatrix")
counts <- t(counts)
colnames(counts) <- barcodes$barcodes
rownames(counts) <- gene_info$unique_name 

so <- CreateSeuratObject(counts = counts,
                            assay = "RNA",
                            meta.data = metadata,
                            project = "tetraploids"
                          )
# An object of class Seurat 
# 63140 features across 3273556 samples within 1 assay 
# Active assay: RNA (63140 features, 0 variable features)
#  1 layer present: counts

saveRDS(so, "SeuratObject_raw.rds")

...

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
head(so@meta.data, 5)


#ploting QC plots 
nFeature_bySampleName <- VlnPlot(so, features = "nFeature_RNA", pt.size = 0, group.by = "sample_name", flip = TRUE, log = TRUE) 
nFeature_bySampleName <- nFeature_bySampleName + theme(
  axis.text.x = element_text(size = 8),
  axis.title.x = element_text(size = 12), 
  plot.title = element_text(size = 14),
  aspect.ratio = 1/3,
  legend.position = "none" 
) + labs(
  x = "Sample Name",  
  y = "log[nFeature]"    
)
ggsave("nFeature_bySampleName_violin.png", plot = nFeature_bySampleName, width = 10, height = 8, dpi = 300)

nCount_bySampleName <- VlnPlot(so, features = "nCount_RNA", pt.size = 0, group.by = "sample_name", flip = TRUE, log = TRUE) 
nCount_bySampleName <- nCount_bySampleName + theme(
  axis.text.x = element_text(size = 8),
  axis.title.x = element_text(size = 12), 
  plot.title = element_text(size = 14),
  aspect.ratio = 1/3,
  legend.position = "none" 
) + labs(
  x = "Sample Name",  
  y = "log[nCount]"    
)
ggsave("nCount_bySampleName_violin.png", plot = nCount_bySampleName, width = 10, height = 8, dpi = 300)

PerMT_bySampleName <- VlnPlot(so, features = "percent.mt", pt.size = 0, group.by = "sample_name", flip = TRUE, log = TRUE) 
PerMT_bySampleName <- PerMT_bySampleName + theme(
  axis.text.x = element_text(size = 8),
  axis.title.x = element_text(size = 12), 
  plot.title = element_text(size = 14),
  aspect.ratio = 1/3,
  legend.position = "none" 
) + labs(
  x = "Sample Name",  
  y = "log[PerMT]"    
)
ggsave("PerMT_bySampleName_violin.png", plot = PerMT_bySampleName, width = 10, height = 8, dpi = 300)


plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt", log = FALSE, jitter = TRUE, group.by = "sample_name")
ggsave("nCounts_vs_PerMT.png", plot = plot1, width = 10, height = 8, dpi = 300)
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample_name")
ggsave("nCounts_vs_nFeature.png", plot = plot2, width = 10, height = 8, dpi = 300)

so@meta.data <- so@meta.data %>%
  mutate(
    batch = ifelse(grepl("_2$", sample_name), NA, batch),
    day = ifelse(grepl("_2$", sample_name), NA, day)
  )

saveRDS(so, "SeuratObject_befroeSubset.rds")

#subseting by removing cells with insufficent feature and transcritp counts and too high mt. fractions and too high complexity. 
so@meta.data$complexity =  so@meta.data$nFeature_RNA/so@meta.data$nCount_RNA
summary(so@meta.data$complexity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00949 0.91667 0.95000 0.93189 0.98276 1.00000

# so 
# An object of class Seurat 
# 63140 features across 3273556 samples within 1 assay 
# Active assay: RNA (63140 features, 0 variable features)
#  1 layer present: counts

###########Subset1##############
so1 <- subset(so, subset = nFeature_RNA > 200 &nFeature_RNA < 18000 & nCount_RNA > 200 & nCount_RNA < 800000 & percent.mt < 20 & complexity < 0.98)
saveRDS(so1, "SeuratObject_Subset1.rds")
# An object of class Seurat 
# 63140 features across 267667 samples within 1 assay 
# Active assay: RNA (63140 features, 0 variable features)
#  1 layer present: counts

#normalisation 
so1 <- NormalizeData(so1, normalization.method = "LogNormalize", scale.factor = 10000)

#identifying highly variable features 
so1 <- FindVariableFeatures(so1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(so1), 10)

top10 
# [1] "CNTN5"           "FLI1"            "COL3A1"          "POSTN"          
#  [5] "FLT1"            "GRID2"           "CDC20B"          "ENSG00000288921"
#  [9] "PTPRT"           "PCAT14" 

plot1 <- FeatureScatter(so1, feature1 = "nCount_RNA", feature2 = "percent.mt", log = FALSE, jitter = TRUE, group.by = "sample_name")
ggsave("nCounts_vs_PerMT_Subset1.png", plot = plot1, width = 10, height = 8, dpi = 300)
plot2 <- FeatureScatter(so1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample_name")
ggsave("nCounts_vs_nFeature_Subset1.png", plot = plot2, width = 10, height = 8, dpi = 300)

plot1 <- VariableFeaturePlot(so1)
ggsave("Unlabled_VariableFeatures_Subset1.png", plot = plot1, width = 10, height = 8, dpi = 300)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("Labled_VariableFeatures_Subset1.png", plot = plot2, width = 10, height = 8, dpi = 300)

#########Subset2##########
so2 <- subset(so, subset = nFeature_RNA > 200 &nFeature_RNA < 18000 & nCount_RNA > 500 & nCount_RNA < 800000 & percent.mt < 20 & complexity < 0.98)
saveRDS(so2, "SeuratObject_Subset2.rds")
# An object of class Seurat 
# 63140 features across 154554 samples within 1 assay 
# Active assay: RNA (63140 features, 0 variable features)
#  1 layer present: counts

#normalisation 
so2 <- NormalizeData(so2, normalization.method = "LogNormalize", scale.factor = 10000)

#identifying highly variable features 
so2 <- FindVariableFeatures(so2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(so2), 10)

top10 
#  [1] "CNTN5"           "FLI1"            "COL3A1"          "POSTN"          
#  [5] "GRID2"           "FLT1"            "CDC20B"          "ENSG00000288921"
#  [9] "PTPRT"           "PCAT14" 

plot1 <- FeatureScatter(so2, feature1 = "nCount_RNA", feature2 = "percent.mt", log = FALSE, jitter = TRUE, group.by = "sample_name")
ggsave("nCounts_vs_PerMT_Subset2.png", plot = plot1, width = 10, height = 8, dpi = 300)
plot2 <- FeatureScatter(so2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample_name")
ggsave("nCounts_vs_nFeature_Subset2.png", plot = plot2, width = 10, height = 8, dpi = 300)

plot1 <- VariableFeaturePlot(so2)
ggsave("Unlabled_VariableFeatures_Subset2.png", plot = plot1, width = 10, height = 8, dpi = 300)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("Labled_VariableFeatures_Subset2.png", plot = plot2, width = 10, height = 8, dpi = 300)

#the top10 for susbet2 is the same as subset1 so i will keep mroe cells in the analysis and continue wiht susbet1. 
so <- so1

all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)
saveRDS(so, "SeuratObject_Norm_Scaled_Subset1.rds")

...

plan("multicore", workers = 4)
so <- RunPCA(so, features = VariableFeatures(object = so))
saveRDS(so, "SeuratObject_Norm_Scaled_PCA_Subset1.rds")
print(so[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(so, dims = 1:5, reduction = "pca", ncol = 5)
pca_SampleName <- DimPlot(so, reduction = "pca", pt.size = 0.5,group.by = "sample_name")
ggsave("PCA_bySampleName_Subset1.png", plot = pca_SampleName, width = 10, height = 8, dpi = 300)

pca_batch <- DimPlot(so, reduction = "pca", pt.size = 0.5,group.by = "batch")
ggsave("PCA_byBatch_Subset1.png", plot = pca_batch, width = 10, height = 8, dpi = 300)

pca_batch <- DimPlot(so, reduction = "pca", pt.size = 0.5,group.by = "day")
ggsave("PCA_byDay_Subset1.png", plot = pca_batch, width = 10, height = 8, dpi = 300)

DimHeatmap(so, dims = 1, cells = 500, balanced = TRUE)
DimHeatmapPCA15 <- DimHeatmap(so, dims = 1:15, cells = 500, balanced = TRUE)
ggsave("DimHeatMap_first_15PCAs.png", plot = DimHeatmapPCA15,  width = 8, height = 12, dpi = 300)

elbow <- ElbowPlot(so)
ggsave("elbow_plot_Subset1.png", plot = elbow,  width = 10, height = 8, dpi = 300)



#Choosing the nuber of PCAs to contimue with to UMAP. 
so <- FindNeighbors(so, dims = 1:8)
so <- FindClusters(so, resolution = 0.2)
head(Idents(so), 8)

so <- RunUMAP(so, dims = 1:8)
umap <- DimPlot(so, reduction = "umap", group.by = "sample_name")
saveRDS(so, "SeuratObject_Norm_UMAP_Subset1.rds")

umap <- DimPlot(so, reduction = "umap", group.by = "day")
ggsave("umap_top5PCAs_byDay.png", plot = umap,  width = 10, height = 8, dpi = 300)

umap <- DimPlot(mso, reduction = "umap", group.by = "batch")
ggsave("umap_top5PCAs_bybatch.png", plot = umap,  width = 10, height = 8, dpi = 300)

umap <- DimPlot(mso, reduction = "umap", group.by = "RNA_snn_res.0.2")
ggsave("umap_top5PCAs_byCluster.png", plot = umap,  width = 10, height = 8, dpi = 300)



#clsuter anlaysis 
so.markers <- FindAllMarkers(so, only.pos = TRUE)
so.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

saveRDS(so.markers, "mso.marker.rds")

VlnPlot(so, features = c("")) 
VlnPlot(so, features = c(""), slot = "counts", log = TRUE)

FeaturePlot(so, features = c())

so.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(so, features = top10$gene) + NoLegend()








