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

so <- readRDS("SeuratObject_Norm_Scaled_Subset1.rds")

plan("multicore", workers = 4) 
so <- RunPCA(so, features = VariableFeatures(object = so))
saveRDS(so, "SeuratObject_Norm_Scaled_PCA_Subset1.rds")
print(so[["pca"]], dims = 1:5, nfeatures = 5)
# PC_ 1 
# Positive:  CACNA1A, NEAT1, KALRN, SDK1, NAV2 
# Negative:  S100A8, ENSG00000278911, ENSG00000273988, ENSG00000254575, OR5D15P 
# PC_ 2 
# Positive:  WNT5B, MECOM, PLCH1, GPC3, PLXDC2 
# Negative:  TRHDE, THSD7B, EBF3, NXPH1, ELAVL3 
# PC_ 3 
# Positive:  WNT5B, RMST, MIAT, NAV3, LINC01414 
# Negative:  GRHL2, CDH1, HAPLN1, ESRP1, ATP8B1 
# PC_ 4 
# Positive:  EPHA3, CREB5, LHX2, NPAS3, DLGAP2 
# Negative:  LINC02315, MIR924HG, TMEM132C, ENSG00000286099, COLEC12 
# PC_ 5 
# Positive:  EDNRA, LINC00511, ITGA4, SOX10, KANK4 
# Negative:  RBFOX1, SHISA9, ENSG00000262801, LINC00278, GRHL2 

VizDimLoadings(so, dims = 1:5, reduction = "pca", ncol = 5)
pca_Sample <- DimPlot(so, reduction = "pca", pt.size = 0.5,group.by = "sample_name")
ggsave("PCA_bySample_Subset1.png", plot = pca_SampleName, width = 10, height = 8, dpi = 300)

pca_batch <- DimPlot(so, reduction = "pca", pt.size = 0.5,group.by = "batch")
ggsave("PCA_byBatch_Subset1.png", plot = pca_batch, width = 10, height = 8, dpi = 300)

pca_batch <- DimPlot(so, reduction = "pca", pt.size = 0.5,group.by = "day")
ggsave("PCA_byDay_Subset1.png", plot = pca_batch, width = 10, height = 8, dpi = 300)

DimHeatmap(so, dims = 1, cells = 500, balanced = TRUE)
DimHeatmapPCA10 <- DimHeatmap(so, dims = 1:10, cells = 500, balanced = TRUE)
ggsave("DimHeatMap_first_10PCAs.png", plot = DimHeatmapPCA15,  width = 8, height = 12, dpi = 300)

elbow <- ElbowPlot(so)
ggsave("elbow_plot_Subset1.png", plot = elbow,  width = 10, height = 8, dpi = 300)



#Choosing the nuber of PCAs to contimue with to UMAP. 
options(future.globals.maxSize = 200 * 1024^3)
so <- FindNeighbors(so, dims = 1:9)
so <- FindClusters(so, resolution = 0.1, future.seed=TRUE)
# UNRELIABLE VALUE: One of the 'future.apply' iterations ('future_lapply-1') unexpectedly generated random numbers without declaring so. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'future.seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via the L'Ecuyer-CMRG method. To disable this check, use 'future.seed = NULL', or set option 'future.rng.onMisuse' to "ignore".
head(Idents(so),9)

so <- RunUMAP(so, dims = 1:9)
umap <- DimPlot(so, reduction = "umap", group.by = "sample_name")
ggsave("umap_top9PCAs_bySample.png", plot = umap,  width = 10, height = 8, dpi = 300)

saveRDS(so, "SeuratObject_Norm_UMAP_Subset1.rds")

umap_day <- DimPlot(so, reduction = "umap", group.by = "day")
ggsave("umap_top9PCAs_byDay.png", plot = umap_day,  width = 10, height = 8, dpi = 300)

umap_batch <- DimPlot(so, reduction = "umap", group.by = "batch")
ggsave("umap_top9PCAs_bybatch.png", plot = umap_batch,  width = 10, height = 8, dpi = 300)

umap_clust <- DimPlot(so, reduction = "umap", group.by = "seurat_clusters")
ggsave("umap_top9PCAs_byCluster.png", plot = umap_clust,  width = 10, height = 8, dpi = 300)

umap_species <- DimPlot(so, reduction = "umap", group.by = "species")
ggsave("umap_top9PCAs_bySpecies.png", plot = umap_species,  width = 10, height = 8, dpi = 300)


#integrationg wiht haarmony 
V <- Embeddings(so, reduction = "pca")
meta_data <- so@meta.data
harmony_embeddings_species <- harmony::RunHarmony(V, meta_data, 'species', verbose=FALSE, max_iter = 50)
harmony_embeddings_sample <- harmony::RunHarmony(V, meta_data, 'sample_name', verbose=FALSE)

so[['harmony.species']] <- CreateDimReducObject(
  embeddings = harmony_embeddings_species,
  key = "harmony_",   
  assay = DefaultAssay(so)
)

so[['harmony.sample']] <- CreateDimReducObject(
  embeddings = harmony_embeddings_sample,
  key = "harmony_",   
  assay = DefaultAssay(so)
)

so <- FindNeighbors(so, dims = 1:10, reduction = "harmony.species")
so <- FindClusters(so, resolution = 0.1, reduction = "harmony.batch")

so <- RunUMAP(so, dims = 1:10, reduction = "harmony.species", reduction.name = "umap.harmony.species")
umap_species <- DimPlot(so, reduction = "umap.harmony.species", group.by = "day")
ggsave("umap_harmony_species_10PCs_byday.png", plot = umap_species,  width = 10, height = 8, dpi = 300)

so <- RunUMAP(so, dims = 1:10, reduction = "harmony.sample", reduction.name = "umap.harmony.sample")
umap_sample <- DimPlot(so, reduction = "umap.harmony.sample", group.by = "day")
ggsave("umap_harmony_sample_10PCs_byDay.png", plot = umap_sample,  width = 10, height = 8, dpi = 300)

saveRDS(so,"SeuratObject_subset_harmony.rds")

#what is drving the variation? 
plot <- FeaturePlot(so, features = "nCount_RNA", 
                    reduction = "umap.harmony.species", 
                    pt.size = 1) +
  scale_color_viridis_c(trans = "log1p") + 
  labs(title = "UMAP with Log-scaled Counts")
ggsave("umap_harmony_species_10PCs_lognCounts.png", plot = plot,  width = 10, height = 8, dpi = 300)

plot <- FeaturePlot(so, features = "nFeature_RNA", 
                    reduction = "umap.harmony.species", 
                    pt.size = 1) +
  scale_color_viridis_c(trans = "log1p") + 
  labs(title = "UMAP with Log-scaled Features")
ggsave("umap_harmony_species_10PCs_lognFeature.png", plot = plot,  width = 10, height = 8, dpi = 300)

plot <- FeaturePlot(so, reduction = "umap.harmony.species", features = "percent.mt") +
  scale_color_viridis_c() + 
  labs(title = "UMAP with Log-scaled Percent MT")
ggsave("umap_harmony_species_10PCs_percent.mt.png", plot = plot,  width = 10, height = 8, dpi = 300)

plot <- FeaturePlot(so, features = "nCount_RNA", 
                    reduction = "umap", 
                    pt.size = 1) +
  scale_color_viridis_c(trans = "log1p") + 
  labs(title = "UMAP with Log-scaled Counts")
ggsave("umap_9PCs_lognCounts.png", plot = plot,  width = 10, height = 8, dpi = 300)

plot <- FeaturePlot(so, features = "nFeature_RNA", 
                    reduction = "umap", 
                    pt.size = 1) +
  scale_color_viridis_c(trans = "log1p") + 
  labs(title = "UMAP with Log-scaled Features")
ggsave("umap_9PCs_lognFeature.png", plot = plot,  width = 10, height = 8, dpi = 300)

plot <- FeaturePlot(so, reduction = "umap", features = "percent.mt") +
  scale_color_viridis_c() + 
  labs(title = "UMAP with Log-scaled Percent MT")
ggsave("umap_9PCs_percent.mt.png", plot = plot,  width = 10, height = 8, dpi = 300)


#clsuter anlaysis 
so.markers <- FindAllMarkers(so, only.pos = TRUE)
so.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

saveRDS(so.markers, "so.marker.rds")
VlnPlot(so, features = c("")) 
VlnPlot(so, features = c(""), slot = "counts", log = TRUE)

#exit form pluripotency (ectodermal/neural fates)
FeaturePlot(so, features = c("NANOG", "OCT4", "SOX2", "PAX6", "SOX1"),reduction = "umap.harmony.sample", cols = c("blue", "red"))
#neuroepithelial-toreadial-glia transition
FeaturePlot(so, features = c("ZEB2",))
#Pan-RG
FeaturePlot(so, Features = c("SOX2","VIM","NES","GLI3","FABP7", "HES1", "GFAP"))
#early radial glia 
FeaturePlot(so, features = c("LIX1", "HMGA2"))
#apical radial glia 
FeaturePlot(so, features = c("FBXO32","PROM1", "PARD3", "CRYAB", "PALLD")
#mitotic/dividing cells  
FeaturePlot(so, features = c("NANOG", "OCT4", "SOX2", "PAX6", "SXO1"))

so.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(so, features = top10$gene) + NoLegend()




#asingin luster identites 
new.cluster.ids <- c("","")
names(new.cluster.ids) <- levels(so)
so <- RenameIdents(so, new.cluster.ids)
DimPlot(so, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(so, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave("UMAP_Clusters_CellTypes.png", height = 7, width = 12, plot = plot, quality = 50)

