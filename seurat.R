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
metadata <- read.csv("/mnt/tetraploids/analysis/preprocessing/metadata.csv")
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
nFeature_bySampleName <- VlnPlot(so, features = "nFeature_RNA", pt.size = 0, group.by = "sample", flip = TRUE, log = TRUE) 
nFeature_bySampleName <- nFeature_bySampleName + theme(
  axis.text.x = element_text(size = 8),
  axis.title.x = element_text(size = 12), 
  plot.title = element_text(size = 14)
) + labs(
  x = "Sample Name",  
  y = "log[nFeature]"    
)
ggsave("nFeature_bySampleName_violin.png", plot = nFeature_bySampleName, width = 10, height = 8, dpi = 300)

