conda activate single_cell
cd /cephfs2/hannas/tetraploids/analysis/preprocessing
python3.10

import scanpy as sc
import pandas as pd
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
plt.close('all')
matplotlib.use('Agg')

#read in the data 
path2data = "/cephfs2/hannas/tetraploids/expdata/data/combined/all-sample/DGE_unfiltered"
counts = scipy.io.mmread(f"{path2data}/count_matrix.mtx")
metadata = pd.read_csv(f"{path2data}/cell_metadata.csv")
genes = pd.read_csv(f"{path2data}/all_genes.csv", header=None)

#column names are resent in 0th row so i want to chnage that to column names and remove it from rows 
genes = genes.drop(0)
genes.columns = ['gene_id','gene_name','genome']
var_names = genes.iloc[:, 1].tolist()

#makign adata
adata = sc.AnnData(X=counts, obs=metadata, var=genes)
adata.obs_names = metadata['bc_wells'] 
adata.var_names = var_names
adata.var_names_make_unique()

from scipy.sparse import csr_matrix
adata.X = csr_matrix(adata.X)

#adata.X.shape
#(3273556, 63140)

#QC
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

adata.write("/cephfs2/hannas/tetraploids/analysis/adata_beforeQC.h5ad")

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
plt.savefig('violin_QC_plots.png', dpi=300)
plt.close()

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig('cell_complexity_coloured_by_prc_mt.png', dpi=300)
plt.close()

#removing cells with low gene and coutns numbers and removign gneses expressed in evry few cells 
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=10)
#adata.shape (677935, 38958)
sc.pp.filter_cells(adata, min_counts=200)
#adata.shape (298557, 38958)

#doublet detection 
sc.pp.scrublet(adata, batch_key="sample")


#Normalisation 

#Diamntionality Reduction 

