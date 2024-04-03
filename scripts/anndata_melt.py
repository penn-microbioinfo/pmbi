import anndata
import argparse
import pandas as pd
import scipy.io

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--h5ad", help = "Path to anndata object stored in h5ad format.")
parser.add_argument("-m", "--metadata_out", help = "Filename (csv) to store the the metadata to for import into Seurat.")
parser.add_argument("-s", "--splits_outdir", help = "Directory to output split matrices and metadata.")
parser.add_argument("-i", "--sampleid_col", help = "Column name in adata.obs that ties barcodes back to original sample id.")
args = parser.parse_args()

adata = anndata.read_h5ad(args.h5ad)
adata = anndata.read_h5ad("arutyunyan_et_al/Organoid_PTO_cellxgene.h5ad")
adata.obs

adata.obs.to_csv(args.metadata_out, sep="\t")
#adata._raw.X[1:10, 1:10].toarray()
adata_raw = adata._raw.to_adata()
del adata

adata_raw.obs["orig_ident"] = adata_raw.obs[args.sampleid_col]
adata_raw.obs["old_index"] = adata_raw.obs.index
adata_raw.var["feature"] = adata_raw.var.index
adata_raw.obs = adata_raw.obs.reset_index()

for oi in adata_raw.obs["orig_ident"].unique():
    which = adata_raw.obs.loc[adata_raw.obs["orig_ident"] == oi].index
    which_adata = adata_raw[which, :]
    scipy.io.mmwrite(f"{args.splits_outdir}/{oi}.mm", which_adata.X.transpose())
    which_adata.obs.to_csv(f"{args.splits_outdir}/{oi}.obs", sep = "\t", header = False, index = False)
    which_adata.var.to_csv(f"{args.splits_outdir}/{oi}.var", sep="\t", header = False, index = False)
