import scvi
import scipy
import sys
import anndata
import mudata
import pandas as pd
import argparse
import torch
import numpy as np

torch.set_float32_matmul_precision("medium")

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-o", "--output", help = "Name of output h5ad file.")
parser.add_argument("-m", "--model", help = "Path to trained, saved model dir.")
parser.add_argument("--proportion_cells", default=1.0, type=float, help = "Proportion of cells to use for calculating normalized expression values, for memory concerns.")
parser.add_argument("--protein_list", help = "List of proteins to pull normalized expression for.")
args = parser.parse_args()

if args.protein_list is not None:
    protein_list=args.protein_list.split(',')
else:
    protein_list=args.protein_list

adata = anndata.read_h5ad(args.adata)
model = scvi.model.SCVI.load(args.model, adata = adata)

cells_to_sample = np.random.choice(np.arange(0, len(adata.obs_names)), size=int(np.ceil(np.multiply(len(adata.obs_names), args.proportion_cells))))
print(f"Sampling from {len(cells_to_sample)} / {len(adata.obs_names)}")

#adata.layers[f"totalvi_normalized_expression_n_latent_20"] = model.get_normalized_expression(adata = adata, n_samples=10, library_size="latent", batch_size=int(np.ceil(np.true_divide(len(adata.obs_names), 50))), protein_list=protein_list, indices = cells_to_sample) 
#print(type(model.get_normalized_expression(adata = adata, n_samples=10, library_size="latent", batch_size=int(np.ceil(np.true_divide(len(adata.obs_names), 50))), protein_list=protein_list, indices = cells_to_sample)))

#'''
chunks = np.array_split(np.arange(0, len(adata.obs_names)), 5)
gex_normexp = scipy.sparse.csc_array((1,len(adata.var_names)), dtype=np.float32)
for chunk in chunks:
    print(chunk)
    normexp_this_chunk = model.get_normalized_expression(adata = adata, n_samples=10, library_size="latent", indices=chunk)
    gex_normexp = scipy.sparse.vstack( (gex_normexp, normexp_this_chunk.astype(pd.SparseDtype("float32",0)).sparse.to_coo().tocsc()) )

adata.layers["scvi_normalized_expression_n_latent_20"] = gex_normexp[1:,:]
adata.write_h5ad(args.output)
