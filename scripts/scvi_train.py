import scvi
import sys
import anndata
import pandas as pd
import argparse
import pathlib
import mudata
import torch

torch.set_float32_matmul_precision("highest")

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-o", "--output", help = "Name of output h5ad file.")
parser.add_argument("-n", "--n_latent_values", help = "Comma-separated list of latent values for trains models for.")
parser.add_argument("-b", "--batch_key", default="orig_ident", help = "Key in adata.obs to use as batch key.")
parser.add_argument("-l", "--layer", default="counts", help = "Layer in adata.layers that contains the raw counts.")
args = parser.parse_args()

adata_path = pathlib.Path(args.adata)
if adata_path.suffix == ".h5mu":
    adata = mudata.read_h5mu(args.adata)["gex"]
elif adata_path.suffix == ".h5ad":
    adata = anndata.read_h5ad(args.adata)
else:
    raise IOError("Unrecognized file extension: {adata_path.suffix}")

scvi.model.SCVI.setup_anndata(adata, layer=args.layer, batch_key=args.batch_key)

n_latent_values = [int(x) for x in args.n_latent_values.split(',')]

models = {}
reps = {}
for nlv in n_latent_values:
    models[nlv] = scvi.model.SCVI(adata, n_latent = nlv)
    models[nlv].train()

    reps[nlv] = models[nlv].get_latent_representation()
    #adata.layers[f"scvi_normalized_expression_n_latent_{nlv}"] = models[nlv].get_normalized_expression(n_samples=10, library_size=5e4)

    repkey = f"X_scvi_n_latent_{nlv}"
    adata.obsm[repkey] = reps[nlv]

    models[nlv].save(f"{repkey}_model", overwrite=True)

    #pd.DataFrame(adata.obsm[repkey]).to_csv(f"simoni_{repkey}.csv")

adata.write_h5ad(args.output)
