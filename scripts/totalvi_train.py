icombined_adatas_by_donor/{pid}.h5ad")mport scvi
import sys
import anndata
import mudata
import pandas as pd
import argparse
import torch

torch.set_float32_matmul_precision("medium")

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-o", "--output", help = "Name of output h5ad file.")
parser.add_argument("-n", "--n_latent_values", help = "Comma-separated list of latent values for trains models for.")
args = parser.parse_args()

adata = mudata.read_h5mu(args.adata)
adata["adt"].X = adata["adt"].X.toarray()
print(type(adata["gex"].X), type(adata["adt"].X))
print(adata.shape, adata["gex"].shape, adata["adt"].shape)
scvi.model.TOTALVI.setup_mudata(adata, modalities={"rna_layer": "gex", "protein_layer": "adt"}, rna_layer = "raw_counts", batch_key="orig_ident")

n_latent_values = [int(x) for x in args.n_latent_values.split(',')]

models = {}
reps = {}
for nlv in n_latent_values:
    models[nlv] = scvi.model.TOTALVI(adata, n_latent = nlv)
    models[nlv].view_anndata_setup()
    print("training")
    models[nlv].train()

    print("done training")
    reps[nlv] = models[nlv].get_latent_representation()
    #astype(pd.SparseDtype("float64",0)).sparse.to_coo().tocsc()adata.layers[f"totalvi_normalized_expression_n_latent_{nlv}"] = models[nlv].get_normalized_expression(n_samples=10, library_size="latent")

    repkey = f"X_totalvi_n_latent_{nlv}"
    adata.obsm[repkey] = reps[nlv]

    models[nlv].save(f"{repkey}_model", overwrite=True)

    #pd.DataFrame(adata.obsm[repkey]).to_csv(f"simoni_{repkey}.csv")

adata.write_h5mu(args.output)

