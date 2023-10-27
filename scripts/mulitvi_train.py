import scvi
import sys
import anndata
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-o", "--output", help = "Name of output h5ad file.")
parser.add_argument("-n", "--n_latent_values", help = "Comma-separated list of latent values for trains models for.")
parser.add_argument("-g", "--n_genes", type = int, help = "Number of genes to use for integration.")
args = parser.parse_args()

adata = anndata.read_h5ad(args.adata)
scvi.model.MULTIVI.setup_anndata(adata, layer="raw_counts", protein_expression_obsm_key="adt", batch_key="orig_ident")

n_latent_values = [int(x) for x in args.n_latent_values.split(',')]

models = {}
reps = {}
for nlv in n_latent_values:
    models[nlv] = scvi.model.MULTIVI(adata, n_latent = nlv, n_genes = args.n_genes)
    models[nlv].view_anndata_setup()
    if args.n_genes is not None:
        models[nlv].train(n_genes = args.n_genes)
    else:
        models[nlv].train()

    reps[nlv] = models[nlv].get_latent_representation()
    adata.layers[f"multivi_normalized_expression_n_latent_{nlv}"] = models[nlv].get_normalized_expression(n_samples=10, library_size="latent")

    repkey = f"X_multivi_n_latent_{nlv}"
    adata.obsm[repkey] = reps[nlv]

    models[nlv].save(f"{repkey}_model", overwrite=True)

    #pd.DataFrame(adata.obsm[repkey]).to_csv(f"simoni_{repkey}.csv")

adata.write_h5ad(args.output)

