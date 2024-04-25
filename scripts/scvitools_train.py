import argparse
import sys
import torch
import numpy as np
import pmbi.wrappers.scvi as pmbint

torch.set_float32_matmul_precision("medium")

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--data", help = "Name of input h5ad/h5mu file.")
parser.add_argument("-b", "--batch", default = None, help = "Batch key in anndata.obs")
parser.add_argument("-l", "--layer", default = None, help = "Layer to access for raw counts.")
parser.add_argument("-o", "--output_prefix", help = "Prefix for output files.")
parser.add_argument("-n", "--n_latent_values", help = "Comma-separated list of latent values for trains models for.")
args = parser.parse_args()

m = pmbint.ScviModeler(batch_key=args.batch)
m.read(args.data)
m.setup_data(layer = args.layer)
for nlv in args.n_latent_values.split(','):
    nlv=int(nlv)
    print(f"Training SCVI with n_latent = {nlv}")
    m.train(nlv)
    lr = m.get_latent_repr(nlv, as_obsm = False)
    #ne = m.get_normalized_expr(nlv, as_layer= False)
    np.save(f"{m.sample_name}_X_scvi_n_latent_{nlv}", lr)
    #np.save(f"{m.sample_name}_normExp_scvi_n_latent_{nlv}", ne)
    m.save_model(nlv) 
