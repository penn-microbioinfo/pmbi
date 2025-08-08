import argparse
import sys
import torch
import numpy as np
import pmbi.wrappers.scvi as pmbint
from pathlib import Path
import scvi

torch.set_float32_matmul_precision("medium")

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--data", help = "Name of input h5ad/h5mu file.")
parser.add_argument("-b", "--batch", default = None, help = "Batch key in anndata.obs")
parser.add_argument("-l", "--layer", default = None, help = "Layer to access for raw counts.")
parser.add_argument("-o", "--output_prefix", help = "Prefix for output files.")
parser.add_argument("-n", "--n_latent_values", help = "Comma-separated list of latent values for trains models for.")
parser.add_argument("--dl_num_workers", type=int, help = "setting for scvi.settings.dl_num_workers")
args = parser.parse_args()

if args.dl_num_workers is not None:
    scvi.settings.dl_num_workers = args.dl_num_workers

m = pmbint.ScviModeler(datafile=Path(args.data), batch_key=args.batch)
m.setup_data(layer = args.layer)
for nlv in args.n_latent_values.split(','):
    nlv=int(nlv)
    print(f"Training SCVI with n_latent = {nlv}")
    m.train(n_latent=nlv)
    lr = m.get_latent_repr(n_latent=nlv, as_obsm=False)
    #ne = m.get_normalized_expr(nlv, as_layer= False)
    np.save(f"{m.sample_name}_X_scvi_n_latent_{nlv}", lr)
    #np.save(f"{m.sample_name}_normExp_scvi_n_latent_{nlv}", ne)
    m.write() 
