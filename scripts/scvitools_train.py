import argparse
import torch
import pmbi.wrappers.scvi as pmbint

torch.set_float32_matmul_precision("medium")

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--data", help = "Name of input h5ad/h5mu file.")
parser.add_argument("-o", "--output_prefix", help = "Prefix for output files.")
parser.add_argument("-n", "--n_latent_values", help = "Comma-separated list of latent values for trains models for.")
args = parser.parse_args()

m = ScviModeler(args.output)
m.read(args.data)
m.setup_data()
m.train(args.n_latent_values)
m.write(f"{args.output_prefix}_latentReprs.h 
