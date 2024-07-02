import anndata
import argparse
import pandas as pd
import scipy.io

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--h5ad", help = "Path to anndata object stored in h5ad format.")
parser.add_argument("-o", "--output_prefix", help = "Filename (csv) to store the the metadata to for import into Seurat.")
args = parser.parse_args()

adata = anndata.read_h5ad(args.h5ad)


