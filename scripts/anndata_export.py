import argparse

import anndata

import pmbi.anndata.io

parser = argparse.ArgumentParser()
parser.add_argument(
    "-a", "--h5ad", help="Path to anndata object stored in h5ad format."
)
parser.add_argument(
    "-o",
    "--output_dir",
    help="Prefix for writing output files.",
)
args = parser.parse_args()

adata = anndata.read_h5ad(args.h5ad)
print(adata.layers)
pmbi.anndata.io.export(adata=adata, output_dir=args.output_dir)

