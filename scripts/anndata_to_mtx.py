import argparse

import anndata

import pmbi.anndata.io

parser = argparse.ArgumentParser()
parser.add_argument(
    "-a", "--h5ad", help="Path to anndata object stored in h5ad format."
)
parser.add_argument(
    "-o",
    "--output_prefix",
    help="Prefix for writing output files.",
)
parser.add_argument(
    "-l",
    "--layer",
    required=False,
    help="Layer in AnnData object to pull counts from. Default: main layer",
)
args = parser.parse_args()

adata = anndata.read_h5ad(args.h5ad)
pmbi.anndata.io.to_mtx(adata=adata, output_prefix=args.output_prefix, layer=args.layer)
