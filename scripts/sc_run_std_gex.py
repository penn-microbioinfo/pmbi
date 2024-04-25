import scanpy as sc
import pmbi.wrappers.scanpy as scp
import argparse

sc._settings.ScanpyConfig.n_jobs = 7

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--matrix", required = True, action = "store", help = "Path to 10X filtered_feature_bc_matrix.h5")
parser.add_argument("-s", "--suffix", required = True, action = "store", help = "Samlple suffix to name output files and plots.")
#parser.add_argument("-", "--", required = True, action = "store", help = "")
args = parser.parse_args()

adata = scp.read_matrix(args.matrix, gex_only = True)
adata.layers["counts"] = adata.X.copy()
adata.var_names_make_unique()
adata_gex_qc = scp.std_qc_gex(adata, sample_suffix = args.suffix)
adata_gex_qc.write(f"{args.suffix}_std_gex_qc.h5ad")
gex = scp.std_gex(adata_gex_qc, sample_suffix = args.suffix)
gex.write(f"{args.suffix}_std_gex.h5ad")
