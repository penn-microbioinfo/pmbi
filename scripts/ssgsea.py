import gseapy
import pmbi.wrappers.scanpy as scp
from pathlib import Path
import pandas as pd
import glob
import gc

gene_sets = pd.read_csv("/home/ubuntu/mnt/ankita/signatures/genesets_scRNA-Seq_Tcell-types.csv")
gene_sets_dict = {}

for col in gene_sets.columns:
    gene_sets_dict[col] = [g for g in gene_sets[col] if not pd.isna(g)]
gene_sets_dict
adatas = scp.read_h5ad_multi([Path(x) for x in glob.glob("/home/ubuntu/mnt/ankita/combined_adatas_by_donor/*leiden*.h5ad")])

for key,adata in adatas.items():
    gc.collect()
    adata.X = adata.layers["normExp_scvi_n_latent_20"]
    res = gseapy.ssgsea(data = adata.to_df().transpose(), gene_sets = gene_sets_dict, outdir = f"/home/ubuntu/mnt/ankita/combined_adatas_by_donor/ssgsea/chenetal2021_scRNA/{key}", threads=12, seed = 55, verbose = True, sample_norm_method = "rank")

gene_sets_dict = {}
gene_sets = pd.read_csv("/home/ubuntu/mnt/ankita/signatures/genesets_bulkRNA-Seq_Tcell-types.csv")
for col in gene_sets.columns:
    gene_sets_dict[col] = [g for g in gene_sets[col] if not pd.isna(g)]
gene_sets_dict
adatas = scp.read_h5ad_multi([Path(x) for x in glob.glob("/home/ubuntu/mnt/ankita/combined_adatas_by_donor/*leiden*.h5ad")])

for key,adata in adatas.items():
    gc.collect()
    adata.X = adata.layers["normExp_scvi_n_latent_20"]
    res = gseapy.ssgsea(data = adata.to_df().transpose(), gene_sets = gene_sets_dict, outdir = f"/home/ubuntu/mnt/ankita/combined_adatas_by_donor/ssgsea/chenetal2021_bulkRNA/{key}", threads=16, seed = 55, verbose = True, sample_norm_method = "rank")
