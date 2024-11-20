# %%
import anndata 
import pandas as pd
from pathlib import Path
import scanpy as sc

PROJ=Path("/home/ubuntu/projmnt/tcsl")

# %%
norm = pd.read_csv(PROJ.joinpath("anndata/parts/286585325/cell_by_gene_norm.csv"), sep = ",")
raw = pd.read_csv(PROJ.joinpath("anndata/parts/286585325/cell_by_gene_raw.csv"), sep = ",")
obs = pd.read_csv(PROJ.joinpath("anndata/parts/286585325/obs.tsv"), sep = "\t")
var = pd.read_csv(PROJ.joinpath("anndata/parts/286585325/var.tsv"), sep = "\t")

# %%
adata = anndata.AnnData(X = norm)
adata.raw = anndata.AnnData(X=raw)
adata.obs = obs
adata.var = var

# %%
adata

