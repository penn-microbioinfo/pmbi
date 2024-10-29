import scanpy as sc
import numpy as np
import squidpy as sq
import pandas as pd
from pathlib import Path
import pmbi.plotting as pmbip
import pmbi.wrappers.scanpy as scp
import importlib
importlib.reload(scp)
importlib.reload(pmbip)

# %%
OUTS={
        "Slide_1": Path("/home/ubuntu/mnt/bhoomi-bcl-dl/spaceranger/Slide_1_HE/outs/binned_outputs"),
        "Slide_2": Path("/home/ubuntu/mnt/bhoomi-bcl-dl/spaceranger/Slide_2_HE/outs/binned_outputs"),
        }

FIGS=Path("/srv/http/bhoomi")
SQUARE_SIZE="square_016um"

# %%

# %%
adatas = {}
for sample,outpath in OUTS.items():
    tissue_pos_path = outpath.joinpath(f"{SQUARE_SIZE}/spatial/tissue_positions.csv")
    if not tissue_pos_path.exists():
        tp = pd.read_parquet(outpath.joinpath(f"{SQUARE_SIZE}/spatial/tissue_positions.parquet"))
        tp.to_csv(tissue_pos_path, sep=",", index = False)
    adata = sq.read.visium(path = outpath.joinpath(SQUARE_SIZE), counts_file="filtered_feature_bc_matrix.h5", library_id = "Slide_1")
    adata.obs["n_umi"] = adata.X.sum(axis=1)[:, 0]
    lil = adata.X.tolil()
    r,c = lil.nonzero()
    lil[r,c]=1
    adata.obs["n_genes"] = lil.sum(axis=1)[:,0]
    adata_filt = adata.copy()
    sc.pp.filter_cells(adata_filt, min_counts = 250)
    adatas[sample] = {
            "raw": adata.copy(),
            "filtered": adata_filt.copy()
            }

# %%
adatas["Slide_2"]["filtered"].obs.n_umi.mean()*16
adatas["Slide_1"]["filtered"].obs.n_umi.mean()*16

# %%
importlib.reload(pmbip)
panel = pmbip.Paneler(2,3,figsize=(12,8), format="png")
max_umis = pd.concat([v["filtered"].obs.n_umi for k,v in adatas.items()]).max()
max_genes = pd.concat([v["filtered"].obs.n_genes for k,v in adatas.items()]).max()
for sample,ads in adatas.items():
    t = pmbip.Theme().xlab("Number of UMIs").ylab("frequency")
    panel.next_ax().hist(ads["filtered"].obs["n_umi"], bins = 100, range=(0,max_umis*1.01))
    t.apply_to(panel.current_ax)
    t = pmbip.Theme().xlab("Number of genes").ylab("frequency").title(sample)
    panel.next_ax().hist(ads["filtered"].obs["n_genes"], bins = 100, range=(0,max_genes*1.01))
    t.apply_to(panel.current_ax)
    t = pmbip.Theme().xlab("Number of UMIs").ylab("Number of genes")
    panel.next_ax().scatter(x=ads["filtered"].obs["n_umi"], y=ads["filtered"].obs["n_genes"], s=0.2)
    t.apply_to(panel.current_ax)

panel.fig.savefig(FIGS.joinpath("umis_genes_dist.png"))
# %%

