import importlib
import logging
import re

import awkward as ak
import pandas as pd
from joblib import Parallel, delayed

from pmbi.airr.schema import get_airr_schema

# %%
importlib.reload(pmbi.airr.schema)
importlib.reload(pmbi.anndata.io)
importlib.reload(pmbi.util)
import os

from pmbi.airr.schema import get_airr_schema
from pmbi.anndata.io import from_airr
from pmbi.util import ak_print

####################################
# %% Read Airr TSVs into AnnData objects, in parallel
####################################
n_jobs = 8
file_list = [
    os.path.join("/home/ubuntu/projmnt/betts/coculture/mixcr/", ff)
    for ff in os.listdir("/home/ubuntu/projmnt/betts/coculture/mixcr/")
    if ff.endswith("clns.airr.tsv")
]
sample_pattern = "^(HPAP[0-9]+[-]rep[0-9]).*$"


# %%
def _compute(f, sample_pattern):
    sample = re.sub(sample_pattern, "\\1", os.path.basename(f))
    logging.warning(f"Working on {sample}")
    mixcr_airr = pd.read_csv(f, sep="\t")
    mixcr_adata = from_airr(mixcr_airr, schema="Rearrangement", sample_key=sample)
    return mixcr_adata


# %%
adatas = Parallel(n_jobs=n_jobs)(
    delayed(_compute)(f, sample_pattern) for f in file_list
)

# %%
concat = anndata.concat(adatas)

# %%
concat.write_h5ad("/home/ubuntu/projmnt/betts/coculture/mixcr/donors_combined.h5ad")
####################################
# %%
####################################
len(anndata.read_h5ad("/home/ubuntu/projmnt/betts/coculture/mixcr/donors_combined.h5ad").obsm["airr"].junction_aa[...,0].to_list())

mixcr_adata.obsm["airr"].junction_aa[..., 0]
ak_print(mixcr_adata.obsm["airr"][5])
mixcr_adata.obsm["airr"][["junction_aa", "junction"]][..., 0].show()


ak_print(mixcr_adata.obsm["airr"])

# %%
for row in mixcr_airr.itertuples():
    print(row)
    print(getattr(row, "duplicate_count"))
    break

# %%
