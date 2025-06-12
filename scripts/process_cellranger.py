# t%%
import subprocess
import importlib
import re
import os
from pathlib import Path
from joblib import Parallel, delayed

import pandas as pd

import pmbi.cellranger.cellranger_command as cc
import pmbi.file_handlers as pfh
import pmbi.plotting as pmbi
from pmbi.config import import_config
from pmbi.file_handlers import CellrangerHandler, CellrangerMultiConfigCsv, LocalBackend

import pmbi.cellranger.cellranger_command as crc

# %%
importlib.reload(pfh)
importlib.reload(cc)
config = import_config(
    "/home/amsesk/super1/t1d-coculture/config/cellranger_config.toml"
)
config.references

h = crc.CellrangerCollection(
    path=Path("/home/amsesk/super1/t1d-coculture/all_fastq_symlinks/"),
    config=config,
    pattern="^HPAP[-][0-9]+[_].+[_]R[0-9][_].+[.]fastq[.]gz$",
    backend=LocalBackend(),
)

hu = h.get_units()
hu[0].table

h.table
# %%
wd = "/home/amsesk/super1/t1d-coculture/cellranger/runs/"
os.chdir(wd)
h_spl = h.split(by="sample_rep")
cmds = [pfh.CellrangerMulti(handler=h, wd=Path("/home/amsesk/super1/t1d-coculture/cellranger/runs"), localcores=8, localmem=32).cmd() for h in h_spl]


# %%
def _run(cmd):
    p = subprocess.run(cmd, capture_output=True)
    return p

out = Parallel(n_jobs=24)(
    delayed(_run)(cmd) for cmd in cmds
)

# %%

importlib.reload(pfh)
sub = h.subset(lambda x: x.sample_rep == "HPAP-135_LO_Run1Day0_1")
csv = pfh.CellrangerMultiConfigCsv(sub)
csv.handler.table["modality"]
print(csv._libraries_section())
csv._feature_type_sections()
csv.write(path=Path("/home/amsesk/super1/t1d-coculture/test_cr.csv"))

# %%

h.table["modality"].unique()
sub._modality_metadata().set_index("cellranger_multi_feature_type").loc[sub.table["modality"].unique()]

# %%
sub._modality_metadata()


sec = pd.DataFrame(
    [{"modality": k, "reference": v} for k, v in config.references.__dict__.items()]
)
sec["section"] = (
    sec["modality"]
    .replace(sub.feature_type_converter)
    .replace(csv._section_header_converter())
)
sec[["section", "reference"]].melt(id_vars="section")


# %%
pd.set_option("display.max_rows", 200)
h.table["modality"] = pd.Categorical(
    values=h.table["modality"], categories=h.table["modality"].unique()
)
modal_counts = (
    h.table.groupby(["sample_rep", "modality"])[["modality"]]
    .value_counts()
    .reset_index()
)
modal_counts.columns
modal_counts.pivot_table(values="count", index="sample_rep", columns="modality")
modal_counts.where(modal_counts["count"] == 0).dropna()

# %%
panel = pmbi.Paneler(1, 1, (8, 10))
pmbi.heatmap(
    matrix=modal_counts.pivot_table(
        values="count", index="sample_rep", columns="modality"
    ),
    ax=panel.next_ax(),
    xlab="modality",
    ylab="Sample_replicate",
)
panel.fig.savefig("/home/amsesk/figures/coculture/modal_cover.png")
# %%

