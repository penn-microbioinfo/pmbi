# %%
import numpy as np
import importlib
import subprocess
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

import pmbi.cellranger.cellranger_command as crc
import pmbi.config as pmbiconf
import pmbi.file_handlers as pfh
import pmbi.plotting as pmbi
from pmbi.config import import_config
from pmbi.file_handlers import LocalBackend

# from pmbi.file_handlers import CellrangerHandler, CellrangerMultiConfigCsv, LocalBackend


# %%
# importlib.reload(pfh)
importlib.reload(crc)
# importlib.reload(pmbiconf)
config = import_config(
    "/home/amsesk/super1/t1d-coculture/config/cellranger_config.toml"
)

h = crc.CellrangerCollection(
    path=Path("/home/amsesk/super1/t1d-coculture/all_fastq_symlinks/"),
    config=config,
    pattern="^HPAP[-][0-9]+[_].+[_]R[0-9][_].+[.]fastq[.]gz$",
    backend=LocalBackend(),
)

# %% Greg says that we are not interested in the...
#       - VDJ-B samples
#       - IO samples
# ... so delete

h.table = h.table[h.table["modality"] != "VDJ-B"]
h.table = h.table[~h.table["sample"].str.contains("[_]IO[_]")]
h.table["sample"]


# %% TODO: Add verify function to figure out which samples are missing according to some other sheet
pd.set_option("display.max_rows", 150)
h.table[["sample", "modality"]]["sample"].drop_duplicates().shape


coculture_meta = pd.read_excel("/home/amsesk/super1/t1d-coculture/CoCulture_metadata.xlsx")
coculture_meta["Sample_Name"]

comp = coculture_meta[["sample", "Library_Type"]]
comp = comp[~comp["sample"].str.fullmatch("HPAP-141_CC_[0-9]+")]
comp = comp[comp["Library_Type"] != "VDJ-B"]
comp = comp[~comp["sample"].str.contains("_IO_")]
comp = comp[~comp["sample"].str.contains("_IsletSupernatant")]


(h.table[["sample"]].value_counts()==24).all()
(h.table[["sample", "modality"]].drop_duplicates()["sample"].value_counts()==3).all()


h.table.shape
s1 = h.table[["sample", "modality"]]["sample"].drop_duplicates()
s1.shape
s2 = comp["sample"].drop_duplicates()
s2.shape

[x for x in s2 if x not in s1.tolist()]
[x for x in s1 if x not in s2.tolist()]


[x for x in h.get_units() if x.sample == "HPAP-141_CC_1a"][0].table
# %%
hur = [x for x in h.get_units() if x.sample == "HPAP-141_CC_1a"][0].create_runner(
    wd=Path("/home/amsesk/super2/cellranger/"),
    localcores=8,
    localmem=32
)
hur._check_output_exists()
dir(hur)

# %%

# %%

def _run(unit, wd, **kwargs):
    unit.create_runner(wd, **kwargs).run()
# %%
_out = Parallel(n_jobs=12)(delayed(_run)(unit,
                                        wd=Path("/home/amsesk/super2/cellranger/"),
                                        localcores=8,
                                        localmem=32
                                        ) for unit in h.get_units())

###########################
# %%
###########################


# %%
hur = hu[8].create_runner(
    wd=Path("/home/amsesk/super1/t1d-coculture/cellranger/runs"),
    localcores=8,
    localmem=32,
)


# %% Run cellranger


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
sub._modality_metadata().set_index("cellranger_multi_feature_type").loc[
    sub.table["modality"].unique()
]

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
