# %%
from pathlib import Path

import anndata
import mudata
import muon as mu
import muon.atac as atac
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse
from scipy.io import mmread
import pmbi.cellranger.cellranger_command as crc
import pmbi.config as pmbiconf

import pmbi.anndata.io as aio
import pmbi.bio.gtf as pgtf
import importlib

# %%
paths=[
    Path("/home/amsesk/super2/jayme_shiv/gex_fastq/H3MY3DSXF/symlinks"),
    Path("/home/amsesk/super2/jayme_shiv/atac_fastq_filtered/symlinks"),
    Path("/home/amsesk/super2/jayme_shiv/adtHto_fastq/symlinks")
]
ccc = crc.CellrangerCollection(
    path=paths,
    config=pmbiconf.import_config("/home/amsesk/super2/jayme_shiv/pmbi_cellranger_config.toml")
)

# %%
wd = Path("/home/amsesk/super2/jayme_shiv")

scc_proc_wd = wd.joinpath("scc-proc_outs")
samples = ccc.table["sample"].drop_duplicates().tolist()

adt_catalog = pd.read_csv("/home/amsesk/super1/ref/adt_catalogs/2023_NHP_catalog_TotalSeqA_antibodies_barcodes_v2.csv")
adt_catalog = adt_catalog.set_index("DNA_ID")

# Prepare list of mitochondrial genes to annotate RNA adatas as such
with open("/home/amsesk/super1/ref/genomes/Mmul_10__SHIV.C.CH505.v2/genes.gtf", 'r') as h:
    gtf = pgtf.FeatureFile.from_gtf(h)

gene_gtf = gtf[gtf["feature"]=="gene"]
gene_gtf["attributes"] = gene_gtf["attributes"].apply(lambda a: pgtf.FeatureFile._parse_attributes_string(a))
mitro_chrom="NC_005943.1"
mito_gene_gtf_attr = gene_gtf[gene_gtf["seqname"]=="NC_005943.1"]["attributes"].tolist()
mito_gene_ids = np.array([ge["gene_id"] for ge in mito_gene_gtf_attr])

# %% Generate ADT AnnData objects from scc_proc outputs 
importlib.reload(aio)
adt_adatas = {}
for s in samples:
    rundir_base = f"{s}_adt_output_count"
    rundir = scc_proc_wd.joinpath(rundir_base)
    mtx_path = rundir.joinpath("output.mtx")
    obs_names_path = rundir.joinpath("output.barcodes.txt")
    var_names_path = rundir.joinpath("output.genes.txt")
    adt_adata = aio.read_mtx(
            mtx_path=mtx_path,
            obs_names_path=obs_names_path,
            var_names_path=var_names_path
    )
    adt_adata.var_names.name="DNA_ID"
    adt_adata.var = adt_catalog.loc[adt_adata.var_names,:]
    # adt_adata.var["DNA_ID"] = adt_adata.var_names
    # adt_adata.var = adt_adata.var.set_index("Gene Name", drop=False)
    # adt_adata.var_names = adt_adata.var["Gene Name"]
    adt_adata.obs_names = pd.Series(adt_adata.obs_names).apply(lambda bc: f"{bc}-1")
    adt_adata.var = adt_adata.var.rename(columns={"Gene Name": "gene_ids"})
    adt_adata.var["feature_types"] = "Antibody"
    adt_adata.var["IsotypeControl"] = adt_adata.var["IsotypeControl"].apply(lambda v: False if np.isnan(v) else True)
    # adt_adata.var["genome"] = "jayme_custom"
    # adt_adata.var["interval"] = np.nan
    adt_adatas[s] = adt_adata

# %% Read in RAW cellranger-arc h5 files
arc_mdatas = {}
for s in samples:
    cr_arc_runs = wd.joinpath("cr_arc_runs")
    arc_raw = mu.read_10x_h5(cr_arc_runs.joinpath(f"{s}/outs/raw_feature_bc_matrix.h5"))
    arc_mdatas[s] = arc_raw

# %% Read in FILTERED cellranger-arc h5 files
arc_mdatas_filt = {}
for s in samples:
    cr_arc_runs = wd.joinpath("cr_arc_runs")
    arc_filt = mu.read_10x_h5(cr_arc_runs.joinpath(f"{s}/outs/filtered_feature_bc_matrix.h5"))
    arc_mdatas_filt[s] = arc_filt

# %% Combine ADT and cellranger-arc adatas into combined mdata
combined_mdatas={}
for s in samples:
    adt_adata = adt_adatas[s]
    arc_raw = arc_mdatas[s]
    shared_adt_raw_arc = adt_adata.obs_names[np.where([x in arc_raw.obs_names for x in adt_adata.obs_names])]
    arc_raw_sub = arc_raw[shared_adt_raw_arc,:]
    adt_adata_sub = adt_adata[shared_adt_raw_arc,:]
    # Add MT annotations to gene_ids for rna and atac adatas
    for modal in ["rna"]:
        arc_raw_sub[modal].var.loc[mito_gene_ids,"gene_ids"] = arc_raw_sub[modal].var.loc[mito_gene_ids,"gene_ids"].apply(lambda gid: f"mt-{gid}")
        arc_raw_sub[modal].var_names = arc_raw_sub[modal].var["gene_ids"]
        arc_raw_sub[modal].var["mito"] = arc_raw_sub[modal].var["gene_ids"].str.startswith("mt-")
    combined = mudata.MuData({
        "rna": arc_raw_sub["rna"],
        "atac": arc_raw_sub["atac"],
        "adt": adt_adata_sub
    })
    combined.obs["cr_is_cell"] = combined.obs_names.isin(arc_mdatas_filt[s].obs_names)
    combined.write(wd.joinpath(f"h5mu/{s}_raw.h5mu"))
    combined_mdatas[s]=combined

# %%
