import argparse
import functools
import importlib
import random
import re
from os import PathLike
from pathlib import Path

import anndata
import numpy as np
import palettable
import pandas as pd
import scanpy as sc
import scipy.stats
import toolz
from joblib import Parallel, delayed
from joblib_progress import joblib_progress
from tqdm import tqdm
import matplotlib.pyplot as plt

import pmbi.anndata.io
import pmbi.plotting as pmbip
import pmbi.scanpy.qc.cutoffs as pqc
import pmbi.wrappers.scanpy as scp
from pmbi.cellranger.globals import \
    ChromiumNextGEMSingleCell5primeHTv2__multiplet_rates

importlib.reload(scp)

sc._settings.ScanpyConfig.n_jobs = 1
tqdm.pandas()

#################
# %% Project specific functions
#################
def experimentName_to_subcategory(en):
    """
    Function specific to the Betts Coculture project 

    Extracts the larger experiment type from the more specific experiment name
    """
    en_suff = en.split("_", 1)[1]
    s = re.search("Day[0-9]", en_suff)
    match s:
        case None:
            if en_suff.endswith("CC"):
                return "CC"
            else:
                return pd.NA
        case _:
            return s.group(0)

#################
# %% ArgumentParser
#################
# parser = argparse.ArgumentParser()
# parser.add_argument(
#     "-m",
#     "--matrix",
#     required=True,
#     action="store",
#     help="Path to 10X filtered_feature_bc_matrix.h5",
# )
# parser.add_argument(
#     "-s",
#     "--suffix",
#     required=True,
#     action="store",
#     help="Samlple suffix to name output files and plots.",
# )
# parser.add_argument(
#     "-o", "--figout", default=".", action="store", help="Directory to output figures."
# )
# # parser.add_argument("-", "--", required = True, action = "store", help = "")
# args = parser.parse_args()


# matrix_dir = Path("/home/amsesk/super2/cellranger/HPAP-135_CC_1__outs/per_sample_outs/HPAP-135_CC_1/count/sample_filtered_feature_bc_matrix.h5")
# %%
importlib.reload(scp)
figout = Path("/home/amsesk/figures/coculture/qc_kde_figures/png")
matrix_dir = Path("/home/amsesk/super2/cellranger/")
mt_prefix = "MT-"
cols = (palettable.cartocolors.qualitative.Pastel_10.mpl_colors) * 10

################################
# %% Read in CoCulture metadata and convert to be easily merge-able with adata.obs
################################
coculture_meta = pd.read_excel(
    "/home/amsesk/super1/t1d-coculture/CoCulture_metadata.xlsx"
)
meta_to_obs = (
    coculture_meta[
        [
            "Experiment_Name",
            "Sample_Name",
            "DonorID",
            "Tissue",
            "Disease_Status",
            "Sequencing_Run",
            "Sex",
            "Age",
            "Library_Type",
        ]
    ]
    .query("Library_Type == 'RNA'")
    .rename(columns={"Sample_Name": "sample_id"})
    .assign(
        Experiment_subcategory=lambda r: r["Experiment_Name"].apply(
            experimentName_to_subcategory
        )
    )
    .drop(columns=["Library_Type"])
    .set_index("sample_id")
)
meta_to_obs = meta_to_obs[~meta_to_obs["Experiment_subcategory"].isna()]

################
# %% Read in adatas and do a little prepocessing, calculate qc metrics, etc
################
fs = list(matrix_dir.iterdir())
adatas = {}
for o in fs:
    sample_id = o.name.split("__")[0]
    matrix_path = o.joinpath(
        f"per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix.h5"
    )
    adata = pmbi.anndata.io.read_matrix(matrix_path, gex_only=True)
    adata.var_names_make_unique()
    adata.obs["sample_id"] = sample_id
    adata.obs["original_barcode"] = adata.obs.index.to_list()
    adata.obs.index = pd.Index([f"{bc}__{sample_id}" for bc in adata.obs.index])
    adata.layers["counts"] = adata.X.copy()
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata.obs = pd.merge(left=adata.obs,
            right=pd.DataFrame(meta_to_obs["Experiment_subcategory"]),
            how="left",
            left_on="sample_id",
            right_index=True)
    adatas[sample_id] = adata

################
# %% QC scatter plots 
################

def _qc_scatter(sample_id, adata):
    panel = pmbip.Paneler(1,1,(3,3))
    mappable = panel.next_ax().scatter(x=adata.obs["total_counts"], 
                                       y=adata.obs["n_genes_by_counts"],
                                       c=adata.obs["pct_counts_mt"],
                                       cmap=palettable.scientific.diverging.Roma_20_r.mpl_colormap,
                                       s=1,
                                       )
    panel.current_ax.set_title(sample_id)
    panel.current_ax.set_xlabel("Total counts")
    panel.current_ax.set_ylabel("Unique features")
    panel.current_ax.set_xticklabels(labels=panel.current_ax.get_xticklabels(), fontsize=4, rotation=90)
    panel.current_ax.set_yticklabels(labels=panel.current_ax.get_yticklabels(), fontsize=4)
    panel.fig.colorbar(mappable=mappable, ax=panel.current_ax, shrink=0.5, label="pMito")
    panel.fig.savefig(f"/home/amsesk/figures/coculture/qc_scatters/{sample_id}.png")

for sample_id,adata in adatas.items():
    _qc_scatter(sample_id, adata)

################
# %% Function for computing cutoff values for the qc metrics
################
def _compute(sample_id, adata, no_lower=None, no_upper=None):
    if no_lower is None:
        no_lower = []
    if no_upper is None:
        no_upper = []
    cutoffs = {}
    cutoffs["sample_id"] = sample_id
    cutoffs["adata"] = adata.copy()
    qckeys = ["pct_counts_mt", "n_genes_by_counts", "total_counts"]
    for qckey in qckeys:
        com = scp.KdeCutoffMaker(
            X=adata.obs[qckey],
            gradient_threshold=1.0,
            dens_x_len=100000,
            endpoint_buffer=3,
        )
        d0 = com.get_curve(derivative_order=0)
        si = com.selection_interval()
        lower = d0.x[si.min()]
        upper = d0.x[si.max()]
        if qckey not in no_lower and qckey not in no_upper:
            limits = (lower, upper)
        else:
            if qckey in no_lower and qckey not in no_upper:
                limits = (-np.inf, upper)
            elif qckey not in no_lower and qckey in no_upper:
                limits = (lower, np.inf)
            else:
                limits = (-np.inf, np.inf)
        adata_filt = adata.copy()
        adata_filt = adata_filt[
            (adata_filt.obs[qckey] >= limits[0]) & (adata_filt.obs[qckey] <= limits[1])
        ]
        cutoffs[qckey] = pd.Series(
            {
                "limits": limits,
                "cells_include": adata_filt.shape[0],
                "cells_total": adata.shape[0],
            }
        )
    #     panel = pmbip.Paneler(ncol=3, nrow=1, figsize=(12, 3))
    #     panel.next_ax().plot(*d0, linestyle="-", c="r")
    #     panel.current_ax.axvspan(xmin=limits[0], xmax=limits[1], color=cols[4])
    #     panel.current_ax.set_title(qckey, loc="center")
    # panel.fig.suptitle(sample_id)
    # panel.fig.savefig(figout.joinpath(f"{sample_id}_qc_kde_diagnostic.png"))
    return cutoffs


# %% Compute cutoffs in parallel and then combine the results into a DataFrame
with joblib_progress("Computing cutoff values...", total=len(adatas)):
    cutoffs_out = Parallel(n_jobs=32)(
        delayed(_compute)(sample_id=si, adata=ad, no_lower=["pct_counts_mt"], no_upper=["n_genes_by_counts", "total_counts"],) for si,ad in adatas.items()
        # delayed(_compute)(sample_id=sample_id, adata=adata, no_lower=["pct_counts_mt"], no_upper=[]) for sample_id,adata in adatas.items()
    )

dfdict = {}
for key in ["sample_id", "adata", "pct_counts_mt", "n_genes_by_counts", "total_counts"]:
    dfdict[key] = [x[key] for x in cutoffs_out]

cutoff_df = pd.DataFrame(dfdict).sort_values("sample_id").set_index("sample_id")

# %% Write cutoffs to disk
importlib.reload(pqc)
cutoff_names = ["pct_counts_mt", "n_genes_by_counts", "total_counts"]
cutoff_sheet = pqc.cutoffs_as_sheet(cutoff_df, cutoff_names)

# cutoff_sheet.to_csv("/home/amsesk/super2/qc_cutoffs/strict.csv", index=False)
cutoff_sheet.to_csv("/home/amsesk/super2/qc_cutoffs/relaxed.csv", index=False)

# %% Read cutoffs back from disk
# cutoffs = pqc.read_cutoff_sheet(path="/home/amsesk/super2/qc_cutoffs/strict.csv")
cutoffs = pqc.read_cutoff_sheet(path="/home/amsesk/super2/qc_cutoffs/relaxed.csv")

# %%
importlib.reload(pqc)
adatas_filt = {}
cutoffs.columns
for sample_id,adata in adatas.items():
    print(sample_id)
    adata_filt = adata.copy()
    for c in cutoffs.columns:
        these_cutoffs = cutoffs.loc[sample_id, c]
        adata_filt = pqc.apply_cell_cutoff_to_adata(adata_filt, qckey=c, min_value=these_cutoffs[0], max_value=these_cutoffs[1])
        print(adata_filt.shape)
    adatas_filt[sample_id] = adata_filt

######################################
# %% Multiplet rate estimates
######################################
mult_rate_df = coculture_meta[coculture_meta["Library_Type"] == "RNA"][
    ["Sample_Name", "Cells_in_Sample_Est"]
].set_index("Sample_Name")

# Estimate multiplet rates based on 2-degree polynomial fit to 10X stated rates
mult_rates = ChromiumNextGEMSingleCell5primeHTv2__multiplet_rates
mult_rate_coef = np.polyfit(
    mult_rates["n_cells_loaded"], mult_rates["multiplet_rate"], deg=2
)
mult_rate_p = np.poly1d(mult_rate_coef)
est_mult_rates = mult_rate_df.assign(
    est_mult_rate=lambda r: mult_rate_p(r["Cells_in_Sample_Est"])
)["est_mult_rate"]
mult_rate_df["est_mult_rate"] = est_mult_rates

assert (
    mult_rate_df["est_mult_rate"] == est_mult_rates.loc[mult_rate_df["est_mult_rate"].index]
).all()


# %%
def always_indiv_preproc(sample_id, adata, est_mult_rate):
    adata = adata.copy()
    sc.pp.scrublet(adata, expected_doublet_rate=est_mult_rate)
    return (sample_id, adata)

def other_preproc(adata, 
                  #min_cells, 
                  meta2merge, 
                  batch_key="sample_id"):
    adata = adata.copy()
    # sc.pp.filter_genes(adata, min_cells=min_cells)
    assert len(meta2merge.index) == len(meta2merge.index.unique())
    adata.obs = pd.merge(
        left=adata.obs,
        right=meta2merge,
        how="left",
        left_on="sample_id",
        right_index=True,
    )
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key=batch_key)
    # adata.layers["scale_data"] = sc.pp.scale(adata, max_value=10, copy=True).X
    # sc.tl.pca(adata)
    # sc.pp.neighbors(adata, use_rep="X_pca", key_added = "X_pca_neighbors")
    # sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=1.0, key_added="X_pca_neighbors_leiden_1.0", neighbors_key="X_pca_neighbors")
    return adata

####################################
# %% Do preprocessing steps that should be done invididually before concating
####################################
with joblib_progress("Doing individual anndata preprocessing...", total=len(adatas_filt)):
    adatas_filt_preproc = Parallel(n_jobs=32)(
        # delayed(_compute)( o, mt_prefix, cols, figout, no_lower=["pct_counts_mt"], no_upper=["n_genes_by_counts", "total_counts"],) for o in fs
        delayed(always_indiv_preproc)(sample_id=sample_id, adata=adata, est_mult_rate=mult_rate_df.loc[sample_id, "est_mult_rate"]) for sample_id,adata in adatas_filt.items()
    )
    adatas_filt_preproc = {si:ad for si,ad in adatas_filt_preproc}

# %% Save individual anndatas
pmbi.anndata.io.write_h5ad_multi(
    adatas_filt_preproc,
    suffix="filt_preproc",
    outdir="/home/amsesk/super2/h5ad/indiv_preproc/relaxedQc",
)
adatas["HPAP-135_CC_2"].obs.columns
#####################
# %% Split up the AnnData's into different combined objects for integration
#####################
CC_Day4_samps = meta_to_obs[meta_to_obs["Experiment_subcategory"].isin(["CC", "Day4"])].index.to_list()
Day0_samps = meta_to_obs[meta_to_obs["Experiment_subcategory"].isin(["Day0"])].index.to_list()

CC_Day4_adatas_filt_preproc = {si:ad for si,ad in adatas_filt_preproc.items() if si in CC_Day4_samps}
Day0_adatas_filt_preproc = {si:ad for si,ad in adatas_filt_preproc.items() if si in Day0_samps}

# %% Split and process all
CC_Day4_concat = anndata.concat(CC_Day4_adatas_filt_preproc.values())
Day0_concat = anndata.concat(Day0_adatas_filt_preproc.values())

# %%
meta_to_obs.columns
CC_Day4_concat_proc = other_preproc(CC_Day4_concat, meta_to_obs.drop(columns=["Experiment_subcategory"]), batch_key="sample_id")
Day0_concat_proc = other_preproc(Day0_concat, meta_to_obs.drop(columns=["Experiment_subcategory"]), batch_key="sample_id")
CC_Day4_concat_proc.obs.columns

# %%
CC_Day4_concat_proc.write_h5ad("/home/amsesk/super2/h5ad/CC_Day4_relaxedQc.h5ad")
Day0_concat_proc.write_h5ad("/home/amsesk/super2/h5ad/Day0_relaxedQc.h5ad")

# %% Check occupancy of genes on gene list
rna_gene_list = pd.read_excel("/home/amsesk/super2/gene_lists/umap_annotation_genes.xlsx", sheet_name="Gene")
rna_gene_list["Gene"].isin(CC_Day4_concat_proc.var_names).all()
rna_gene_list["Gene"].isin(Day0_concat_proc.var_names).all()

rna_gene_list["Gene"].isin(sc.pp.filter_genes(CC_Day4_concat_proc, min_cells=3, inplace=False).var_names).all()
np.where(rna_gene_list["Gene"].isin(CC_Day4_concat_proc[:,sc.pp.filter_genes(CC_Day4_concat_proc, min_cells=3, inplace=False)[0]].var_names))[0].shape
np.where(rna_gene_list["Gene"].isin(Day0_concat_proc[:,sc.pp.filter_genes(Day0_concat_proc, min_cells=3, inplace=False)[0]].var_names))[0].shape

# %% Combined adatas with gene filter
# CC_Day4_concat_proc_geneFilt = CC_Day4_concat_proc[:,sc.pp.filter_genes(CC_Day4_concat_proc, min_cells=3, inplace=False)[0]]
# Day0_concat_proc_geneFilt = Day0_concat_proc[:,sc.pp.filter_genes(Day0_concat_proc, min_cells=3, inplace=False)[0]]
#
# CC_Day4_concat_proc_geneFilt.write_h5ad("/home/amsesk/super2/h5ad/CC_Day4_relaxedQc_geneFilt.h5ad")
# Day0_concat_proc_geneFilt.write_h5ad("/home/amsesk/super2/h5ad/Day0_relaxedQc_geneFilt.h5ad")



######################################
# %% Concatenate filtered AnnDatas and output for integration
######################################
for adata in adata_preproc:
    assert not any([vn.startswith("Hu.") for vn in adata.var_names])


pmbi.anndata.io.write_h5ad_multi(
    {x.obs["sample_id"][0]: x for x in adata_preproc},
    suffix="preproc",
    outdir="/home/amsesk/super2/h5ad/indiv_preproc",
)
combined = anndata.concat(adata_preproc, axis=0, join="outer")
combined.obs["Experiment_subcategory"].unique()

combined
adata_cc_day4 = combined[
    combined.obs["Experiment_subcategory"].isin(["CC", "Day4"])
].copy()
sc.pp.filter_genes(adata_cc_day4, min_cells=3)
adata_cc_day4.shape
adata_cc_day4.write_h5ad("/home/amsesk/super2/h5ad/combined_CC_Day4_raw_counts.h5ad")

adata_day0 = combined[combined.obs["Experiment_subcategory"].isin(["Day0"])].copy()
sc.pp.filter_genes(adata_day0, min_cells=3)
adata_day0.shape
adata_day0.write_h5ad("/home/amsesk/super2/h5ad/combined_Day0_raw_counts.h5ad")

###################
# %% ^^^ To Integration ^^^
###################

# %% Combined multiplet rate scatter
adata_cc_day4 = sc.read_h5ad(
    "/home/amsesk/super2/h5ad/combined_CC_Day4_raw_counts.h5ad"
)
adata_day0 = sc.read_h5ad("/home/amsesk/super2/h5ad/combined_Day0_raw_counts.h5ad")


# %%
def observed_multiplet_rate(grp):
    return len(np.where(grp["predicted_doublet"])[0]) / grp.shape[0]


cc_day4_mr = pd.DataFrame(
    {
        "multiplet_rate": adata_cc_day4.obs.groupby("sample_id").apply(
            observed_multiplet_rate
        )
    }
)
cc_day4_mr["object"] = "CC_Day4"
day0_mr = pd.DataFrame(
    {
        "multiplet_rate": adata_day0.obs.groupby("sample_id").apply(
            observed_multiplet_rate
        )
    }
)
day0_mr["object"] = "Day0"
combmr = pd.concat([cc_day4_mr, day0_mr], axis=0)
combmr = pd.merge(
    left=combmr,
    right=coculture_meta[coculture_meta["Library_Type"] == "RNA"].set_index(
        "Sample_Name"
    )["Cells_in_Sample_Est"],
    how="left",
    left_index=True,
    right_index=True,
)
cc_day4_mr = combmr[combmr["object"] == "CC_Day4"]
day0_mr = combmr[combmr["object"] == "Day0"]
cc_day4_mr
# %%
panel = pmbip.Paneler(nrow=1, ncol=1, figsize=(4, 4))
panel.next_ax().scatter(
    cc_day4_mr["Cells_in_Sample_Est"],
    cc_day4_mr["multiplet_rate"],
    marker=".",
    s=5,
    c="blue",
)
panel.current_ax.scatter(
    day0_mr["Cells_in_Sample_Est"], day0_mr["multiplet_rate"], marker=".", s=5, c="red"
)
panel.fig.savefig("/home/amsesk/figures/coculture/scrub_multiplet_rate.pdf")

# %% Just making a subset for working with on macbook air
test_adata = sc.read_h5ad("/home/amsesk/super2/h5ad/combined_Day0_raw_counts.h5ad")
coculture_meta[coculture_meta["Cells_in_Sample_Est"] == 0]

# %%

sc.pp.normalize_total(test_adata)
sc.pp.log1p(test_adata)
sc.pp.highly_variable_genes(test_adata, n_top_genes=500, batch_key="sample_id")


# %%
def random_subset_adata(adata, n_cells, n_genes, from_most_variable=True):
    cells_total, genes_total = adata.shape
    cell_samp = np.random.randint(0, cells_total - 1, size=(n_cells,))
    if from_most_variable:
        var_genes_idx = np.where(test_adata.var["highly_variable"])[0]
        if n_genes > len(var_genes_idx):
            n_genes = len(var_genes_idx)
            gene_samp = var_genes_idx
        else:
            gene_samp = var_genes_idx[np.random.randint(0, n_genes, size=(n_genes,))]
    else:
        gene_samp = np.random.randint(0, genes_total - 1, size=(n_genes,))
    return adata[cell_samp, gene_samp]


random_subset_adata(test_adata, 50000, 500).write_h5ad(
    "/home/amsesk/super2/h5ad/combined_Day0_raw_counts_50000x500_SUB.h5ad"
)


# %%
# cutoff_df = cutoff_df.apply(lambda r: other_preproc(row=r, gene_filt_min_cells=3), axis=1)
# cutoff_df.iloc[0,5].shape


# %%
cutoff_df.apply(lambda r: r["adata_filt"].shape[0] / r["adata"].shape[0], axis=1).mean()
cutoff_df.loc["HPAP-135_CC_1"]["adata"].obs["n_genes_by_counts"].max()
cutoff_df.loc["HPAP-135_CC_1"]["adata"].obs[
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
].sort_values("total_counts", ascending=False).iloc[0:250, :]
cutoff_df.loc["HPAP-135_CC_1"]["adata_filt"].obs["n_genes_by_counts"].max()

# %%
est_mult_rate = ncells_mult_rate_est.loc["HPAP-135_CC_1"]["est_multiplet_rate"]
assert isinstance(est_mult_rate, np.floating)
adata = cutoff_df.loc["HPAP-135_CC_1"]["adata_filt"].copy()
sc.pp.filter_genes(adata, min_cells=3)
scrub_sim = sc.pp.scrublet_simulate_doublets(adata, layer="counts")
sc.pp.scrublet(adata, adata_sim=scrub_sim)
adata.obs[
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "predicted_doublet"]
].sort_values("total_counts", ascending=False).iloc[0:250, :]
adata.obs[adata.obs["predicted_doublet"]]

sc.pp.scrublet(adata, expected_doublet_rate=est_mult_rate)
adata.obs[adata.obs["predicted_doublet"]].shape
adata.shape
min([x.obs["n_genes_by_counts"].min() for x in cutoff_df["adata_filt"]])


# %%
cutoffs_out_2 = []
for sample_id, co in zip(matrix_dir.iterdir(), cutoffs_out):
    these_cutoffs = {}
    these_cutoffs["sample_id"] = sample_id.name.replace("__outs", "")
    for k, v in co.items():
        for kk, vv in v.items():
            if kk == "limits":
                these_cutoffs[f"{k}__{kk}_min"] = vv[0]
                these_cutoffs[f"{k}__{kk}_max"] = vv[1]
            else:
                these_cutoffs[f"{k}__{kk}"] = vv
    cutoffs_out_2.append(these_cutoffs)


# %%
hmdf = pd.DataFrame(cutoffs_out_2)
hmdf["pct_counts_mt__cells_include_p"] = hmdf["pct_counts_mt__cells_include"] / (
    hmdf["pct_counts_mt__cells_include"] + hmdf["pct_counts_mt__cells_exclude"]
)
hmdf["n_genes_by_counts__cells_include_p"] = hmdf[
    "n_genes_by_counts__cells_include"
] / (
    hmdf["n_genes_by_counts__cells_include"] + hmdf["n_genes_by_counts__cells_exclude"]
)
hmdf["total_counts__cells_include_p"] = hmdf["total_counts__cells_include"] / (
    hmdf["total_counts__cells_include"] + hmdf["total_counts__cells_exclude"]
)
hmdf["merged__cells_include_p"] = hmdf["merged__cells_include"] / (
    hmdf["merged__cells_include"] + hmdf["merged__cells_exclude"]
)
hmdf = hmdf.set_index("sample_id")

# %%
importlib.reload(pmbip)
t = pmbip.Theme()
hmcol = palettable.scientific.sequential.Acton_20_r.mpl_colormap
mat = hmdf[
    [
        "pct_counts_mt__cells_include_p",
        "n_genes_by_counts__cells_include_p",
        "total_counts__cells_include_p",
        "merged__cells_include_p",
    ]
]
mat = mat.sort_values("sample_id")
panel = pmbip.Paneler(nrow=1, ncol=1, figsize=(8, 16))
im = pmbip.heatmap(
    mat, ax=panel.next_ax(), xlab="p(cells kept) by metric", ylab="Sample", cmap=hmcol
)
t.apply_to(panel.current_ax)
panel.current_ax.set_xticklabels(labels=panel.current_ax.get_xticklabels(), rotation=65)
panel.fig.colorbar(im, ax=panel.current_ax)
panel.fig.savefig("/home/amsesk/figures/coculture/qc_hm.png")


# %%
importlib.reload(scp)
scp.Curve.from_xy(x=np.array([1, 2, 3]), y=np.array([2, 3, 4]))
print(sample_name)
dydx_2_eq0 = com.crosses_zero_at(V=com.dens_rel_dydx_second[1])
com.threshold_intervals(V=com.dens_rel_dydx_second.y, threshold=0.0)


# %%
importlib.reload(scp)
importlib.reload(pmbip)
# x = np.linspace(-30,30,1000)
# y = np.array([xx**3 for xx in x])
# test_curve = scp.Curve(x,y)
panel = pmbip.Paneler(ncol=1, nrow=4, figsize=(8, 6))

d0 = com.get_relative_curve(derivative_order=0)
d1 = com.get_relative_curve(derivative_order=1)
d2 = com.get_relative_curve(derivative_order=2)

si = com.selection_interval(threshold=1.0, derivative_order=1)

panel.next_ax().plot(*d0, linestyle="-", c="r")
panel.current_ax.axvspan(xmin=d0.x[si.min()], xmax=d0.x[si.max()], color=cols[4])
# scp.axvspans_from_intervals(panel.current_ax, [i.as_index_of(d0.x) for i in intervals], colors= cols)
panel.next_ax().plot(*d1, linestyle="-", c="r")
panel.next_ax().plot(*d2, linestyle="-", c="r")

panel.fig.savefig(figout.joinpath(f"{sample_name}_qc_kde_diagnostic.png"))

# %%

np.diff(np.linspace(0, 10, 11)).shape
panel.next_ax().plot(*com._dydx_xy(order=1), marker=".", markersize=0.5, c="r")
panel.next_ax().plot(*com._dydx_xy(order=2), marker=".", markersize=0.5, c="r")
# com.threshold_intervals(V=com.dens_rel_dydx_second.y, threshold=0)

# panel.current_ax.axvline(x=com.dens_rel_x[dydx_2_eq0], linewidth=0.2, c="black")
# panel.next_ax().plot(*com.dens_rel_dydx_first, marker=".", markersize=0.5)
# pmbip.axvlines(panel.current_ax, com.dens_rel_x[dydx_2_eq0], linewidth = "0.25", color="blue")
# panel.next_ax().plot(*com.dens_rel_dydx_second, marker=".", markersize=0.5)
# pmbip.axvlines(panel.current_ax, com.dens_rel_x[dydx_2_eq0], linewidth = "0.25", color="blue")
# panel.fig.suptitle(sample_name)
# adata_gex_qc = scp.std_qc_gex(adata.copy(), sample_suffix = sample_name, mt_prefix = "MT-", plot_save_path=figout, automatic_cutoff_opts = {"threshold": 3.0})

# %%
adata_gex_qc.obs.columns
key = "n_genes_by_counts"
key = "pct_counts_mt"
adata_gex_qc.obs[key].min()
adata_gex_qc.obs[key].max()

# %%
[x.name.split("__")[0] for x in matrix_dir.iterdir()]
adata = pmbi.anndata.io.read_matrix()

# %%
# %%
# %%
np.where(adata_gex_qc.var.mt)
adata_gex_qc.write(f"{args.suffix}_std_gex_qc.h5ad")
gex = scp.std_gex(adata_gex_qc, sample_suffix=args.suffix)
gex.write(f"{args.suffix}_std_gex.h5ad")
# %%

[0, 1, 2, 3] * 3
