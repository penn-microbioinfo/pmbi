import argparse
import functools
import importlib
import re
from pathlib import Path

import numpy as np
from tqdm import tqdm
import palettable
import pandas as pd
import scanpy as sc
import scipy.stats
from joblib import Parallel, delayed
import anndata

import pmbi.anndata.io
import pmbi.plotting as pmbip
import pmbi.wrappers.scanpy as scp
from pmbi.cellranger.globals import ChromiumNextGEMSingleCell5primeHTv2__multiplet_rates

importlib.reload(scp)

sc._settings.ScanpyConfig.n_jobs = 1
tqdm.pandas()

# %% ArgumentParser
parser = argparse.ArgumentParser()
parser.add_argument(
    "-m",
    "--matrix",
    required=True,
    action="store",
    help="Path to 10X filtered_feature_bc_matrix.h5",
)
parser.add_argument(
    "-s",
    "--suffix",
    required=True,
    action="store",
    help="Samlple suffix to name output files and plots.",
)
parser.add_argument(
    "-o", "--figout", default=".", action="store", help="Directory to output figures."
)
# parser.add_argument("-", "--", required = True, action = "store", help = "")
args = parser.parse_args()


# matrix_dir = Path("/home/amsesk/super2/cellranger/HPAP-135_CC_1__outs/per_sample_outs/HPAP-135_CC_1/count/sample_filtered_feature_bc_matrix.h5")
# %%
importlib.reload(scp)
figout = Path("/home/amsesk/figures/coculture/qc_kde_figures/png")
matrix_dir = Path("/home/amsesk/super2/cellranger/")
mt_prefix = "MT-"
cols = (palettable.cartocolors.qualitative.Pastel_10.mpl_colors) * 10


# %%
def _compute(o, mt_prefix, cols, figout, no_lower=None, no_upper=None):
    if no_lower is None:
        no_lower = []
    if no_upper is None:
        no_upper = []
    sample_name = o.name.split("__")[0]
    matrix_path = o.joinpath(
        f"per_sample_outs/{sample_name}/count/sample_filtered_feature_bc_matrix.h5"
    )
    adata = pmbi.anndata.io.read_matrix(matrix_path, gex_only=True)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    cutoffs = {}
    cutoffs["sample_id"] = sample_name
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
    # panel.fig.suptitle(sample_name)
    # panel.fig.savefig(figout.joinpath(f"{sample_name}_qc_kde_diagnostic.png"))
    return cutoffs


# %%

fs = list(matrix_dir.iterdir())
cutoffs_out = Parallel(n_jobs=16)(
    delayed(_compute)(o, mt_prefix, cols, figout, no_lower=["pct_counts_mt"], no_upper=["n_genes_by_counts", "total_counts"]) for o in fs
)

# %%
dfdict = {}
for key in ["sample_id", "adata", "pct_counts_mt", "n_genes_by_counts", "total_counts"]:
    dfdict[key] = [x[key] for x in cutoffs_out]

cutoff_df = pd.DataFrame(dfdict).sort_values("sample_id").set_index("sample_id")
pd.set_option("display.width", 1000)


######################################
# %% Multiplet rate estimates
######################################

# %% Get sample metadata for giving scruble estimated multiplet rates
coculture_meta = pd.read_excel(
    "/home/amsesk/super1/t1d-coculture/CoCulture_metadata.xlsx"
)
ncells_est = coculture_meta[coculture_meta["Library_Type"] == "RNA"][
    ["Sample_Name", "Cells_in_Sample_Est"]
].set_index("Sample_Name")

# %% Estimate multiplet rates based on 2-degree polynomial fit to 10X stated rates
mult_rates = ChromiumNextGEMSingleCell5primeHTv2__multiplet_rates
mult_rate_coef = np.polyfit(
    mult_rates["n_cells_loaded"], mult_rates["multiplet_rate"], deg=2
)
mult_rate_p = np.poly1d(mult_rate_coef)
est_mult_rates = ncells_est.assign(
    est_mult_rate=lambda r: mult_rate_p(r["Cells_in_Sample_Est"])
)["est_mult_rate"]
cutoff_df["est_mult_rate"] = est_mult_rates

assert (
    cutoff_df["est_mult_rate"] == est_mult_rates.loc[cutoff_df["est_mult_rate"].index]
).all()


# %%
def apply_cutoffs(row, cutoff_columns):
    adata_filt = row["adata"]
    for qckey in cutoff_columns:
        adata_filt = adata_filt[
            (adata_filt.obs[qckey] >= row[qckey].loc["limits"][0])
            & (adata_filt.obs[qckey] <= row[qckey].loc["limits"][1])
        ].copy()
        adata_filt.uns[f"cutoffs__{qckey}"] = row[qckey].to_dict()
    row["adata_filt"] = adata_filt
    return row

# %%
def expertimentName_to_subcategory(en):
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

# %%
def other_preproc(sample_id, adata, est_mult_rate, min_cells, meta2merge):
    assert len(meta2merge.index) == len(meta2merge.index.unique())
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.scrublet(adata, expected_doublet_rate=est_mult_rate)
    adata.obs["original_barcode"] = adata.obs.index.to_list()
    adata.obs.index = pd.Index([f"{bc}__{sample_id}" for bc in adata.obs.index])
    adata.obs["sample_id"] = sample_id
    adata.obs = pd.merge(
        left=adata.obs,
        right=meta2merge,
        how="left",
        left_on="sample_id",
        right_index=True,
    )
    return adata


# %%
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
            expertimentName_to_subcategory
        )
    )
    .drop(columns=["Library_Type"])
    .set_index("sample_id")
)
meta_to_obs = meta_to_obs[~meta_to_obs["Experiment_subcategory"].isna()]

# %%
cutoff_df = cutoff_df.progress_apply(apply_cutoffs, axis=1, args=(["pct_counts_mt", "n_genes_by_counts", "total_counts"],))

# %%
adata_preproc = Parallel(n_jobs=32, verbose=10)(
    delayed(other_preproc)(sample_id, adata, est_mult_rate, min_cells, meta_to_obs)
    for sample_id, adata, est_mult_rate, min_cells in [
        (r.Index, r.adata_filt, r.est_mult_rate, 3)
        for r in cutoff_df.itertuples()
    ]
)

######################################
# %% Concatenate filtered AnnDatas and output for integration
######################################
for adata in adata_preproc:
    assert not any([vn.startswith("Hu.") for vn in adata.var_names])


pmbi.anndata.io.write_h5ad_multi({x.obs["sample_id"][0]: x for x in adata_preproc}, suffix="preproc", outdir="/home/amsesk/super2/h5ad/indiv_preproc")
adata_preproc[0].uns["cutoffs__pct_counts_mt"].to_dict()
combined = anndata.concat(adata_preproc, axis=0, join="outer")
combined.write_h5ad("/home/amsesk/super2/h5ad/combined_raw_counts.h5ad")
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
