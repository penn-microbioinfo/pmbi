import copy
import gc
import os
import pathlib
import re
from typing import Optional

import anndata
import matplotlib as mpl
import matplotlib.axis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scirpy as ir
import seaborn as sns

# local imports
import pmbi.io

global MT_CUT
global NFEAT_CUT
global NRNA_CUT

MT_CUT = (0.0, 0.95)
NFEAT_CUT = (0.025, 0.90)
NRNA_CUT = (0.025, 0.90)


class Paneler(object):
    def __init__(
        self,
        nrow,
        ncol,
        output_prefix=None,
        figsize=(3, 3),
        layout="tight",
        format="tiff",
        dpi=400,
    ):
        Paneler.theme()
        self.nrow = nrow
        self.ncol = ncol
        self.output_prefix = output_prefix
        self.figsize = figsize
        self.layout = layout
        self.format = format
        self.dpi = dpi
        self.fig, self.axs = self._new_fig()
        self.panel_idx = 0
        self.image_num = 1
        self.current_ax = self._get_ax((0, 0))

    def _new_fig(self):
        return plt.subplots(
            nrows=self.nrow, ncols=self.ncol, figsize=self.figsize, layout=self.layout
        )

    @staticmethod
    def theme():
        plt.rcdefaults()
        plt.rcParams.keys()
        plt.rcParams.update(
            {  # "font.sans-serif": "Arial",
                "font.size": 5,
                "figure.dpi": 300,
                "axes.titlesize": 6,
            }
        )
        plt.rcParams["legend.frameon"] = False

    @staticmethod
    def subplot_idx_to_pos(nrow, ncol, idx):
        r = int(np.floor(np.true_divide(idx, ncol)))
        rem = idx % ncol
        c = rem
        if r >= nrow:
            r = r % nrow
        return (r, c)

    @staticmethod
    def pos_to_subplot_idx(nrow: int, ncol: int, coords: tuple[int, int]) -> int:
        row = coords[0]
        col = coords[1]
        if row > nrow - 1 or col > ncol - 1:
            raise IndexError(
                f"coordinates {coords} outside of limits of panel with shape ({nrow}, {ncol})"
            )
        idx = (row * ncol) + (col)
        return idx

    def _get_ax(self, coords: tuple[int, int]):
        if self.panel_idx_out_of_bounds():
            raise IndexError("paneler.panel_idx >= panel_dimensions")
        else:
            if self.nrow * self.ncol == 1:
                ax = self.axs
            elif self.nrow == 1 or self.ncol == 1:
                if self.nrow == 1:
                    ax = self.axs[coords[1]]
                else:
                    ax = self.axs[coords[0]]
            else:
                r, c = coords
                ax = self.axs[r, c]
            return ax

    def panel_idx_out_of_bounds(self):
        if self.panel_idx > (self.nrow * self.ncol) - 1:
            return True
        else:
            return False

    def next_ax(self):
        if self.panel_idx_out_of_bounds():
            self._advance_image()
        coords = Paneler.subplot_idx_to_pos(self.nrow, self.ncol, self.panel_idx)
        print(self.image_num, coords)
        self.current_ax = self._get_ax(coords)
        self.panel_idx += 1
        return self.current_ax

    """
    def insert(self, plt_function, **kwargs):
        if self.nrow*self.ncol == 1:
            plt_function(ax=self.ax, **kwargs)
        elif self.nrow or self.ncol == 1:
            plt_function(ax=self.ax[panel_idx], **kwargs)
        else:
            row,col = Paneler.subplot_idx_to_pos(self.nrow, self.ncol, self.panel_idx)
            plt_function(ax=self.ax[row,col], **kwargs)
    """

    def subplots_adjust(self, **kwargs):
        self.fig.subplots_adjust(**kwargs)

    def _advance_image(self):
        self.savefig(f"{self.output_prefix}_{self.image_num}.{self.format}")
        self.fig, self.axs = self._new_fig()
        self.image_num += 1
        self.panel_idx = 0

    def savefig(self, filename: Optional[str] = None):
        if filename is None:
            filename = f"{self.output_prefix}_{self.image_num}.{self.format}"
        self.subplots_adjust(hspace=0.5, wspace=0.5)
        self.fig.savefig(filename, bbox_inches="tight")
        plt.close(self.fig)
        plt.clf()
        gc.collect()


def adata_add_gct(adata: anndata.AnnData, gct_path: os.PathLike, rsuffix: str):
    gct = pmbi.io.read_gct(gct_path)
    adata.obs = adata.obs.join(gct, rsuffix=rsuffix)


def adata_to_gct(adata, outpath, layer=None):
    adata = adata.copy()
    if layer is not None:
        adata.X = adata.layers[layer]
    adata_df = adata.to_df().transpose().reset_index()
    adata_df = adata_df.rename(columns={"index": "NAME"})
    adata_df.insert(1, "Description", adata_df["NAME"])
    with open(outpath, "w") as out:
        out.write("#1.2\n")
        out.write(f"{len(adata.var_names)}\t{len(adata.obs_names)}\n")
        adata_df.to_csv(out, sep="\t", index=False)


def obs_add_orig_barcode(adata):
    p = re.compile("([ATCG]+)([-][0-9])([-][0-9])*")
    orig_barcodes = []
    for obsname in adata.obs_names:
        m = re.match(p, obsname)
        if m is not None:
            orig_barcodes.append(f"{m.groups()[0]}{m.groups()[1]}")
    adata.obs["orig_barcode"] = orig_barcodes
    return adata


def obs_names_unique(adata: anndata.AnnData) -> bool:
    return len(adata.obs_names) == len(adata.obs_names.unique())


def var_names_unique(adata: anndata.AnnData) -> bool:
    return len(adata.var_names) == len(adata.var_names.unique())


def combine_adatas(adatas: dict) -> anndata.AnnData:
    idented = {}
    for sample, adata in adatas.items():
        print(sample)
        adata.obs["orig_ident"] = sample
        adata.obs_names = pd.Series(adata.obs_names).str.cat(
            others=[str(sample)] * len(adata.obs_names), sep="_"
        )
        idented[sample] = adata
        if len(adata.obs_names) != len(adata.obs_names.unique()):
            raise
    return anndata.concat(idented.values(), axis=0, join="outer")


def calc_umi_per_bc(adata: anndata.AnnData) -> pd.DataFrame:
    umi_per_bc = adata.X.sum(axis=1)[:, 0].A1
    umi_per_bc.sort()
    umi_per_bc = pd.DataFrame(
        {
            "rank": pd.Series(range(1, len(umi_per_bc) + 1)),
            "umi_count": pd.Series(umi_per_bc[::-1] + 1),
        }
    )
    return umi_per_bc


def ax_umi_per_bc(adata: anndata.AnnData, axs_title="barcode umi counts", xmax=5e4):
    umi_per_bc = calc_umi_per_bc(adata)
    ax = sns.lineplot(umi_per_bc, x="rank", y="umi_count")
    ax.set(yscale="log", ylim=(0, max(umi_per_bc["umi_count"]) * 2.5), xlim=(0, xmax))
    plt.title(axs_title)
    return ax


def subplot_idx_to_pos(nrow, ncol, idx):
    r = int(np.floor(np.true_divide(idx, ncol)))
    rem = idx % ncol
    c = rem
    if r >= nrow:
        r = r % nrow
    return (r, c)


def percent_expressing(
    adata: anndata.AnnData, groupby: str, var_names: list, layer: str
):
    subadata = adata[:, var_names]
    frame = pd.DataFrame(subadata.layers[layer].toarray())
    frame.index = adata.obs_names
    frame[groupby] = subadata.obs[groupby]
    print(frame)
    print(frame.groupby(groupby).agg(lambda x: len([v for v in x if v != 0.0])))


def cluster_leiden(adata, resolutions):
    for r in resolutions:
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}")
        # adata.obs[f"leiden_{r}"] = adata.obs[f"leiden_{r}"].astype(str) + "__" + adata.obs[f"orig_ident"].astype(str)


# implementation of sc.pp.normalize_geometric(adt) from https://github.com/scverse/scanpy/pull/1117
def clr_normalize(adata=anndata.AnnData):
    counts_nda = adata.X.toarray()
    adata.X = scipy.sparse.csr_matrix(
        np.divide(counts_nda, scipy.stats.mstats.gmean(counts_nda + 1, axis=0))
    )
    sc.pp.log1p(adata)
    return adata


def subset(adata: anndata.AnnData, obs_name: str, bounds: tuple):
    lower, upper = np.quantile(adata.obs[obs_name], bounds)
    return adata[(adata.obs[obs_name] >= lower) & (adata.obs[obs_name] <= upper)]


def add_var_column_by_index(adata, convert, new_colname):
    adata = adata.copy()
    adata.var[new_colname] = [convert[n] for n in adata.var.index]
    return adata


def std_qc_gex(
    adata: anndata.AnnData,
    sample_suffix: str = "sample",
    min_cells: int = 3,
    min_genes: int = 200,
    plot: bool = True,
    plot_save_path: os.PathLike = pathlib.Path("qc_violins.pdf"),
    plot_point_size: int = 1,
    mt_prefix: str = "mt-",
    manual_cutoffs: Optional[dict] = None,
):

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    if manual_cutoffs is None:
        cuts = {
            "pct_counts_mt": [
                np.quantile(adata.obs["pct_counts_mt"], MT_CUT[0]),
                np.quantile(adata.obs["pct_counts_mt"], MT_CUT[1]),
            ],
            "n_genes_by_counts": [
                np.quantile(adata.obs["n_genes_by_counts"], NFEAT_CUT[0]),
                np.quantile(adata.obs["n_genes_by_counts"], NFEAT_CUT[1]),
            ],
            "total_counts": [
                np.quantile(adata.obs["total_counts"], NRNA_CUT[0]),
                np.quantile(adata.obs["total_counts"], NRNA_CUT[1]),
            ],
        }
    else:
        cuts = manual_cutoffs

    print(cuts)

    adata.uns["qc_filter_ranges"] = copy.deepcopy(cuts)

    if plot:
        panel = Paneler(ncol=3, nrow=1, figsize=(9, 3))
        for key in cuts:
            sc.pl.violin(
                adata=adata,
                keys=key,
                jitter=0.4,
                multi_panel=False,
                save=False,
                size=plot_point_size,
                ax=panel.next_ax(),
            )
            panel.current_ax.axhline(y=cuts[key][0], color="r")
            panel.current_ax.axhline(y=cuts[key][1], color="r")
        panel.fig.suptitle(sample_suffix)
        panel.savefig(plot_save_path)

    adata = adata[
        (adata.obs["pct_counts_mt"] <= cuts["pct_counts_mt"][1])
        & (adata.obs["pct_counts_mt"] >= cuts["pct_counts_mt"][0])
    ]
    adata = adata[
        (adata.obs["n_genes_by_counts"] <= cuts["n_genes_by_counts"][1])
        & (adata.obs["n_genes_by_counts"] >= cuts["n_genes_by_counts"][0])
    ]
    adata = adata[
        (adata.obs["total_counts"] <= cuts["total_counts"][1])
        & (adata.obs["total_counts"] >= cuts["total_counts"][0])
    ]

    return adata


def manual_qc_cutoffs_to_dict(csv, sample_column):
    mancuts = pd.read_csv(csv, sep=",")
    keys = set([re.sub("([_]lower)|([_]upper)", "", c) for c in mancuts.columns])
    keys.remove(sample_column)
    d = {}
    for idx, row in mancuts.iterrows():
        d[row[sample_column]] = {}
        for k in keys:
            d[row[sample_column]][k] = [row[f"{k}_lower"], row[f"{k}_upper"]]
    return d


def std_gex(adata: anndata.AnnData, sample_suffix: str = "sample", plot: bool = True):

    adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)

    adata = adata[:, adata.var.highly_variable]

    adata.layers["totalcounts_pctmt_regressed_out"] = sc.pp.regress_out(
        adata, ["total_counts", "pct_counts_mt"], copy=True
    ).X

    adata.layers["scale_data"] = sc.pp.scale(adata, max_value=10, copy=True).X

    sc.tl.pca(adata, svd_solver="arpack")

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    sc.tl.umap(adata)
    print("making names unique")

    if plot:
        sc.pl.highly_variable_genes(adata, save=f"_{sample_suffix}")
        sc.pl.pca(adata, save=f"_{sample_suffix}")
        sc.pl.umap(adata, save=f"_{sample_suffix}")

    return adata


def std_qc_adt(
    adata: anndata.AnnData,
    sample_suffix: str = "sample",
    plot: bool = True,
    plot_save_path: os.PathLike = ".",
):

    sc.pp.filter_genes(adata, min_counts=1)

    adata.var["control"] = adata.var_names.str.startswith("Isotype")
    sc.pp.calculate_qc_metrics(
        adata,
        percent_top=(5, 10, 15),
        var_type="antibodies",
        qc_vars=("control",),
        inplace=True,
    )
    if plot:
        panel = Paneler(ncol=1, nrow=1, figsize=(3, 3))
        panel.fig = sns.jointplot(
            x="log1p_total_counts",
            y="n_antibodies_by_counts",
            data=adata.obs,
            kind="hex",
            norm=mpl.colors.LogNorm(),
        ).figure
        panel.savefig(os.path.join(plot_save_path, f"{sample_suffix}_adt1.png"))
        panel = Paneler(ncol=1, nrow=1, figsize=(3, 3))
        panel.fig = sns.jointplot(
            x="log1p_total_counts",
            y="log1p_total_counts_control",
            data=adata.obs,
            kind="hex",
            norm=mpl.colors.LogNorm(),
        ).figure
        panel.savefig(os.path.join(plot_save_path, f"{sample_suffix}_adt2.png"))
        plt.clf()
        plt.close()
    return adata


def std_adt_prepocess(adata: anndata.AnnData, sample_suffix: str = "sample"):

    adata = std_qc_adt(adata, sample_suffix=sample_suffix, plot=False)
    adata = clr_normalize(adata)
    adata.layers["scale_data"] = sc.pp.scale(adata, max_value=10, copy=True).X

    return adata


def ranked_genes_to_df(adata):
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    result.keys()
    ranked_ds = pd.DataFrame(
        {
            group + "_" + key: result[key][group]
            for group in groups
            for key in ["names", "pvals", "pvals_adj", "logfoldchanges"]
        }
    )
    return ranked_ds


# %%
def get_counts(
    adata: anndata.AnnData, layer: str | None = None
) -> scipy.sparse.csc_matrix:
    if layer is None:
        return adata.X
    else:
        return adata.layers[layer]

def get_count_frame(
    adata: anndata.AnnData, layer: str | None = None
) -> pd.DataFrame:
    if layer is None:
        return pd.DataFrame(adata.X.toarray(), index = adata.obs_names, columns = adata.var_names)
    else:
        return pd.DataFrame(adata.layers[layer].toarray(), index = adata.obs_names, columns = adata.var_names)

def mean_expression(
    adata: anndata.AnnData, groupby: str | None = None, layer: str = None
) -> pd.DataFrame:
    if groupby is None:
        counts = get_counts(adata, layer)
        meanexpr = counts.mean(axis=0).transpose().getA()[:, 0]
        meanexpr = pd.DataFrame({"all_cells": meanexpr}, index=adata.var_names)
    else:
        groups = adata.obs[groupby].unique()
        group_means = []
        for group in groups:
            group_means.append(
                get_counts(adata[adata.obs[groupby] == group], layer=layer)
                .mean(axis=0)
                .transpose()
                .getA()
            )
        meanexpr = pd.DataFrame(
            np.column_stack(group_means), index=adata.var_names, columns=groups
        )
    return meanexpr


if __name__ == "__main__":
    COLOR_GENES = ["CD3D", "CD4", "CD8A", "MS4A1", "CD19", "TRCF", "TCF7"]
    COLOR_AB = [
        "Hu.CD3_UCHT1",
        "Hu.CD4_RPA.T4",
        "Hu.CD8",
        "Hu.CD20_2H7",
        "Hu.CD19",
        "Hu.CD71",
    ]

    adata = sc.read_10x_h5(
        "/home/ubuntu/mnt/betts/aggr/135/outs/count/filtered_feature_bc_matrix.h5",
        gex_only=False,
    )
    adata.layers["counts"] = adata.X.copy()
    adata.var_names_make_unique()

    adata.var.feature_types.value_counts()
    gex = std_gex(adata[:, adata.var.feature_types == "Gene Expression"], plot=False)
    adt = std_adt(adata[:, adata.var.feature_types == "Antibody Capture"], plot=False)

    gex.x.std()
    adt.x.std()

    [x for x in gex.var_names if x == "cd19"]

    adt
    gex
    [x for x in list(adt.var_names)[0:130] if x.split(".")[1] in color_genes]

    import itertools

    comb = anndata.concat([gex, adt], axis=1, merge="first")
    comb
    sc.pl.umap(
        comb,
        color=[
            g
            for g in list(itertools.chain(*[x for x in zip(color_genes, color_ab)]))
            if g in comb.var_names
        ],
        save=true,
        ncols=2,
    )

    ad.var.feature_types

    np.quantile(ad.obs.pct_counts_mt, [1.0, 0.95, 0.90])
    min(ad.obs.pct_counts_mt)

    p13_final.obs.to_csv("p13_metadata.tsv", sep="\t")
    p13_final._raw.x[1:10, 1:10].toarray()
    p13_raw = p13_final._raw.to_adata()
    del p13_final

    len(p13_raw.var.columns)
    p13_raw.obs.columns
    p13_raw.obs.barcode_sample

    p13_raw.obs["orig_ident"] = p13_raw.obs.index
    p13_raw.obs["orig_ident"] = p13_raw.obs["orig_ident"].apply(
        lambda x: re.sub("[atcg]+[-][0-9]+[_]", "", x)
    )
    p13_raw.obs["old_index"] = p13_raw.obs.index
    p13_raw.var["feature"] = p13_raw.var.index
    p13_raw.obs = p13_raw.obs.reset_index()

    p13_raw.var.shape
    for oi in p13_raw.obs["orig_ident"].unique():
        which = p13_raw.obs.loc[p13_raw.obs["orig_ident"] == oi].index
        which_adata = p13_raw[which, :]
        scipy.io.mmwrite(f"p13_split/{oi}.mm", which_adata.x.transpose())
        which_adata.obs.to_csv(
            f"p13_split/{oi}.obs", sep="\t", header=false, index=false
        )
        which_adata.var.to_csv(
            f"p13_split/{oi}.var", sep="\t", header=false, index=false
        )
# %%
