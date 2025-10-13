from __future__ import annotations

import copy
import gc
import os
import pathlib
import re
from pathlib import Path
from typing import Callable, Optional, Tuple, Type, TypeVar

import anndata
import matplotlib as mpl
import matplotlib.axes
import matplotlib.axis
import matplotlib.pyplot as plt
import numpy as np
import palettable
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.signal
import scipy.stats
import scirpy as ir
import seaborn as sns

# local imports
import pmbi.io
import pmbi.plotting as pmbip
from pmbi.logging import streamLogger

global MT_CUT
global NFEAT_CUT
global NRNA_CUT

MT_CUT = (0.0, 0.95)
NFEAT_CUT = (0.025, 0.90)
NRNA_CUT = (0.025, 0.90)


# %%
def qc_kde(adata, key, x_len=1000):
    if key not in adata.obs.columns:
        raise ValueError(f"key not in obs: {key}")

    X = adata.obs[key]
    kde = scipy.stats.gaussian_kde(X)
    # x_values = np.linspace(min(X), max(X), int(np.ceil(max(X)-min(X))))
    x_values = np.linspace(min(X), max(X), x_len)
    y_values = kde(x_values)
    return (x_values, y_values)


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
    umi_per_bc = adata.X.sum(axis=1)[:, 0]
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


def cluster_leiden(adata, resolutions, **kwargs):
    for r in resolutions:
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", **kwargs)
        # adata.obs[f"leiden_{r}"] = adata.obs[f"leiden_{r}"].astype(str) + "__" + adata.obs[f"orig_ident"].astype(str)


# implementation of sc.pp.normalize_geometric(adt) from https://github.com/scverse/scanpy/pull/1117
def clr_normalize(adata: anndata.AnnData) -> anndata.AnnData:
    adata = adata.copy()
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


class Interval:
    def __init__(self, values: np.ndarray):
        self.values = values

    def min(self):
        return self.values.min()

    def max(self):
        return self.values.max()

    def as_index_of(self, X):
        return X[self.values]


T = TypeVar("T", bound="Points")


class Points:
    def __init__(self, xy: np.ndarray):
        assert xy.ndim == 2, "xy should be a 2-dimensional array"
        self.xy = xy

    def __getitem__(self, key):
        return Points(self.xy[key])

    def __setitem__(self, key, value):
        self.xy[key] = value

    @classmethod
    def from_xy(cls: Type[T], x: np.ndarray, y: np.ndarray) -> T:
        return cls(xy=np.column_stack((x, y)))

    @property
    def x(self):
        return self.xy[:, 0]

    @property
    def y(self):
        return self.xy[:, 1]

    def __iter__(self):
        yield self.x
        yield self.y


class Curve(Points):
    def __init__(self, xy: np.ndarray):
        super().__init__(xy)

    def axis_values(self, axis):
        if axis not in ["x", "y"]:
            raise ValueError(f"Axis must be in: {axis}")
        return getattr(self, axis)

    def critical_points(self):
        first_derivative = self.dydx(order=1)
        return first_derivative.crosses_zero_at()

    def maxima(self):
        dens_crits = self.critical_points()
        second_derivative = self.dydx(order=2)
        return np.array(
            [xx for xx in dens_crits if second_derivative.y[xx] < 0], dtype=np.int64
        )

    def minima(self):
        dens_crits = self.critical_points()
        second_derivative = self.dydx(order=2)
        return np.array(
            [xx for xx in dens_crits if second_derivative.y[xx] > 0], dtype=np.int64
        )

    def dydx(self, order=1) -> Curve:
        if order == 0:
            return self
        else:
            dydx = np.gradient(self.y, self.x)

            return Curve.from_xy(self.x, dydx).dydx(order=(order - 1))

    def crosses_threshold_at(self, threshold, V_ignore_sign=True, axis="y"):
        V = self.axis_values(axis=axis)
        if V_ignore_sign:
            V = abs(V)
        return np.where(np.diff(np.sign(V - threshold)))[0]

    def crosses_zero_at(self, axis="y"):
        if axis not in ["x", "y"]:
            raise ValueError(f"Axis must be in: {axis}")
        return self.crosses_threshold_at(threshold=0.0, V_ignore_sign=False, axis=axis)

    def threshold_intervals(self, threshold, V_ignore_sign=True, axis="y"):
        V = self.axis_values(axis=axis)
        if threshold == 0.0:
            return [
                Interval(s)
                for s in np.split(range(0, len(V)), self.crosses_zero_at(axis=axis))
            ]
        else:
            return [
                Interval(s)
                for s in np.split(
                    range(0, len(V)),
                    self.crosses_threshold_at(
                        threshold, V_ignore_sign=V_ignore_sign, axis=axis
                    ),
                )
            ]


class CutoffMaker:
    def __init__(self, X):
        self.X = X


class KdeCutoffMaker(CutoffMaker):
    def __init__(
        self, X, gradient_threshold=0.25, dens_x_len=100000, endpoint_buffer=3
    ):
        super().__init__(X)
        self.gradient_threshold = gradient_threshold
        self.dens_x_len = dens_x_len
        self.endpoint_buffer = endpoint_buffer
        self.dens = self._make_density_curve()
        self.dens_rel = Curve.from_xy(
            x=self.dens.x / self.dens.x.max(), y=self.dens.y / self.dens.y.max()
        )

    def _make_density_curve(self):
        kde = scipy.stats.gaussian_kde(self.X)

        # Buffer endpoints in preparation for taking derivatives with np.gradient, if needed
        if self.endpoint_buffer > 0:
            dens_x_unbuffered = np.linspace(min(self.X), max(self.X), self.dens_x_len)
            dens_x_increment = np.diff(dens_x_unbuffered)[0]
            buffered_min = min(self.X) - (self.endpoint_buffer * dens_x_increment)
            buffered_max = max(self.X) + (self.endpoint_buffer * dens_x_increment)
            dens_x = np.linspace(
                buffered_min,
                buffered_max,
                (self.dens_x_len + (2 * self.endpoint_buffer)),
            )
        else:
            dens_x = np.linspace(min(self.X), max(self.X), self.dens_x_len)

        dens_y = kde(dens_x)
        return Curve.from_xy(dens_x, dens_y)

    # Returns a Curve based the relativized density curve on the derivative order. Buffered endpoints are truncated if they exist
    def get_relative_curve(self, derivative_order=0):
        x, y = self.dens_rel.dydx(order=derivative_order)
        if self.endpoint_buffer > 0:
            ret = Curve.from_xy(
                x[self.endpoint_buffer : -(self.endpoint_buffer)],
                y[self.endpoint_buffer : -(self.endpoint_buffer)],
            )
        else:
            ret = Curve.from_xy(x, y)

        assert (
            ret.x.shape[0] == self.dens_x_len
        ), "Shape of gotten curve is not as expected"
        return ret

    # Returns a Curve based the absolute density curve on the derivative order. Buffered endpoints are truncated if they exist
    def get_curve(self, derivative_order=0):
        x, y = self.dens.dydx(order=derivative_order)
        if self.endpoint_buffer > 0:
            ret = Curve.from_xy(
                x[self.endpoint_buffer : -(self.endpoint_buffer)],
                y[self.endpoint_buffer : -(self.endpoint_buffer)],
            )
        else:
            ret = Curve.from_xy(x, y)

        assert (
            ret.x.shape[0] == self.dens_x_len
        ), "Shape of gotten curve is not as expected"
        return ret

    def selection_interval(self):
        relcurve_d0 = self.get_relative_curve(derivative_order=0)
        relcurve_d1 = self.get_relative_curve(derivative_order=1)
        t_intervals = relcurve_d1.threshold_intervals(threshold=self.gradient_threshold)
        y_max = 1.0
        most_prominent_interval_idx = np.where(
            [y_max in ti.as_index_of(relcurve_d0.y) for ti in t_intervals]
        )[0]
        assert (
            len(most_prominent_interval_idx) == 1
        ), "More than 1 interval that includes the max y"
        center_int_idx = most_prominent_interval_idx[0]
        if len(t_intervals[center_int_idx].values) < 3:
            raise ValueError(
                "Center interval around max y has length < 3, this can screw up interval selection. Increase dens_x_len."
            )
        # print(center_int_idx)
        intervals_to_cat = [
            ti.values
            for ti in t_intervals[(center_int_idx - 1) : (center_int_idx + 1) + 1]
        ]
        # print(intervals_to_cat)
        # print(intervals_to_cat.shape)
        selection_interval = np.concatenate(intervals_to_cat)
        return selection_interval


def axvspans_from_intervals(
    ax: matplotlib.axes.Axes, intervals: list[Interval], colors: list[str], **kwargs
):
    for idx, interval in enumerate(intervals):
        ax.axvspan(
            xmin=interval.min(), xmax=interval.max(), color=colors[idx], **kwargs
        )


def std_qc_gex(
    adata: anndata.AnnData,
    sample_suffix: str = "sample",
    min_cells: int = 3,
    min_genes: int = 200,
    plot: bool = True,
    plot_save_path: os.PathLike = pathlib.Path("qc_violins.pdf"),
    plot_point_size: int = 1,
    mt_prefix: str = "mt-",
    automatic_cutoff_opts=None,
    manual_cutoffs: Optional[dict] = None,
):

    AUTOMATIC_CUTOFF_OPTS_DEFAULT = {"threshold": 0.50, "x_len": 100000}
    AUTOMATIC_CUTOFF_OPTS = {}

    logger = streamLogger("std_qc_gex")
    plot_save_path = Path(plot_save_path)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    AUTOMATIC_CUTOFF_OPTS = dict(AUTOMATIC_CUTOFF_OPTS_DEFAULT)
    if automatic_cutoff_opts is not None:
        AUTOMATIC_CUTOFF_OPTS.update(automatic_cutoff_opts)

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
        panel = pmbip.Paneler(ncol=3, nrow=4, figsize=(8, 6))
        threshold = AUTOMATIC_CUTOFF_OPTS["threshold"]
        x_len = AUTOMATIC_CUTOFF_OPTS["x_len"]
        for keyidx, key in enumerate(cuts):
            x, y = qc_kde(adata, key, x_len=x_len)
            y_rel = y / y.max()
            x_rel = x / x.max()
            dy_dx = np.gradient(y_rel, x_rel)
            dy_dx_eq0 = np.where(np.diff(np.sign(dy_dx)))[0]
            dy_dx_crosses_thresh = np.where(np.diff(np.sign(abs(dy_dx) - threshold)))[0]
            dy_dx_2 = np.gradient(dy_dx, x_rel)
            inflect_idx = np.where(np.diff(np.sign(dy_dx_2)) != 0)[0]
            maxima = [xx for xx in dy_dx_eq0 if dy_dx_2[xx] < 0]
            minima = [xx for xx in dy_dx_eq0 if dy_dx_2[xx] > 0]
            # x_unit_width = (x.max()-x.min())/x_len
            # height = 1/(x_len*10)
            # width = np.ceil( x_len/100 )

            # These are intervals separated by the points where dy/dx crosses the thershold value
            intervals = np.split(range(0, x_len), dy_dx_crosses_thresh)

            # intervals_idx = np.split(range(0, x_len), maxima_minima_idx)
            max_y = y_rel.max()
            assert max_y == 1.0, "Maximum of relativized y values is not 1.0"

            most_prominent_interval_idx = np.where(
                [max_y in y_rel[xx] for xx in intervals]
            )[0]
            assert (
                len(most_prominent_interval_idx) == 1
            ), "More than 1 interval between inflection points that include the max y"
            center_int_idx = most_prominent_interval_idx[0]
            if len(intervals[center_int_idx]) < 3:
                logger.warning(
                    "Center interval around max y has length < 3, this can screw up interval selection. Increase x_len."
                )
            selection_interval = np.concatenate(
                intervals[(center_int_idx - 1) : (center_int_idx + 1) + 1]
            )
            co = (selection_interval.min(), selection_interval.max())
            print(x[co[0]], x[co[1]])
            # print(center_int_idx)

            cols = (palettable.cartocolors.qualitative.Pastel_10.mpl_colors) * 10
            select_col = palettable.cartocolors.qualitative.Antique_10.mpl_colors

            currax = panel._get_ax([0, keyidx])
            currax.plot(x_rel, y_rel, marker=".", markersize=1)
            currax.axvspan(
                xmin=x_rel[selection_interval.min()],
                xmax=x_rel[selection_interval.max()],
                color=select_col[7],
            )
            currax.set_title(key)

            for ii, yval in enumerate([y_rel, dy_dx, dy_dx_2]):
                currax = panel._get_ax([ii + 1, keyidx])
                for iii, interval in enumerate(intervals):
                    currax.plot(
                        x_rel[interval],
                        yval[interval],
                        marker=".",
                        markersize=1,
                        color=cols[iii],
                    )
                    currax.axhline(y=0, linewidth=0.5, c="black", linestyle="--")
                    # currax.axvspan(xmin = x_rel[interval.min()], xmax=x_rel[interval.max()], color=cols[iii])
                    currax.axvline(x=x_rel[interval.min()], linewidth=0.2, c="black")
                for m in minima:
                    currax.axvline(x=x_rel[m], linewidth=0.2, c="red")
                for m in maxima:
                    currax.axvline(x=x_rel[m], linewidth=0.2, c="blue")

        panel.fig.suptitle(sample_suffix)
        panel.fig.savefig(
            plot_save_path.joinpath(f"{sample_suffix}_qc_kde_diagnostic.png")
        )

    # adata = adata[
    #     (adata.obs["pct_counts_mt"] <= cuts["pct_counts_mt"][1])
    #     & (adata.obs["pct_counts_mt"] >= cuts["pct_counts_mt"][0])
    # ]
    # adata = adata[
    #     (adata.obs["n_genes_by_counts"] <= cuts["n_genes_by_counts"][1])
    #     & (adata.obs["n_genes_by_counts"] >= cuts["n_genes_by_counts"][0])
    # ]
    # adata = adata[
    #     (adata.obs["total_counts"] <= cuts["total_counts"][1])
    #     & (adata.obs["total_counts"] >= cuts["total_counts"][0])
    # ]

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


def get_count_frame(adata: anndata.AnnData, layer: str | None = None) -> pd.DataFrame:
    if layer is None:
        return pd.DataFrame(
            adata.X.toarray(), index=adata.obs_names, columns=adata.var_names
        )
    else:
        return pd.DataFrame(
            adata.layers[layer].toarray(),
            index=adata.obs_names,
            columns=adata.var_names,
        )


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
