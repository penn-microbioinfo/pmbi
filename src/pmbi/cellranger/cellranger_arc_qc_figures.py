from typing import Optional, Union

import matplotlib.axes
import mudata
import muon.atac as atac
import numpy as np
import pandas as pd
import pysam
import snapatac2.genome
from sklearn.neighbors import KernelDensity

import pmbi.plotting as pmbip
from pmbi.wrappers.scanpy import calc_umi_per_bc


# %%
def tss_enrichment(
    mdata: mudata.MuData,
    n_tss: int = 2000,
    extend_upstream: int = 1000,
    extend_downstream: int = 1000,
):
    fragments = pysam.TabixFile(
        mdata["atac"].uns["files"]["fragments"], parser=pysam.asBed()
    )
    feats = atac.tl.get_gene_annotation_from_rna(mdata)
    feats_filt = feats[feats.Chromosome.isin(fragments.contigs)]
    feats_with_window = []
    for row in feats_filt.itertuples():
        try:
            _ = fragments.fetch(
                row.Chromosome, row.Start - extend_upstream, row.End + extend_downstream
            )
            feats_with_window.append(row.gene_id)
        except ValueError:
            pass
    feats_filt = feats_filt.loc[feats_with_window, :]
    return atac.tl.tss_enrichment(
        mdata,
        features=feats_filt,
        n_tss=n_tss,
        extend_upstream=extend_upstream,
        extend_downstream=extend_downstream,
    )


def tss_enrichment_plot(
    mdata: mudata.MuData,
    ax: matplotlib.axes.Axes,
    n_tss: int = 2000,
    extend_upstream: int = 1000,
    extend_downstream: int = 1000,
    color="red",
    fill_alpha=0.2,
):
    tsse = tss_enrichment(
        mdata=mdata,
        n_tss=n_tss,
        extend_upstream=extend_upstream,
        extend_downstream=extend_downstream,
    )
    x = tsse.var["TSS_position"]
    means = tsse.X.mean(axis=0)
    ax.plot(x, means, linewidth=1, color=color)
    ax.fill_between(x, means, alpha=fill_alpha, color=color)
    ax.axvline(x=0.0, linestyle="dashed", linewidth=0.5, color="black")
    ax.set_xlabel("Position relative to TSS")
    ax.set_ylabel("Enrichment")
    ax.set_title("TSS Enrichment")


# %%
def fragment_size_distr(
    mdata: mudata.MuData,
    genome: snapatac2.genome.Genome,
    sample_size: int = 50000,
    region_str: Optional[str] = None,
    fragments: Optional[pd.DataFrame] = None,
    modal_key: str = "atac",
    kde_bandwidth: float = 1.0,
) -> KernelDensity:
    if fragments is None:
        if region_str is None:
            chrom_sizes = pd.Series(genome.chrom_sizes).sort_values(ascending=False)
            largest_chrom = chrom_sizes.index[0]
            region_str = f"{largest_chrom}:1-{chrom_sizes.loc[largest_chrom]}"
        frags = atac.tl.fetch_regions_to_df(
            fragment_path=mdata[modal_key].uns["files"]["fragments"],
            features=region_str,
        )
    else:
        frags = fragments
    frags["fragment_size"] = frags["End"] - frags["Start"]
    sidx = np.random.choice(
        np.arange(0, frags.shape[0]), size=sample_size, replace=False
    )
    fs = frags["fragment_size"].to_numpy()[sidx, np.newaxis]
    return KernelDensity(kernel="gaussian", bandwidth=kde_bandwidth).fit(fs)


# %%
def fragment_size_distr_plot(
    mdata: mudata.MuData,
    genome: snapatac2.genome.Genome,
    ax: matplotlib.axes.Axes,
    sample_size: int = 50000,
    region_str: Union[str, None] = None,
    modal_key: str = "atac",
    kde_bandwidth: float = 1.0,
    xmax: float = 1000.0,
    xlabel: str = "fragment size (bp)",
    ylabel: str = "density",
    **kwargs,
) -> None:
    kde = fragment_size_distr(
        mdata=mdata,
        genome=genome,
        sample_size=sample_size,
        region_str=region_str,
        modal_key=modal_key,
        kde_bandwidth=kde_bandwidth,
    )
    x_plot = np.arange(0, xmax)[:, np.newaxis]
    y_plot = kde.score_samples(x_plot)
    ax.plot(x_plot[:, 0], np.exp(y_plot), **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Fragment size distribution")


# %%
def barcode_totals(mdata: mudata.MuData, modal_key: str):
    per_bc = np.asarray(mdata[modal_key].X.sum(axis=1)).reshape(-1)
    return pd.Series(per_bc, index=mdata[modal_key].obs_names, name=modal_key)


# %%
def arc_counts_per_bc_plot(
    mdata: mudata.MuData,
    ax: matplotlib.axes.Axes,
    ccm: pmbip.CategoricalColormap,
    is_cell_obs_key: str = "cr_is_cell_str",
    **kwargs,
):
    pltdata = pd.concat(
        [
            barcode_totals(mdata, "rna"),
            barcode_totals(mdata, "atac"),
            mdata.obs[is_cell_obs_key],
        ],
        axis=1,
    )
    pmbip.scatter(
        pltdata["atac"],
        pltdata["rna"],
        ax=ax,
        c=ccm.encoded(),
        cmap=ccm.cmap(),
        norm=ccm.norm(),
        s=1.5,
        **kwargs,
    )
    ax.set(
        yscale="log",
        xscale="log",
        xlabel="transpo events per bc",
        ylabel="rna umi per bc",
    )
    scatter_leg = pmbip.Legend.from_mpl_color_dict(ccm.mpl_color_dict())
    scatter_leg.on_ax(ax)


# %%
def adt_ab_total_v_detected_plot(mdata: mudata.MuData, ax: matplotlib.axes.Axes):
    ab_total = np.asarray(mdata["adt"].X.sum(axis=1)).reshape(-1)
    ab_nonzero = np.count_nonzero(mdata["adt"].X.toarray(), axis=1)
    pmbip.scatter(ab_total, ab_nonzero, ax=ax, c="red", s=1.5, alpha=0.1)
    ax.set(xscale="log", xlabel="Total Ab counts", ylabel="Ab detecteded")


# %%
def adt_ab_total_v_control_counts_plot(mdata: mudata.MuData, ax: matplotlib.axes.Axes):
    control_idx = np.where(mdata["adt"].var["IsotypeControl"].to_numpy())[0]
    control_adt_adata = mdata["adt"][:, control_idx]
    ab_total = np.asarray(mdata["adt"].X.sum(axis=1)).reshape(-1)
    control_ab = np.asarray(control_adt_adata.X.sum(axis=1)).reshape(-1)
    pmbip.scatter(ab_total, control_ab, ax=ax, c="red", s=1.5, alpha=0.1)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Total Ab counts",
        ylabel="Total Isotype Ab counts",
    )


# %%
