from typing import Optional, Union

import matplotlib.axes
import mudata
import muon.atac as atac
import muon as mu
import numpy as np
import pandas as pd
import pysam
import snapatac2.genome
from sklearn.neighbors import KernelDensity
import snapatac2

import pmbi.plotting as pmbip
from pmbi.wrappers.scanpy import calc_umi_per_bc


# %%
def tss_enrichment(
    mdata: mudata.MuData,
    n_tss: int = 2000,
    extend_upstream: int = 2000,
    extend_downstream: int = 2000,
    features: Union[pd.DataFrame, None] = None,
    **kwargs
):
    """
    Calculates enrichment scores for positions up and downstream of gene TSSs.

    Args:
        mdata (muon.Mdata):
        n_tss (int):
        extend_upstream (int):
        extend_downstream (int):
        features (pd.DataFrame): DataFrame with (atleast) the columns Chromosome, Start, and End; 
            feature ids (e.g., `gene_id` from GTF) should be in the index

    Raises:
        ValueError: For a variety of argument checks
        KeyError: If path to the fragments file is missing from MuData["atac"].uns["files"]

    Returns:
        anndata.AnnData: contains enrichment scores

    Example:

    """
    if "atac" not in mdata.mod_names:
        raise ValueError("expected MuData with 'atac' modality")

    if "rna" in mdata.mod_names:
        features = atac.tl.get_gene_annotation_from_rna(mdata)
    else:
        if features is None:
            raise ValueError(
                "`features` must be provided if there is no 'rna' modalitiy in MuData"
            )
        else:
            expected_cols = ["Chromosome", "Start", "End"]
            if not all([c in features.columns for c in expected_cols]):
                raise ValueError(
                    f"features dataframe must have columns: {expected_cols}"
                )

    try:
        fragments = pysam.TabixFile(
            mdata["atac"].uns["files"]["fragments"], parser=pysam.asBed()
        )
    except KeyError:
        raise KeyError(
            "expected path to fragments file in AnnData.uns['files']['fragments']"
        )

    features_filt = features[features.Chromosome.isin(fragments.contigs)]
    features_with_window = []
    for row in features_filt.itertuples():
        try:
            _ = fragments.fetch(
                row.Chromosome, row.Start - extend_upstream, row.End + extend_downstream
            )
            features_with_window.append(row.Index)
        except ValueError:
            pass
    features_filt = features_filt.loc[features_with_window, :]
    return atac.tl.tss_enrichment(
        mdata,
        features=features_filt,
        n_tss=n_tss,
        extend_upstream=extend_upstream,
        extend_downstream=extend_downstream,
        **kwargs,
    )


def tss_enrichment_plot(
    mdata: mudata.MuData,
    ax: matplotlib.axes.Axes,
    n_tss: int = 2000,
    extend_upstream: int = 2000,
    extend_downstream: int = 2000,
    features: Union[pd.DataFrame, None] = None,
    color="red",
    fill_alpha=0.2,
    **kwargs
):
    tsse = tss_enrichment(
        mdata=mdata,
        n_tss=n_tss,
        extend_upstream=extend_upstream,
        extend_downstream=extend_downstream,
        features=features,
        **kwargs,
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


def peak_counts_scatter_with_dist(mdata, axs):
    ic_pal = {"Is cell": "#118dff", "Not cell": "#ae12a7"}
    adata = mdata["atac_raw"]
    fop_per_bc = pd.Series(np.asarray(adata.X.sum(axis=1)).reshape(-1)).to_numpy()
    fop_per_bc_plt = (
        pd.DataFrame(
            {
                "fop_per_bc": fop_per_bc,
                "cr_is_cell_str": adata.obs["cr_is_cell_str"].to_numpy(),
            },
        )
        .sort_values("fop_per_bc", ascending=False)
        .pipe(lambda df: df.assign(sort_idx=np.arange(0, df.shape[0])))
    )
    ccm = pmbip.CategoricalColormap(
        values=fop_per_bc_plt["cr_is_cell_str"], colors=ic_pal
    )
    nsteps = 100
    xmax = adata.shape[0]
    xmax_log = np.log10(xmax)
    idx_step = xmax_log / nsteps
    si_plt = {}
    for c in ["Is cell", "Not cell"]:
        c_si = fop_per_bc_plt.loc[
            fop_per_bc_plt["cr_is_cell_str"] == c, "sort_idx"
        ].to_numpy()
        c_si_log = np.log10(c_si + 1)
        kde = KernelDensity(kernel="gaussian", bandwidth=idx_step).fit(
            c_si_log[:, np.newaxis]
        )
        x = np.arange(0, xmax_log, step=idx_step * (1e-1))[:, np.newaxis]
        si_plt[c] = {"x": (10 ** (x)) - 1, "y": np.exp(kde.score_samples(x))}
    axs[0].plot(
        si_plt["Is cell"]["x"],
        si_plt["Is cell"]["y"],
        color=ic_pal["Is cell"],
        linewidth=0.8,
    )
    axs[0].plot(
        si_plt["Not cell"]["x"],
        si_plt["Not cell"]["y"],
        color=ic_pal["Not cell"],
        linewidth=0.8,
    )
    axs[0].set(
        xlim=[1e0, xmax],
        xscale="log",
        ylabel="density",
    )
    axs[1].scatter(
        fop_per_bc_plt["sort_idx"],
        fop_per_bc_plt["fop_per_bc"],
        c=ccm.encoded(),
        cmap=ccm.cmap(),
        norm=ccm.norm(),
        s=1.5,
        alpha=0.5,
    )
    axs[1].set(
        ylim=[1e0, fop_per_bc_plt["fop_per_bc"].max()],
        xlim=[1e0, xmax],
        xscale="log",
        yscale="log",
        xlabel="barcodes",
        ylabel="Total peak counts",
    )
    scatter_leg = pmbip.Legend.from_mpl_color_dict(ccm.mpl_color_dict())
    scatter_leg.on_ax(axs[1])

# %%

def asap_qc_panel(mdata: mu.MuData, genome: snapatac2.genome.Genome, features: pd.DataFrame):
    design = [
        ["bc_counts_distr", "bc_counts_distr"],
        ["bc_counts", "bc_counts"],
        ["bc_counts", "bc_counts"],
        ["fragment_size", "tss_enrichment"],
        ["adt1", "adt2"]
    ]
    mosaic=pmbip.MosaicPaneler(design=design, figsize=(5,8))
    peak_counts_scatter_with_dist(mdata=mdata, axs=list(mosaic.get_axs(keys=["bc_counts_distr", "bc_counts"]).values()))
    fragment_size_distr_plot(mdata=mdata, genome=genome, ax=mosaic.get_ax("fragment_size"), linewidth=0.5, c="green")
    tss_enrichment_plot(mdata=mdata, ax=mosaic.get_ax("tss_enrichment"), color="green", features=features, barcodes="original_barcode")
    adt_ab_total_v_detected_plot(mdata, mosaic.get_ax("adt1"))
    adt_ab_total_v_control_counts_plot(mdata, mosaic.get_ax("adt2"))
    return mosaic

