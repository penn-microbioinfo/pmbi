import os

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import palettable
import pmbi.plotting as pmbip
import scanpy as sc
# import scanpy.plotting._ridgeplot as spr


def expr_summary_figure(adata, gene, outdir=".", figsize=(6, 8), format="png", layer="scale_data", umap_key="X_umap", clust_key="leiden_1.0"):
    panel = pmbip.Paneler(1, 2, figsize=figsize, format=format, width_ratios=[2, 3])
    gridspec = panel.axs[0].get_subplotspec().get_gridspec()
    panel.axs[0].remove()
    gridspec[0]
    subfig = panel.fig.add_subfigure(gridspec[0])
    axsLeft = subfig.subplots(3, 1, sharex=True)
    colormap = plt.get_cmap("gist_ncar")
    nclust = len(adata.obs[clust_key].unique())
    colors = colormap(np.linspace(0, 1, nclust))
    sc.pl.embedding(
        adata,
        basis=umap_key,
        color=clust_key,
        ax=axsLeft[0],
        show=False,
        layer=layer,
        legend_loc="on data",
        title="",
        # palette = colors,
        legend_fontsize=5,
    )
    # clust_col = np.array(
    #     [matplotlib.colors.to_rgba(c) for c in adata.uns["leiden_1.0_colors"]]
    # )
    gene_idx = np.where([x == gene for x in adata.var_names])[0]
    top_1_perc = np.quantile(adata.layers[layer][:, gene_idx], 0.99)
    sc.pl.embedding(
        adata,
        basis=umap_key,
        color=gene,
        ax=axsLeft[1],
        show=False,
        layer=layer,
        title="",
        vmax=top_1_perc,
        cmap=palettable.colorbrewer.diverging.RdYlBu_11_r.mpl_colormap,
        legend_fontsize=5,
    )
    # spr.ridgeplot(
    #     adata,
    #     gene,
    #     "vaeda_calls",
    #     bandwidth=0.01,
    #     save=False,
    #     return_fig=False,
    #     ax=axsLeft[2],
    #     title="",
    #     figsize=figsize,
    # )
    # panel.fig.get_axes()[5].set_yticks(ticks=[], labels=[])
    # panel.fig.get_axes()[5].set_ylabel("")

    # spr.ridgeplot(
    #     adata,
    #     gene,
    #     "leiden_1.0",
    #     bandwidth=0.1,
    #     save=False,
    #     return_fig=False,
    #     # cmap=palettable.colorbrewer.diverging.RdYlBu_11_r.mpl_colormap,
    #     ax=panel.axs[1],
    #     figsize=figsize,
    #     palette=clust_col,
    #     title="",
    # )

    panel.fig.suptitle(gene, y = 0.99, fontsize = 12)
    panel.fig.savefig(os.path.join(outdir, f"{gene}_summary.{format}"))
