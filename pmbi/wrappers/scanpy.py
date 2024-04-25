import anndata
import pathlib
import re
import pandas as pd
import scipy.io
import scanpy as sc
import scirpy as ir
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

global MT_CUT
global NFEAT_CUT
global NRNA_CUT

MT_CUT = (0.0, 0.925)
NFEAT_CUT = (0.10, 0.95)
NRNA_CUT = (0.10, 0.95)

class Paneler(object):
    def __init__(self, nrow, ncol, figsize=(3,3), layout="tight", format = "tiff", dpi = 400):
        self.nrow = nrow
        self.ncol = ncol
        self.figsize = figsize
        self.layout = layout
        self.format = format
        self.dpi = dpi
        self.fig, self.axs = plt.subplots(nrows = self.nrow, ncols = self.ncol, figsize = self.figsize, layout = self.layout)
        self.panel_idx = 0
        self.image_num = 1

    def subplot_idx_to_pos(nrow, ncol, idx):
        r = int(np.floor(np.true_divide(idx, ncol))) 
        rem = idx % ncol
        c = rem
        if r >= nrow:
            r = r % nrow
        return (r,c)

    def insert(self, plt_function, **kwargs):
        if self.nrow*self.ncol == 1:
            plt_function(ax=self.ax, **kwargs)
        elif self.nrow or self.ncol:
            plt_function(ax=self.ax[panel_idx], **kwargs)
        else:
            row,col = Paneler.subplot_idx_to_pos(self.nrow, self.ncol, self.panel_idx)
            plt_function(ax=self.ax[row,col], **kwargs)



def read_matrix(matpath, **kwargs):
    matpath = pathlib.Path(matpath)
    if matpath.suffix == ".h5":
        return sc.read_10x_h5(str(matpath), **kwargs)
    elif matpath.suffix == ".h5ad":
        return anndata.read_h5ad(str(matpath), **kwargs)

def obs_names_unique(adata: anndata.AnnData) -> bool:
    return len(adata.obs_names) == len(adata.obs_names.unique())

def var_names_unique(adata: anndata.AnnData) -> bool:
    return len(adata.var_names) == len(adata.var_names.unique())

def combine_adatas(adatas: dict) -> anndata.AnnData:
    idented = {}
    for sample,adata in adatas.items():
        print(sample)
        adata.obs["orig_ident"] = sample
        adata.obs_names = pd.Series(adata.obs_names).str.cat(others=[str(sample)]*len(adata.obs_names), sep="_") 
        idented[sample] = adata 
        if len(adata.obs_names) != len(adata.obs_names.unique()):
            raise
    return anndata.concat(idented.values(), axis = 0, join="outer")

def calc_umi_per_bc(adata: anndata.AnnData) -> pd.DataFrame:
    umi_per_bc = adata.X.sum(axis = 1)[:,0].A1
    umi_per_bc.sort()
    umi_per_bc = pd.DataFrame({
        "rank": pd.Series(range(1,len(umi_per_bc)+1)),
        "umi_count": pd.Series(umi_per_bc[::-1]+1)
                 })
    return umi_per_bc

def ax_umi_per_bc(adata: anndata.AnnData, axs_title = "barcode umi counts"):
    umi_per_bc = calc_umi_per_bc(adata)
    ax = sns.lineplot(umi_per_bc, x="rank", y="umi_count")
    ax.set(yscale="log", ylim=(0,max(umi_per_bc["umi_count"])*2.5), xlim=(0,5e4))
    plt.title(axs_title)
    return ax

def subplot_idx_to_pos(nrow, ncol, idx):
    r = int(np.floor(np.true_divide(idx, ncol))) 
    rem = idx % ncol
    c = rem
    if r >= nrow:
        r = r % nrow
    return (r,c)

def cluster_leiden(adata, resolutions):
    for r in resolutions:
        sc.tl.leiden(adata, resolution=r, key_added = f"leiden_{r}")
        #adata.obs[f"leiden_{r}"] = adata.obs[f"leiden_{r}"].astype(str) + "__" + adata.obs[f"orig_ident"].astype(str)

# implementation of sc.pp.normalize_geometric(adt) from https://github.com/scverse/scanpy/pull/1117
def clr_normalize(adata = anndata.AnnData):
    counts_nda = adata.X.toarray()
    adata.X = scipy.sparse.csr_matrix(np.divide(counts_nda, scipy.stats.mstats.gmean(counts_nda+1, axis=0)))
    sc.pp.log1p(adata)
    return adata

def subset(adata: anndata.AnnData, obs_name: str, bounds: tuple):
    lower,upper = np.quantile(adata.obs[obs_name], bounds)
    return adata[(adata.obs[obs_name] >= lower) & (adata.obs[obs_name] <= upper)]

def std_qc_gex(adata: anndata.AnnData, sample_suffix: str = "sample", plot: bool = True):
    
    print("making names unique")
    adata.var_names_make_unique()

    print("save raw")
    adata.layers["raw_counts"] = adata.X.copy()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    adata = subset(adata, "pct_counts_mt", MT_CUT)
    adata = subset(adata, "n_genes_by_counts", NFEAT_CUT)
    adata = subset(adata, "total_counts", NRNA_CUT)

    if plot:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save = f"_{sample_suffix}")

    return adata

def std_gex(adata: anndata.AnnData, sample_suffix: str = "sample", plot: bool = True):

    adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    #sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)

    adata = adata[:, adata.var.highly_variable]

    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

    sc.pp.scale(adata,  max_value = 10)

    sc.tl.pca(adata, svd_solver="arpack")

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    sc.tl.umap(adata)

    if (plot):
        sc.pl.highly_variable_genes(adata, save = f"_{sample_suffix}")
        sc.pl.pca(adata, save = f"_{sample_suffix}")
        sc.pl.umap(adata, save = f"_{sample_suffix}")

    return adata

def std_adt(adata: anndata.AnnData, sample_suffix: str = "sample", plot: bool = True):

    sc.pp.filter_genes(adata, min_counts=1)

    adata.var["control"] = adata.var_names.str.startswith("Isotype")
    sc.pp.calculate_qc_metrics(
            adata,
            percent_top=(5, 10, 15),
            var_type="antibodies",
            qc_vars=("control",),
            inplace=True,
            )

    adata = clr_normalize(adata)

    sc.pp.scale(adata,  max_value = 10)
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=20)
    sc.tl.leiden(adata, key_added = "adt_leiden")

    sc.tl.umap(adata)

    if plot:

        sns.jointplot(x="log1p_total_counts", y="n_antibodies_by_counts", data=adata.obs, kind="hex", norm=mpl.colors.LogNorm()).fig.savefig(f"adt_qc1_{sample_suffix}.pdf")
        sns.jointplot(x="log1p_total_counts", y="log1p_total_counts_control", data=adata.obs, kind="hex", norm=mpl.colors.LogNorm()).fig.savefig(f"adt_qc2_{sample_suffix}.pdf")
        sc.pl.umap(adata, color = "adt_leiden", size=10, save=f"_{sample_suffix}")

    return adata


if __name__ == "__main__":
    COLOR_GENES = ["CD3D", "CD4", "CD8A", "MS4A1", "CD19", "TRCF", "TCF7"]
    COLOR_AB = ["Hu.CD3_UCHT1", "Hu.CD4_RPA.T4", "Hu.CD8", "Hu.CD20_2H7", "Hu.CD19", "Hu.CD71"]

    adata = sc.read_10x_h5("/home/ubuntu/mnt/betts/aggr/135/outs/count/filtered_feature_bc_matrix.h5", gex_only = False)
    adata.layers["counts"] = adata.X.copy()
    adata.var_names_make_unique()

    adata.var.feature_types.value_counts()
    gex = std_gex(adata[:,adata.var.feature_types == "Gene Expression"], plot=False)
    adt = std_adt(adata[:,adata.var.feature_types == "Antibody Capture"], plot = False)

    gex.x.std()
    adt.x.std()

    [x for x in gex.var_names if x == "cd19"]

    adt
    gex
    [x for x in list(adt.var_names)[0:130] if x.split('.')[1] in color_genes]

    import itertools


    comb = anndata.concat([gex, adt], axis = 1, merge = "first")
    comb
    sc.pl.umap(comb, color=[g for g in list(itertools.chain(*[x for x in zip(color_genes, color_ab)])) if g in comb.var_names], save = true, ncols = 2)

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
    p13_raw.obs["orig_ident"] = p13_raw.obs["orig_ident"].apply(lambda x: re.sub("[atcg]+[-][0-9]+[_]", "", x))
    p13_raw.obs["old_index"] = p13_raw.obs.index
    p13_raw.var["feature"] = p13_raw.var.index
    p13_raw.obs = p13_raw.obs.reset_index()

    p13_raw.var.shape
    for oi in p13_raw.obs["orig_ident"].unique():
        which = p13_raw.obs.loc[p13_raw.obs["orig_ident"] == oi].index
        which_adata = p13_raw[which, :]
        scipy.io.mmwrite(f"p13_split/{oi}.mm", which_adata.x.transpose())
        which_adata.obs.to_csv(f"p13_split/{oi}.obs", sep = "\t", header = false, index = false)
        which_adata.var.to_csv(f"p13_split/{oi}.var", sep="\t", header = false, index = false)
