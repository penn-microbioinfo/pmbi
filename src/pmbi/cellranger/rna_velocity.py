# %%
import os
from pathlib import Path
from typing import Union

import anndata
import loompy
import scipy.io
import pandas as pd
import numpy as np
import scanpy as sc

import pmbi.anndata.io as aio
import pmbi.plotting as pmbip

# %%
def read_star_solo_mtx(
    matrix_dir: os.PathLike,
    matrix_fn: str = "matrix.mtx",
    barcodes_fn: str = "barcodes.tsv",
    features_fn: str = "features.tsv",
    barcodes_as_rows=True,
) -> anndata.AnnData:
    matrix_dir = Path(matrix_dir)
    adata = aio.read_mtx(
        mtx_path=matrix_dir.joinpath(matrix_fn),
        var_names_path=matrix_dir.joinpath(barcodes_fn),
        obs_names_path=matrix_dir.joinpath(features_fn),
    )
    if barcodes_as_rows:
        adata = adata.transpose()
        adata.X = adata.X.tocsr()
        return adata
    else:
        return adata


# %%
def add_star_solo_velocyto_layers(
    adata: anndata.AnnData,
    velocyto_outdir: os.PathLike,
    barcodes_as_rows=True,
    inplace=True,
) -> Union[anndata.AnnData, None]:
    if not inplace:
        adata = adata.copy()
    velocyto_outdir = Path(velocyto_outdir)
    expected_fn = ["spliced.mtx", "unspliced.mtx", "ambiguous.mtx"]
    for fn in expected_fn:
        layer = fn.replace(".mtx", "")
        mat = scipy.io.mmread(velocyto_outdir.joinpath(fn)).tocsr()
        if barcodes_as_rows:
            mat = mat.T.tocsr()
        adata.layers[layer] = mat
    if not inplace:
        return adata


# %% NOTE: Probably unnecessay... just use anndata.AnnData.write_loom() instead
# def adata_write_loom(adata: anndata.AnnData, outfile: os.PathLike):
#     layer_dict = {"": adata.X}
#     layer_dict.update(dict(adata.layers))
#     loompy.create(
#         filename=outfile,
#         layers=layer_dict,
#         row_attrs={"barcode": adata.obs.index.to_numpy()},
#         col_attrs={"feature": adata.var.index.to_numpy()},
#     )

if __name__ == "__main__":
    pass
# %% Read from Seurat
seu_adata = aio.read_mtx(
        mtx_path = "/home/amsesk/super1/montse/object_parts/counts.mtx",
        obs_names_path = "/home/amsesk/super1/montse/object_parts/features.tsv",
        var_names_path = "/home/amsesk/super1/montse/object_parts/barcodes.tsv"
).transpose()
obs = pd.read_csv("/home/amsesk/super1/montse/object_parts/barcode_metadata.tsv", sep="\t")
var = pd.read_csv("/home/amsesk/super1/montse/object_parts/feature_metadata.tsv", sep="\t")


# %%
seu_adata.obs = pd.concat([seu_adata.obs, obs.set_index("barcode", drop=True)], axis=1)
seu_adata.var = pd.concat([seu_adata.var, var.set_index("feature", drop=True)], axis=1)

seu_adata.obsm["umap"] = pd.read_csv("/home/amsesk/super1/montse/object_parts/umap.tsv", sep="\t").set_index("barcode", drop=True).to_numpy()

# mtx = scipy.io.mmread("/home/amsesk/super1/montse/object_parts/counts.mtx").tocsr()
# seu_adata.layers["counts"]=mtx.T.tocsr()
# seu_adata.layers["scale.data"] = pd.read_csv("/home/amsesk/super1/montse/object_parts/scale.data.tsv", sep="\t").set_index("feature", drop=True).transpose()

# %% Read from STARsolo 
solo_adatas = {}
wells = ["well_1", "well_2", "well_3", "well_4", "well_5", "well_6"]
for well in wells:
    adata = read_star_solo_mtx(
        f"/home/amsesk/super1/montse/STAR_runs/{well}Solo.out/Gene/filtered",
        features_fn="feature_names.tsv"
    )
    add_star_solo_velocyto_layers(
        adata, f"/home/amsesk/super1/montse/STAR_runs/{well}Solo.out/Velocyto/filtered"
    )
    adata.obs_names = list(map(lambda n: f"{n}-1", adata.obs_names))
    adata = adata[pd.Series(adata.obs_names).isin(seu_adata.obs_names),:]
    adata = adata[:,pd.Series(adata.var_names).isin(seu_adata.var_names)]
    adata = adata[:,pd.Series(adata.var_names).drop_duplicates().index]
    solo_adatas[well]=adata

# %%

solo_concat=anndata.concat(solo_adatas, axis=0)
solo_concat = solo_concat[pd.Series(solo_concat.obs_names).drop_duplicates().index,:]
(solo_concat.var_names.value_counts()>1).any()
(solo_concat.obs_names.value_counts()>1).any()
solo_concat.shape

seu_adata=seu_adata[solo_concat.obs_names,solo_concat.var_names]
seu_adata.obs["clusters"] = pd.Series(seu_adata.obs["seurat_clusters"], dtype="category")

# %%
seu_adata.layers["spliced"] = solo_concat.layers["spliced"].copy()
seu_adata.layers["unspliced"] = solo_concat.layers["unspliced"].copy()
seu_adata.layers["spliced"]

scv.pl.proportions(seu_adata, save="/home/amsesk/figures/montse/velo/prop.png")
# %%
# scv.pp.filter_genes(seu_adata, min_shared_counts=20)
# scv.pp.normalize_per_cell(seu_adata)
# scv.pp.filter_genes_dispersion(seu_adata, n_top_genes=2000)
# sc.pp.log1p(seu_adata)
scv.pp.filter_and_normalize(seu_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(seu_adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(seu_adata)
scv.tl.velocity_graph(seu_adata)

scv.pl.velocity_embedding_stream(seu_adata, basis='umap', save="/home/amsesk/figures/montse/velo/stream.png")


seu_adata

panel = pmbip.MosaicPaneler()
panel._figure.savefig("/home/amsesk/figures/montse/velo/umap.png")

