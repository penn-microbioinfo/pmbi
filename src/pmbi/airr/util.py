import numpy as np
import pandas as pd
import scirpy as ir
import scirpy.util
import pmbi.util
import awkward as ak

from pmbi.anndata.io import pickle_piece
from anndata import AnnData

# %%
def n_overlapping_junction_aa_identical(df1, df2, subject_key="donor"):
    matrix = []
    g1_subjects = sorted(df1[subject_key].unique())
    g2_subjects = sorted(df2[subject_key].unique())
    for g1_member in g1_subjects:
        row = []
        g1_sub = df1[df1[subject_key] == g1_member].junction_aa.to_numpy()
        for g2_member in g2_subjects:
            g2_sub = df2[df2[subject_key] == g2_member].junction_aa.to_numpy()
            n_overlaps = np.intersect1d(g1_sub, g2_sub).shape
            row.append(n_overlaps[0])
        matrix.append(row)
    return pd.DataFrame(np.array(matrix), index=g1_subjects, columns=g2_subjects)


# %%
def n_unique_junction_aa(
    df1, df2, subject_key="donor", count_key="junction_aa", axis=0
):
    matrix = []
    index = sorted(df1[subject_key].unique())
    columns = sorted(df2[subject_key].unique())
    matrix = pd.DataFrame(0, index=index, columns=columns)
    if axis == 0:
        for g1_member in index:
            g1_total = df1[df1[subject_key] == g1_member][count_key].unique().shape[0]
            matrix.loc[g1_member, :] = g1_total
    elif axis == 1:
        for g2_member in index:
            g2_total = df2[df2[subject_key] == g2_member][count_key].unique().shape[0]
            matrix.loc[g2_member, :] = g2_total
    else:
        raise pd.errors.IndexingError
    return matrix


# %%
def pull_airr_field_conditional(data, conditional, airr_mod="airr", airr_key="airr"):
    airr = ir.util.DataHandler(data, airr_mod=airr_mod, airr_key=airr_key).airr
    return airr[conditional(airr)]


# %%
def write_ir_dist(data, uns_key, output_dir, airr_mod="airr", airr_key="airr"):
    adata = ir.util.DataHandler(data, airr_mod=airr_mod, airr_key=airr_key).adata
    pickle_piece(adata=adata, outer_key="uns", inner_key=uns_key, output_dir=output_dir)


# %%
def clonotypes_to_obs(data, uns_key, key_added, inner_key="cell_indices", airr_mod = "airr", airr_key="airr", inplace = True):
    dh = scirpy.util.DataHandler(data, airr_mod=airr_mod, airr_key=airr_key)
    cell_clonotypes = pmbi.util.invert_dict(dh.adata.uns[uns_key][inner_key])
    cell_clonotypes.name = key_added
    dh.adata.obs = pd.merge(dh.adata.obs, cell_clonotypes, how = "left", left_index = True, right_index=True)
    if inplace:
        return None
    else:
        return dh.data.copy()

# %%
def clonotype_sizes(data, obs_key, key_added, airr_mod = "airr", airr_key="airr", inplace = True): 
    dh = scirpy.util.DataHandler(data, airr_mod=airr_mod, airr_key=airr_key)
    clonotype_size=pd.Series(dh.adata.obs.groupby(obs_key).size(), name = key_added)
    dh.adata.obs = pd.merge(dh.adata.obs, clonotype_size, how="left", left_on=obs_key, right_index=True)
    if inplace:
        return None
    else:
        return dh.data.copy()

# %%
def airr_locus_df(adata: AnnData, locus: str, columns: list[str]|None = None):
    if columns is None:
        return ak.to_dataframe(pull_airr_field_conditional(adata, lambda x: x.locus == locus)).reset_index()
    else:
        return ak.to_dataframe(pull_airr_field_conditional(adata, lambda x: x.locus == locus)[columns]).reset_index()

# %%
