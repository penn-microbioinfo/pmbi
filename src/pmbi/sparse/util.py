import Levenshtein as Lev
import itertools
import pandas as pd
import numpy as np
import muon
import awkward as ak

# %%
def ir_cell_indices_to_pd(cellind_dict):
    cell_indices = pd.Series(cellind_dict.values(), index = cellind_dict.keys())
    cell_indices = cell_indices.apply(lambda r: r.astype(str))
    return cell_indices

# %%
def ir_airr_juncAA_to_pd(airr, cell_id_col = "barcode", locus_col = "locus", seq_col = "junction_aa"):
    seqs = ak.to_dataframe(airr[[cell_id_col, locus_col, seq_col]])
    seqs = seqs.reset_index()[[cell_id_col, locus_col, seq_col]]
    seqs.columns = ["cell_id", locus_col, seq_col]
    seqs["cell_id"] = seqs["cell_id"].apply(lambda r: r.astype(str))
    return seqs

# %%
query = mdata["airr"].uns["cc_aa_vdjdb_lev_cutoff2_query"]

# %% CHUNK: Make pandas DataFrame from airr awkward array with locus names and junction_aa seqs
mdata["airr"].obsm["airr"]["barcode"] = mdata["airr"].obs.index
our_seqs = ir_airr_juncAA_to_pd(mdata["airr"].obsm["airr"])
vdjdb_seqs = ir_airr_juncAA_to_pd(vdjdb.obsm["airr"], cell_id_col = "cell_id")

# %%
our_seqs
vdjdb_seqs

# %%
ci_rows = ir_cell_indices_to_pd(query["cell_indices"])
ci_cols = ir_cell_indices_to_pd(query["cell_indices_reference"])

# %%
s = query["distances"]
nz = s.nonzero()

levs = np.empty(nz[0].shape)
x = 0
for i,j in zip(*s.nonzero()):
    seqs1 = our_seqs[our_seqs["cell_id"].isin(ci_rows.iloc[i])]["junction_aa"].to_numpy()
    seqs2 = vdjdb_seqs[vdjdb_seqs["cell_id"].isin(ci_cols.iloc[j])]["junction_aa"].to_numpy()
    dists = list()
    for c in np.array(np.meshgrid(seqs1, seqs2)).T.reshape(-1,2):
        d = Lev.distance(c[0], c[1])
        dists.append(d)
    levs[x] = min(dists)
    x += 1

levs
print(min(levs), max(levs), len(levs))

# %% 
s[0,249]

# %%
seqs1 = seqs1["junction_aa"].to_numpy()
seqs2 = seqs2["junction_aa"].to_numpy()

# %%
l = LevenshteinDistanceCalculator(cutoff = 100)
dmat = l.calc_dist_mat(seqs1["junction_aa"].to_numpy(), seqs2["junction_aa"].to_numpy()).toarray()
dmat

# %%


# %%
# %%
l.
# ss = s[rows,:]
ss = ss[:,cols]

s_arr = s.toarray()
s_arr[0,249]
ss_arr = ss.toarray()
print(s_arr.shape, ss_arr.shape)

s.indices

cols
np.where(s_arr[:,cols] != ss_arr)

