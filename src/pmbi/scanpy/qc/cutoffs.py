from pathlib import Path
from os import PathLike
import pandas as pd
import toolz

# %%
def _is_cutoff_limit(column_name, min_suff="min", max_suff="max"):
    """
    Determine if a string (column_name) matches the expected format for a QC cutoff limit.
    That format is either: "XXXXX_<min_suff>" or "XXXXX_<max_suff>"

    Allows filtering of DataFrame column names based on whether or not they match this format.
    Assists with filtering input into _cutoff_names, which does no such check.
    """
    if column_name.endswith(f"_{min_suff}") or column_name.endswith(f"_{max_suff}"):
        return True
    else:
        return False


# %%
def _cutoff_names(columns: pd.Index, min_suff="min", max_suff="max"):
    """
    Condenses column names of a cutoff limits to only include the actual name of the cutoff.
    """
    return set(
        [
            x.replace(f"_{min_suff}", "").replace(f"_{max_suff}", "")
            for x in columns.to_list()
        ]
    )


# %%
def read_cutoff_sheet(
    path, sample_id_column="sample_id", min_suff="min", max_suff="max"
):
    """
    Read in a sheet from disk derived from cutoff generation as part of pmbi.scanpy
    Enables saving of cutoffs in sheets to avoid redundantly having to rerun the long-running
    steps that generate such cutoffs.

    Pairs with write_cutoff_sheet.
    """
    df = pd.read_csv(path)
    cutoff_limit_cols = toolz.groupby(_is_cutoff_limit, df.columns.to_list())[True]
    cutoff_names = _cutoff_names(pd.Index(cutoff_limit_cols))
    tupled_cols = []
    for c in cutoff_names:
        tupled_cols.append(
            pd.Series(list(zip(df[f"{c}_{min_suff}"], df[f"{c}_{max_suff}"])), name=c)
        )
    new_df = pd.concat(tupled_cols, axis=1)
    new_df = new_df.set_index(df[sample_id_column])
    return new_df


# %%
def cutoffs_as_sheet(df: pd.DataFrame, cutoff_names: list):
    """
    Convert cutoff dataframe to writable format and write.

    Input format represents cutoffs as series of series with limits and cell counts before/after (it's really messy).

    Output format extracts the limits, splits them into 2 columns (ie min/max) for each cutoff name and writes to disk
    along with the Index reset to the leading column.

    Pairs with write_cutoff_sheet.
    """
    df = df[cutoff_names].apply(lambda x: x.apply(lambda xx: xx["limits"]))
    for c in cutoff_names:
        df[[f"{c}_min", f"{c}_max"]] = df[c].apply(pd.Series)
    df = df.drop(columns=cutoff_names).reset_index()
    return df

# %%
def apply_cell_cutoff_to_adata(adata, qckey, min_value, max_value):
    adata_filt = adata[
        (adata.obs[qckey] >= min_value)
        & (adata.obs[qckey] <= max_value)
    ]
    adata_filt.uns[f"cutoffs__{qckey}"] = [min_value, max_value]
    return adata_filt

