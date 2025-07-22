from pprint import pprint
from pathlib import Path
import os
import re
from traceback import format_list

import awkward as ak
import pandas as pd
import numpy as np

# %% CHUNK: This function wraps re.search to pull a substring from a string - substring should be under capture group 1 {{{
def get_substring(string: str, pattern: str) -> str:
    s = re.search(pattern, os.path.basename(string))
    if s is not None:
        if len(s.groups()) > 1:
            raise ValueError(
                f"More than one matching substring for pattern `{pattern}` in string: {string}"
            )
        else:
            return s.group(1)
    else:
        raise ValueError(
            f"No match found for pattern `{pattern}` in string: `{string}`"
        )
# }}}

def ak_print(arr: ak.Array) -> None:
    pprint(ak.to_list(arr))


def invert_dict(d, return_type="series") -> pd.Series:
    s = pd.DataFrame(
        {"clonotype_id": pd.Series(d.keys()), "barcode": pd.Series(d.values())}
    )
    if return_type == "series":
        return s.explode("barcode").set_index("barcode").loc[:, "clonotype_id"]
    elif return_type == "dict":
        raise NotImplementedError(f"Unimplemented return_type: {return_type}")
    else:
        raise NotImplementedError(f"Unimplemented return_type: {return_type}")


def unpack_nested_dict(d, keys=[], depth:int = 0, stop_at:int|None = None):
    ele = []
    for k, v in d.items():
        if depth is not None and depth<=stop_at:
            if isinstance(v, dict):
                depth += 1
                ele.extend(unpack_nested_dict(v, keys + [k], depth, stop_at))
            else:
                ele.append(tuple(keys + [k, v]))
        else:
            ele.append(tuple(keys + [k, v]))
    return tuple(ele)


def n_intersect_table(
    df: pd.DataFrame, group_by: str, intersect_by: str
) -> pd.DataFrame:
    """
    Computes a table showing the number of intersections between groups in a DataFrame.

    This function takes a DataFrame and two column names, grouping the DataFrame by the `group_by` column,
    and then calculates the number of common occurrences in the `intersect_by` column for each group.
    The result is returned as a DataFrame where rows represent groups and columns represent the count of intersections.

    Parameters:
    df (pd.DataFrame): The DataFrame.
    group_by (str): The name of the column to group the data by.
    intersect_by (str): The name of the column to find intersections within each group.

    Returns:
    pd.DataFrame: A DataFrame with the groups as rows/columns and the intersection counts in the cells.
    """
    groups = df.groupby(group_by).groups.keys()
    matrix = []
    for g1 in groups:
        g1_values = df[df[group_by] == g1][intersect_by].unique()
        row = []
        for g2 in groups:
            g2_values = df[df[group_by] == g2][intersect_by].unique()
            n_overlap = np.intersect1d(g1_values, g2_values).shape[0]
            row.append(n_overlap)
        matrix.append(row)
    index = pd.Index(groups)
    matrix = pd.DataFrame(matrix, index=index, columns=index)
    return matrix

# %% FUNC: Check if we are in an interactive python session or not
def is_ipython():
    try:
        _ipy = get_ipython()
        return True
    except NameError:
        return False

# %% FUNC: Check for broken symlinks
def symlink_is_valid(path: os.PathLike) -> bool:
    path = Path(path)
    if path.exists(follow_symlinks=False):
        if path.is_symlink():
            return path.exists(follow_symlinks=True)
        else:
            raise ValueError(f"Path is not a symlink: {path}")
    else:
        raise ValueError(f"Path does not exist: {path}")


# %%
symlink_is_valid("/home/amsesk/super1/t1d-coculture/all_fastq_symlinks/HPAP-135_CC_1_ADT_S49_L001_R1_001.fastq.gz")
symlink_is_valid("/home/amsesk/super1/t1d-coculture/blahblah")

Path(Path(Path("/home/am")))


# %%

