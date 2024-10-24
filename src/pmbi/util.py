from pprint import pprint

import awkward as ak
import pandas as pd
import numpy as np


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


def unpack_nested_dict(d, keys=[]):
    ele = []
    for k, v in d.items():
        if isinstance(v, dict):
            ele.extend(unpack_nested_dict(v, keys + [k]))
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
