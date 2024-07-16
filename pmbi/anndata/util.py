# %%
import functools
import re

import anndata
import pandas as pd

BARCODE_PATTERN = re.compile("([ATCG]+)([-][0-9])([-][0-9])*")


# %%
def obs_names_unique(adata: anndata.AnnData) -> bool:
    """
    Check if obs_names in the AnnData object are unique.

    Parameters:
    -----------
    adata : anndata.AnnData
        The AnnData object to check for unique obs_names.

    Returns:
    --------
    bool
        True if obs_names are unique, False otherwise.
    """
    return len(adata.obs_names) == len(adata.obs_names.unique())


# %%
def var_names_unique(adata: anndata.AnnData) -> bool:
    """
    Check if var_names in the AnnData object are unique.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to check for unique var_names.

    Returns
    -------
    bool
        True if the var_names are unique, False otherwise.
    """
    return len(adata.var_names) == len(adata.var_names.unique())


# %%
def _original_barcode(barcode: str, pattern: re.Pattern = BARCODE_PATTERN) -> str:
    """
    Extracts the original barcode from the given string based on a specific pattern.

    Args:
        barcode (str): The string containing the barcode to be extracted.
        pattern (re.Pattern, optional): The regex pattern used to identify the barcode. Defaults to "([ATCG]+)([-][0-9])([-][0-9])*", defined in pmbi.anndata.util.

    Returns:
        str: The extracted barcode if the pattern matches, otherwise returns an empty string.
    """
    m = re.match(pattern, barcode)
    if m is not None:
        return f"{m.groups()[0]}{m.groups()[1]}"
    else:
        raise ValueError(f"Barcode does not match pattern: {pattern}")


# %%
def original_barcodes(current_barcodes: pd.Series | pd.Index) -> pd.Series:
    """
    Processes the current barcodes and returns the original barcodes.

    Args:
        current_barcodes (pd.Series|pd.Index): A pandas Series/Index object containing the current barcodes.

    Returns:
        pd.Series|pd.Index: A pandas Series object containing the original barcodes.
    """
    if isinstance(current_barcodes, pd.Series):
        return current_barcodes.apply(_original_barcode)
    elif isinstance(current_barcodes, pd.Index):
        return pd.Series(current_barcodes.to_series().apply(_original_barcode))
    else:
        raise ValueError(f"Expected pd.Series, got: {type(original_barcodes)}")


# %%
def shared_barcodes(adatas: list) -> pd.Series:
    """
    Identify shared barcodes among a list of AnnData objects.

    Parameters:
    adatas (list): A list of AnnData objects from which to identify shared barcodes.

    Returns:
    pd.Series: A pandas Series containing the barcodes that are shared among all input AnnData objects.
    """
    if any([not obs_names_unique(adata) for adata in adatas]):
        raise ValueError(
            "Member of adatas given to shared_barcodes has non-unique barcodes. Fix this and canonicalize between members of list before calling again."
        )
    indices = [pd.Series(adata.obs.index) for adata in adatas]
    merged = functools.reduce(lambda x, y: x[x.isin(y)], indices)
    return merged


# %%
def get_barcode_mapper(adata: anndata.AnnData, batch_key: str) -> pd.DataFrame:
    """
    Generate a DataFrame mapping original barcodes to unique barcodes and an batch_key that
    specifies what batch the observations come from. This can be used to canonicalize unique
    and common barcodes.

    Parameters
    ----------
    adata : anndata.AnnData
        An AnnData object upon which make_obs_unique has already been run.
    batch_key : str
        The key in `adata.obs` that corresponds to batch annotations.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing original barcodes, unique barcodes, and batch annotations as columns.
    """
    if not obs_names_unique(adata):
        logging.warning("Given adata has non-unique barcodes. 'unique barcodes' in returned mapper are not unique")
        # raise ValueError(
        #     "Given adata has non-unique barcodes. Run anndata.make_obs_names_unique"
        # )
    return pd.DataFrame(
        {
            k: v.to_numpy()
            for k, v in {
                "unique_barcode": adata.obs.index,
                "original_barcode": original_barcodes(adata.obs.index),
                batch_key: adata.obs[batch_key],
            }.items()
        }
    )


# %%
def obs_canonicalize_barcodes(
    adata: anndata.AnnData, based_on: anndata.AnnData, batch_key: str
) -> anndata.AnnData:
    """
    Canonicalizes the barcodes in the observation (obs) of an AnnData object based on the unique obs_names defined
    in another object AnnData object. These unique names often come from running make_obs_unique on an object.

    Barcodes unique in `based_on` are assigned to non-unique barcodes in `adata` by...
    1) deriving the original barcode with `original_barcode` from this module
    2) constucting a DataFrame containing with columns: `unique_barcode`, `original_barcode`, and `batch_key` from 
    `based_on`
    3) merging that DataFrame with one constructed from `adata.obs.index` and `batch_key` from `adata`

    Values in `unique_barcode` in the merge product represent the canonicalized barcodes since original barcodes 
    should not be duplicated within batches.

    Barcodes present in the index of `adata` with no corresponding index in `based_on` will be `nan` in the merge product.
    First, those are replaced back with their original barcode from `adata`.
    Finally, the entire index is deduplicated with anndata.utils.`make_index_unique`. NOTE: Because of this final step,
    barcodes are likely to change between `based_on` and the Index of the returned AnnData.

    Parameters:
    adata (anndata.AnnData): The AnnData object whose observation barcodes need to be canonicalized with another.
    based_on (anndata.AnnData): The AnnData object that serves as the reference for canonicalizing the barcodes.
    batch_key (str): The key used to identify batches in the observations of the AnnData objects.

    Returns:
    anndata.AnnData: The AnnData object with canonicalized observation barcodes.
    """
    adata = adata.copy()
    adata.obs = adata.obs.reset_index(names="original_barcode")
    mapper = get_barcode_mapper(adata=based_on, batch_key=batch_key)
    mapper = pd.merge(
        mapper,
        adata.obs[["original_barcode", batch_key]],
        how="right",
        on=["original_barcode", batch_key],
    )
    if any(
        [
            len(v) > 1
            for k, v in mapper.groupby(["original_barcode", batch_key]).groups.items()
        ]
    ):
        raise ValueError(
            "There are common original barcodes within batches. Are you sure you selected the correct batch_key?"
        )
    mapper["unique_barcode"] = mapper.apply(
        lambda row: (
            row["original_barcode"]
            if pd.isnull(row["unique_barcode"])
            else row["unique_barcode"]
        ),
        axis=1,
    )
    new_index = make_index_unique(pd.Index(mapper["unique_barcode"]), join = "-")
    adata.obs = adata.obs.set_index(new_index)
    assert obs_names_unique(adata), "Obs names are not unique. Something has gone wrong."
    return adata


# %%
