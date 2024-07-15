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
        raise ValueError(
            "Given adata has non-unique barcodes. Run anndata.make_obs_names_unique"
        )
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

    Parameters:
    adata (anndata.AnnData): The AnnData object whose observation barcodes need to be canonicalized with another.
    based_on (anndata.AnnData): The AnnData object that serves as the reference for canonicalizing the barcodes.
    batch_key (str): The key used to identify batches in the observations of the AnnData objects.

    Returns:
    anndata.AnnData: The AnnData object with canonicalized observation barcodes.
    """
    adata = adata.copy()
    adata.obs["original_barcode"] = adata.obs.index.to_series()
    mapper = get_barcode_mapper(adata=based_on, batch_key=batch_key)
    mapper = pd.merge(
        mapper,
        adata.obs[["original_barcode", batch_key]],
        how="right",
        on=["original_barcode", batch_key],
    )
    mapper["unique_barcode"] = mapper.apply(
        lambda row: (
            row["original_barcode"]
            if pd.isnull(row["unique_barcode"])
            else row["original_barcode"]
        ),
        axis=1,
    )
    if any(mapper["unique_barcode"].isnull()):
        raise ValueError(
            "Unable to canonicalize some barcodes. NAs would be introduced"
        )
    adata.obs.set_index(mapper["unique_barcode"])
    return adata


# %%
