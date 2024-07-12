# %%
import re

import anndata
import pandas as pd

BARCODE_PATTERN = re.compile("([ATCG]+)([-][0-9])([-][0-9])*")


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
def original_barcodes(current_barcodes: pd.Series | pd.Index) -> pd.Series | pd.Index:
    """
    Processes the current barcodes and returns the original barcodes.

    Args:
        current_barcodes (pd.Series|pd.Index): A pandas Series or Index object containing the current barcodes.

    Returns:
        pd.Series|pd.Index: A pandas Series or Index object (depending on input) containing the original barcodes.
    """
    if isinstance(current_barcodes, pd.Series):
        return current_barcodes.apply(_original_barcode)
    elif isinstance(current_barcodes, pd.Index):
        return pd.Index(current_barcodes.to_series().apply(_original_barcode))
    else:
        raise ValueError(f"Expected pd.Series, got: {type(original_barcodes)}")


# %%
def get_barcode_mapper(adata: anndata.AnnData, batch_key: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            k: v.to_numpy()
            for k, v in {
                "unique_barcode": adata.obs.index.to_numpy(),
                "original_barcode": original_barcodes(adata.obs.index),
                "batch_key": adata.obs[batch_key],
            }.items()
        }
    )

# %%
def canonicalize_barcodes(
    adata: anndata.AnnData, based_on: anndata.AnnData, batch_key: str
) -> anndata.AnnData:
    mapper = based_on.get_barcode_mapper()
    adata.obs["current_barcode"] = adata.obs.index

