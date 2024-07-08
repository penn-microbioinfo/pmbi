import os
from pathlib import Path

import anndata
import mudata
import numpy as np
import pandas as pd
import scanpy as sc


def get_key_default(path: Path) -> str:
    """
    Extract a key from a file path by removing the suffix and splitting by underscore.

    Parameters:
        path (Path): The file path from which to extract the key.

    Returns:
        str: The extracted key.

    Example:
        # Assuming 'file_path' is a pathlib.Path object
        key = get_key_default(file_path)
    """

    return path.name.replace(f"{path.suffix}", "").split("_")[0]

def read_gct(gctpath: os.PathLike) -> pd.DataFrame:
    """
    Read data from a GCT file and return as a pandas DataFrame.

    Parameters:
        gctpath (os.PathLike): The path to the GCT file.

    Returns:
        pd.DataFrame: The data from the GCT file as a pandas DataFrame.

    Example:
        # Assuming 'gct_file_path' is a path to a GCT file
        dataframe = read_gct(gct_file_path)
    """

    with open(gctpath, "r") as gct:

        # Skip the GCT header
        for _ in range(0, 2):
            gct.readline()

        df = pd.read_csv(gct, sep="\t").transpose()
        df.columns = df.iloc[0, :]
        df = df.iloc[3:, :]
        df = df.loc[df.index.str.match("[ACGT]+[-][0-9]")]
        df = df.astype(np.float64)

        return df

