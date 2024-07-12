import logging
import os

import anndata
import numpy
import scanpy as sc
import scipy.io
from pathlib import Path
from pmbi.io import get_key_default

import pmbi.anndata.get


def to_mtx(
    adata: anndata.AnnData, output_prefix: str, layer: str | None = None
) -> None:
    counts_out = f"{output_prefix}_counts.mtx"
    obs_out = f"{output_prefix}_obs.csv"
    var_out = f"{output_prefix}_var.csv"
    logging.warning("`pmbi.convert.adata_write_mtx` currently ignores `adata.uns`")
    with open(counts_out, "wb") as cout:
        scipy.io.mmwrite(cout, pmbi.anndata.get.counts(adata, layer=layer))

    adata.obs.to_csv(path_or_buf=obs_out, sep=",", index_label="barcode")
    adata.var.to_csv(path_or_buf=var_out, sep=",", index_label="feature")

def write_counts(adata: anndata.AnnData, file: str, sep = ',', layer: str|None = None):
    if layer is None:
        df = pd.DataFrame(adata.X.toarray())

    else:
        df = pd.DataFrame(adata.layers[layer].toarray())

    df.index = adata.obs_names
    df.columns = adata.var_names

    df.to_csv(file, sep = sep)


def read_matrix(path: Path, **kwargs) -> anndata.AnnData:
    """
    Read matrix data from file and return as an AnnData object.
    Currently supported formats:
        h5ad
        10X h5
        10X mtx (directory)

    Parameters:
        path (Path): The path to the matrix directory or file.
        **kwargs: Additional keyword arguments to pass to the reader function.

    Returns:
        anndata.AnnData: The matrix data as an AnnData object.

    Raises:
        ValueError: If the file type is not supported.

    Example:
        # Assuming 'file_path' is a pathlib.Path object
        matrix_data = read_matrix(file_path)
    """
    if os.path.isdir(path):
        return sc.read_10x_mtx(path, **kwargs)
    else:
        if path.suffix == ".h5":
            return sc.read_10x_h5(str(path), **kwargs)

        elif path.suffix == ".h5ad":
            return anndata.read_h5ad(str(path), **kwargs)
        else:
            raise ValueError(f"Unsupported filetype: {path.suffix}")


def read_matrix_multi(
    paths: list[Path], getkey=get_key_default
) -> dict[str, anndata.AnnData]:
    """
    Read matrix data from multiple files and return as a dictionary of AnnData objects.
    Supports the same formats as `read_matrix`

    Parameters:
        paths (list[Path]): A list of paths to the matrix files or directories.
        getkey (callable): A function to extract keys from file paths. Defaults to `get_key_default`.

    Returns:
        dict[str, anndata.AnnData]: A dictionary mapping keys to AnnData objects.

    Example:
        # Assuming 'file_paths' is a list of pathlib.Path objects
        data_dict = read_matrix_multi(file_paths)
    """

    adatas = {}

    for path in paths:
        adatas[getkey(path)] = read_matrix(path)

    return adatas


def read_h5ad_multi(
    paths: list[Path], getkey=get_key_default
) -> dict[str, anndata.AnnData]:
    """
    Read h5ad files from multiple paths and return as a dictionary of AnnData objects.

    Parameters:
        paths (list[Path]): A list of paths to the h5ad files.
        getkey (callable): A function to extract keys from file paths. Defaults to get_key_default.

    Returns:
        dict[str, anndata.AnnData]: A dictionary mapping keys to AnnData objects.

    Example:
        # Assuming 'file_paths' is a list of pathlib.Path objects
        data_dict = read_h5ad_multi(file_paths)
    """

    adatas = {}

    for path in paths:
        adatas[getkey(path)] = sc.read_h5ad(path)

    return adatas


def write_h5ad_multi(
    adatas: dict[str, anndata.AnnData], suffix: str, outdir: str
) -> None:
    """
    Write multiple AnnData objects to h5ad files with specified suffix and output directory.

    Parameters:
        adatas (dict[str, anndata.AnnData]): A dictionary mapping keys to AnnData objects.
        suffix (str): The suffix to append to the filenames.
        outdir (str): The output directory to save the files.

    Returns:
        None

    Example:
        # Assuming 'data_dict' is a dictionary mapping keys to AnnData objects
        write_h5ad_multi(data_dict, suffix="normalized", outdir="/path/to/output_directory")
    """

    for key, adata in adatas.items():
        adata.write_h5ad(Path(os.path.join(outdir, f"{key}_{suffix}.h5ad")))


