import logging
import h5py
import os
import pickle
from pathlib import Path

import anndata
import numpy
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse
from scirpy.io._convert_anndata import from_airr_cells
from scirpy.io._datastructures import AirrCell

import pmbi.anndata.get
from pmbi.airr.schema import calls_to_imgt_locus_names, get_airr_schema
from pmbi.io import get_key_default


def write_mtx(
    adata: anndata.AnnData, output_dir: os.PathLike, layer: str | None = None, write_obs = True, write_var = True
) -> None:
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents = True, exist_ok = False) # Will raies FileExistsError if exists
    obs_out = output_dir.joinpath("obs.csv")
    var_out = output_dir.joinpath("var.csv")
    mat = pmbi.anndata.get.counts(adata, layer=layer)
    if scipy.sparse.issparse(mat):
        counts_out = output_dir.joinpath("counts.mtx")
        with open(counts_out, "wb") as cout:
            scipy.io.mmwrite(cout, mat)
    else:
        counts_out = output_dir.joinpath("counts.h5")
        h5f = h5py.File(counts_out, 'w')
        h5f.create_dataset("counts", data=mat)
        h5f.close()


    if write_obs:
        adata.obs.to_csv(path_or_buf=obs_out, sep=",", index_label="barcode")

    if write_var:
        adata.var.to_csv(path_or_buf=var_out, sep=",", index_label="feature")

def export(
    adata: anndata.AnnData, output_dir: os.PathLike
    ):
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents = True, exist_ok = False) # Will raies FileExistsError if exists

    # Write count mtx from anndata.X
    write_mtx(adata, output_dir)

    # Wrtie out layers as mtx
    layers_out = output_dir.joinpath("layers")
    layers_out.mkdir()
    for layer in adata.layers:
        write_mtx(adata, output_dir = layers_out.joinpath(layer), layer=layer, write_obs=False, write_var=False)

    # Wrtie out obsm items into on h5 file
    h5f = h5py.File(output_dir.joinpath("obsm.h5"), 'w')
    for key,value in adata.obsm.items():
        h5f.create_dataset(key, data=value)
    h5f.close()


def write_counts(adata: anndata.AnnData, file: str, sep=",", layer: str | None = None):
    if layer is None:
        df = pd.DataFrame(adata.X.toarray())

    else:
        df = pd.DataFrame(adata.layers[layer].toarray())

    df.index = adata.obs_names
    df.columns = adata.var_names

    df.to_csv(file, sep=sep)


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
    path = Path(path)
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


def pickle_piece(
    adata: anndata.AnnData, outer_key: str, inner_key: str, output_dir: os.PathLike
) -> None:
    pkl_path = os.path.join(output_dir, f"{outer_key}__{inner_key}.pkl")
    with open(pkl_path, "wb") as pkl_f:
        pickle.dump(getattr(adata, outer_key)[inner_key], pkl_f)


def pickle_pieces(
    adata: anndata.AnnData,
    pickle_what: dict[str, list[str]],
    output_dir: os.PathLike,
) -> None:
    for outer_key, inner_key_list in pickle_what.items():
        for inner_key in inner_key_list:
            pickle_piece(
                adata=adata,
                outer_key=outer_key,
                inner_key=inner_key,
                output_dir=output_dir,
            )


def load_pickle_piece(path: os.PathLike) -> object:
    # "alignment_nt_annot.csv.gz" %%
    with open(path, "rb") as pkl_f:
        piece = pickle.load(pkl_f)
    return piece


def from_airr(
    airr_df: pd.DataFrame, schema: str, sample_key: str | None = None
) -> anndata.AnnData:
    cells = []
    airr_df["productive"] = [
        False if p.upper() == "F" else True for p in airr_df["productive"]
    ]
    for row in airr_df.itertuples():
        cell = AirrCell(cell_id=str(row.sequence_id))
        ac = AirrCell.empty_chain_dict()
        ac["locus"] = calls_to_imgt_locus_names([row.v_call, row.d_call, row.j_call])
        ac.update(
            {
                req_field: getattr(row, req_field)
                for req_field in get_airr_schema(schema=schema)["required"]
            }
        )
        cell.add_chain(ac)
        cells.append(cell)
    adata = from_airr_cells(cells)
    if sample_key is not None:
        adata.obs["sample_key"] = sample_key
    return adata
