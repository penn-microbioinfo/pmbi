import re
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from anndata import AnnData
from muon import MuData
from toolz import reduce

import pmbi.anndata.io as aio
from pmbi.bio.dna import revcomp


def read_scc_proc_mtx(mtxdir: Path) -> AnnData:
    adata = aio.read_mtx(
        mtx_path=mtxdir.joinpath("output.mtx"),
        obs_names_path=mtxdir.joinpath("output.barcodes.txt"),
        var_names_path=mtxdir.joinpath("output.genes.txt"),
    )
    return adata


def construct_mudata_from_matrices(
    paths: dict[str, Union[str, dict]],
    filter_barcodes_for: Union[list[str], None] = None,
    reverse_complement_barcodes: Union[list[bool], None] = None
):
    adatas = {}
    if filter_barcodes_for is not None and any([f not in paths for f in filter_barcodes_for]):
        raise ValueError(f"unexpected key in `filter_barcodes_for`")
    if reverse_complement_barcodes is None:
        reverse_complement_barcodes = [False] * len(paths)
    else:
        assert len(reverse_complement_barcodes) == len(paths)
    for (key, val), rc in zip(paths.items(), reverse_complement_barcodes):
        if isinstance(val, Path):
            adata = aio.read_matrix(val, gex_only=False)
        elif isinstance(val, dict):
            assert all(
                [
                    k in val.keys()
                    for k in ["mtx_path", "obs_names_path", "var_names_path"]
                ]
            )
            assert len(val) == 3
            adata = aio.read_mtx(**val)
        else:
            raise ValueError(f"unexpected type in kwargs: {type(val)}")
        trailing_number_barcode_pattern = "(?:[-][0-9]+)+$"
        if any([re.search(trailing_number_barcode_pattern, b) is not None for b in adata.obs_names]):
            adata.obs = adata.obs.assign(original_barcode=lambda x: x.index)
            adata.obs_names = [re.sub(trailing_number_barcode_pattern, "", b) for b in adata.obs_names]
        if rc:
            adata.obs_names = [revcomp(x) for x in adata.obs_names.to_list()]
        adatas[key] = adata
    if filter_barcodes_for is not None:
        all_barcodes = [a.obs_names.to_numpy() for k,a in adatas.items() if k in filter_barcodes_for]
        # keep_barcodes = reduce(np.intersect1d, all_barcodes, adatas["atac"].obs_names.to_numpy())
        keep_barcodes = reduce(np.intersect1d, all_barcodes)
        for key in filter_barcodes_for:
            adatas[key] = adatas[key][adatas[key].obs_names.isin(keep_barcodes), :]
        assert len(set([a.shape[0] for k,a in adatas.items() if k in filter_barcodes_for]))==1
    return MuData(adatas)

# %%
