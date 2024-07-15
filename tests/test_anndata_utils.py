import importlib

import anndata
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import random

import pmbi.anndata.util

importlib.reload(pmbi.anndata.util)


# %%
def test_original_barcode_leaves_single_trailing():
    current_barcode = "TGAGCATAGTTCCACA-1"
    original_barcode = pmbi.anndata.util._original_barcode(
        current_barcode, pmbi.anndata.util.BARCODE_PATTERN
    )
    assert original_barcode == "TGAGCATAGTTCCACA-1"


# %%
def test_original_barcode_removes_second_trailing():
    current_barcode = "TGAGCATAGTTCCACA-1-1"
    original_barcode = pmbi.anndata.util._original_barcode(
        current_barcode, pmbi.anndata.util.BARCODE_PATTERN
    )
    assert original_barcode == "TGAGCATAGTTCCACA-1"


# %%
def test_original_barcode_removes_multi_trailing():
    current_barcode = "TGAGCATAGTTCCACA-1-1-1-1"
    original_barcode = pmbi.anndata.util._original_barcode(
        current_barcode, pmbi.anndata.util.BARCODE_PATTERN
    )
    assert original_barcode == "TGAGCATAGTTCCACA-1"


# %%
def test_original_barcodes_correctly_replace_series():
    current_barcodes = pd.Series(
        ["TGAGCATAGTTCCACA-1-1-1-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1-1"]
    )
    original_barcodes = pmbi.anndata.util.original_barcodes(current_barcodes)
    assert original_barcodes.equals(
        pd.Series(["TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1"])
    )
    assert isinstance(original_barcodes, pd.Series)


# %%
def test_original_barcodes_correctly_replace_index():
    current_barcodes = pd.Index(
        ["TGAGCATAGTTCCACA-1-1-1-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1-1"]
    )
    original_barcodes = pmbi.anndata.util.original_barcodes(current_barcodes)
    assert original_barcodes.equals(
        pd.Series(
            ["TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1"],
            index=current_barcodes,
        )
    )
    assert isinstance(original_barcodes, pd.Series)


# %%
def test_shared_barcodes_correctly_id_shared():
    bc1 = [
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACG-1",
        "TGAGCATAGTTCCACC-1",
        "TGAGCATAGTTCCAGA-1",
    ]
    bc2 = [
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACG-1",
        "TGAGCATAGTTCCACC-1",
        "TGAGCATAGTTCCAGA-1",
    ]
    adata1 = anndata.AnnData(
        X=random(m=4, n=10), obs=pd.DataFrame(np.random.random(size=(4, 3)), index=bc1)
    )
    adata2 = anndata.AnnData(
        X=random(m=4, n=10), obs=pd.DataFrame(np.random.random(size=(4, 3)), index=bc2)
    )
    pmbi.anndata.util.shared_barcodes([adata1, adata2])
    assert pmbi.anndata.util.shared_barcodes([adata1, adata2]).equals(
        pd.Series(
            [
                "TGAGCATAGTTCCACT-1",
                "TGAGCATAGTTCCACG-1",
                "TGAGCATAGTTCCACC-1",
                "TGAGCATAGTTCCAGA-1",
            ]
        )
    )


# %%
def test_shared_barcodes_correctly_id_unshared():
    bc1 = [
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACG-1",
        "TGAGCATAGTTCCACC-1",
        "TGAGCATAGTTCCAGA-1",
    ]
    bc2 = [
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACG-1",
        "TGAGCATAGTTCCACC-1",
        "TGAGCATAGTTCCATT-1",
    ]
    adata1 = anndata.AnnData(
        X=random(m=4, n=10), obs=pd.DataFrame(np.random.random(size=(4, 3)), index=bc1)
    )
    adata2 = anndata.AnnData(
        X=random(m=4, n=10), obs=pd.DataFrame(np.random.random(size=(4, 3)), index=bc2)
    )
    assert pmbi.anndata.util.shared_barcodes([adata1, adata2]).equals(
        pd.Series(["TGAGCATAGTTCCACT-1", "TGAGCATAGTTCCACG-1", "TGAGCATAGTTCCACC-1"])
    )


# %%
def test_shared_barcodes_correctly_id_shared_many_anndata():
    adatas = []
    bcs = [
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACG-1",
        "TGAGCATAGTTCCACC-1",
        "TGAGCATAGTTCCAGA-1",
    ]
    for _ in range(0, 30):
        adatas.append(
            anndata.AnnData(
                X=random(m=4, n=10),
                obs=pd.DataFrame(np.random.random(size=(4, 3)), index=bcs),
            )
        )
    assert pmbi.anndata.util.shared_barcodes(adatas).equals(
        pd.Series(
            [
                "TGAGCATAGTTCCACT-1",
                "TGAGCATAGTTCCACG-1",
                "TGAGCATAGTTCCACC-1",
                "TGAGCATAGTTCCAGA-1",
            ]
        )
    )


def test_shared_barcodes_correctly_raise_value_exception_on_non_unique_bcs():
    bc1 = [
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACC-1",
        "TGAGCATAGTTCCAGA-1",
    ]
    bc2 = [
        "TGAGCATAGTTCCACT-1",
        "TGAGCATAGTTCCACG-1",
        "TGAGCATAGTTCCACC-1",
        "TGAGCATAGTTCCATT-1",
    ]
    adata1 = anndata.AnnData(
        X=random(m=4, n=10), obs=pd.DataFrame(np.random.random(size=(4, 3)), index=bc1)
    )
    adata2 = anndata.AnnData(
        X=random(m=4, n=10), obs=pd.DataFrame(np.random.random(size=(4, 3)), index=bc2)
    )
    with pytest.raises(ValueError):
        pmbi.anndata.util.shared_barcodes([adata1, adata2])


# %%
def test_get_barcode_mapper_produces_correct_map():
    adata = pmbi.anndata.random.random_adata(shape=(10, 10))
    adata.obs["batch_key"] = "batch1"
    bcs = adata.obs.index
    mapper = pmbi.anndata.util.get_barcode_mapper(adata=adata, batch_key="batch_key")
    should_be = pd.DataFrame(
        {
            "unique_barcode": bcs,
            "original_barcode": bcs,
            "batch_key": adata.obs["batch_key"],
        }
    ).reset_index(drop=True)
    assert mapper.equals(should_be)


def test_get_barcode_mapper_raises_value_erros_when_obs_names_not_unique():
    adata = pmbi.anndata.random.random_adata(shape=(5, 5))
    adata.obs.index = [list(adata.obs.index)[i] for i in [0, 0, 0, 0, 0]]
    adata.obs["batch_key"] = "batch1"
    with pytest.raises(ValueError):
        pmbi.anndata.util.get_barcode_mapper(adata=adata, batch_key="batch_key")


def test_obs_canonicalize_barcodes():
    adata1 = pmbi.anndata.random.random_adata(shape=(20, 20))
    adata1.obs["batch_key"] = ["batch1"] * 10 + ["batch2"] * 10
    adata1.obs.index = [list(adata1.obs.index)[n] for n in list(range(0, 10)) * 2]
    adata2 = adata1.copy()
    adata1.obs_names_make_unique()
    assert adata1.obs.index.equals(
        pmbi.anndata.util.obs_canonicalize_barcodes(
            adata2, adata1, batch_key="batch_key"
        ).obs.index
    )
