import pandas as pd

import pmbi.anndata.util


def test_original_barcode_leaves_single_trailing():
    current_barcode = "TGAGCATAGTTCCACA-1"
    original_barcode = pmbi.anndata.util._original_barcode(
        current_barcode, pmbi.anndata.util.BARCODE_PATTERN
    )
    assert original_barcode == "TGAGCATAGTTCCACA-1"


def test_original_barcode_removes_second_trailing():
    current_barcode = "TGAGCATAGTTCCACA-1-1"
    original_barcode = pmbi.anndata.util._original_barcode(
        current_barcode, pmbi.anndata.util.BARCODE_PATTERN
    )
    assert original_barcode == "TGAGCATAGTTCCACA-1"


def test_original_barcode_removes_multi_trailing():
    current_barcode = "TGAGCATAGTTCCACA-1-1-1-1"
    original_barcode = pmbi.anndata.util._original_barcode(
        current_barcode, pmbi.anndata.util.BARCODE_PATTERN
    )
    assert original_barcode == "TGAGCATAGTTCCACA-1"


def test_original_barcodes_correctly_replace_series():
    current_barcodes = pd.Series(
        ["TGAGCATAGTTCCACA-1-1-1-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1-1"]
    )
    original_barcodes = pmbi.anndata.util.original_barcodes(current_barcodes)
    assert original_barcodes.equals(
        pd.Series(["TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1"])
    )
    assert isinstance(original_barcodes, pd.Series)

def test_original_barcodes_correctly_replace_index():
    current_barcodes = pd.Index(
        ["TGAGCATAGTTCCACA-1-1-1-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1-1"]
    )
    original_barcodes = pmbi.anndata.util.original_barcodes(current_barcodes)
    assert original_barcodes.equals(
        pd.Index(["TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1", "TGAGCATAGTTCCACA-1"])
    )
    assert isinstance(original_barcodes, pd.Index)
