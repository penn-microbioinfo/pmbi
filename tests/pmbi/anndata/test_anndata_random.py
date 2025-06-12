import pmbi.anndata.random
import anndata
from pmbi.anndata.util import BARCODE_PATTERN
import re

# %%
def test_random_barcode_matches_pattern():
    bc = pmbi.anndata.random.random_barcode()
    assert re.match(BARCODE_PATTERN, bc)

# %%
def test_random_barcodes_assert_unique():
    bcs = pmbi.anndata.random.random_barcodes(n = 1024, length = 5, assert_unique = True)
    assert len(bcs) == len(set(bcs))

# %%
def test_random_barcodes_all_match_pattern():
    bcs = pmbi.anndata.random.random_barcodes(n = 1024, length = 5, assert_unique = True)
    assert all([re.match(BARCODE_PATTERN, bc) for bc in bcs])

def test_random_adata():
    adata = pmbi.anndata.random.random_adata(shape = (10, 10), nobs_cols=5, nvar_cols=5)
    assert adata.shape == (10,10)
    assert adata.obs.shape == (10,5)
    assert adata.var.shape == (10,5)
    assert isinstance(adata, anndata.AnnData)
