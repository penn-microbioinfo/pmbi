import pmbi.anndata.random
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

