import numpy as np


# %%
def random_barcode(length: int = 16) -> str:
    """
    Generate a random 10X-style barcode of a specified length.

    Parameters:
    length (int): The length of the barcode string to generate. Default is 16.

    Returns:
    str: A randomly generated barcode of the specified length.
    """
    nucl = "ATCG"
    randidx = np.random.randint(0, 4, size=(length,))
    randbc = "".join([nucl[i] for i in randidx])
    return f"{randbc}-1"


# %%
def random_barcodes(n: int, length: int = 16, assert_unique=False) -> list:
    """
    Generate a list of random 10X-style barcodes.

    Parameters:
    n (int): The number of barcodes to generate.
    length (int, optional): The length of each barcode. Defaults to 16.
    assert_unique (bool, optional): If True, ensure all generated barcodes are unique. Defaults to False.

    Returns:
    list: A list of randomly generated barcodes.
    """
    if assert_unique:
        bcs = []
        for _ in range(0, n):
            while True:
                bc = random_barcode(length)
                if bc in bcs:
                    continue
                else:
                    bcs.append(bc)
                    break
        return bcs
    else:
        return [random_barcode(length) for _ in range(0, n)]


# %%
def random_adata(shape: tuple[int, int], nobs_cols=5, nvar_cols=5, assert_unique_barcodes = False) -> anndata.AnnData:
    return anndata.AnnData(
        X=scipy.sparse.random(m=shape[0], n=shape[1]),
        obs=pd.DataFrame(
            np.random.random(size=(shape[0], nobs_cols)),
            index=random_barcodes(n=shape[0], assert_unique = assert_unique_barcodes),
        ),
        var=pd.DataFrame(np.random.random(size=(shape[1], nvar_cols))),
    )

