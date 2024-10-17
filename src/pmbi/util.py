import awkward as ak
from pprint import pprint
import pandas as pd

def ak_print(arr: ak.Array) -> None:
    pprint(ak.to_list(arr))

def invert_dict(d, return_type = "series") -> pd.Series:
    s = pd.DataFrame({"clonotype_id": pd.Series(d.keys()), "barcode": pd.Series(d.values())})
    if return_type == "series":
        return s.explode("barcode").set_index("barcode").loc[:,"clonotype_id"]
    elif return_type == "dict":
        raise NotImplementedError(f"Unimplemented return_type: {return_type}")
    else:
        raise NotImplementedError(f"Unimplemented return_type: {return_type}")

def dict_try_add(d, k, v, f = lambda x: x.):
    if k in d:
        d[k] = v
    else:

