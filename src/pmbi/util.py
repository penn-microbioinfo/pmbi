import awkward as ak
from pprint import pprint

def ak_print(arr: ak.Array) -> None:
    pprint(ak.to_list(arr))

