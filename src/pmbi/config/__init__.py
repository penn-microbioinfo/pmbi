from munch import Munch
import pandas as pd
import tomllib

# %%
def import_config(path: str) -> Munch:
    with open(path, "rb") as c:
        config = Munch.fromDict(tomllib.load(c))
    return config

def table_to_dataframe(l: list[Munch]):
    if not isinstance(l, list) or not all([isinstance(x, Munch) for x in l]):
        raise ValueError("l must be list[Munch]")
    return pd.DataFrame(map(lambda x: x.__dict__, l))

