from munch import Munch
import tomllib

# %%
def import_config(path: str) -> Munch:
    with open(path, "rb") as c:
        config = Munch.fromDict(tomllib.load(c))
    return config

