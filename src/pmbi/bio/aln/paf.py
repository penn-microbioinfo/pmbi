# %%
import numpy as np
import pandas as pd
import os

# %%
class Tag:
    def __init__(self, key, type_, value):
        self.key = key
        self.type_ = type_
        self.value = value
    
    def __repr__(self):
        return f"{self.key}:{self.type_}:{self.value}"


# %%
PAF_COLNAMES = [
    "qname",
    "qlen",
    "qs",
    "qe",
    "strand",
    "tname",
    "tlen",
    "ts",
    "te",
    "nmatching",
    "nbases",
    "mapq"
]

# %%
def read_paf_to_df(paf_path: os.PathLike) -> pd.DataFrame:
    paf = pd.read_csv(paf_path, sep="\t")
    req_cols = paf.iloc[:,list(range(0,12))]
    req_cols.columns=PAF_COLNAMES
    tags = paf.iloc[:,list(range(12,paf.shape[1]))]
    tags = tags.apply(split_tags, axis=1)
    comb = req_cols
    comb = comb.assign(tags=tags)
    return comb

# %%
def split_tags(row):
    tag_dict = {}
    for t in row.reset_index(drop=True):
        if pd.isnull(t):
            continue
        spl = t.split(":")
        if len(spl) == 3:
            spl = tuple(spl)
            tag_dict[spl[0]] = Tag(*spl)
        else:
            raise ValueError(spl)
    return tag_dict

