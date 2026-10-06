import muon as mu
from typing import Union

# %%
def add_is_cell_annotation(
    mdata: mu.MuData,
    raw_key: str,
    filtered_key: str,
    key_added: str = "cr_is_cell_str",
    inplace=True,
) -> Union[mu.MuData, None]:
    if not inplace:
        mdata = mdata.copy()
    mdata[raw_key].obs[key_added] = (
        mdata[raw_key]
        .obs_names.to_series()
        .isin(mdata[filtered_key].obs_names.to_series())
        .apply(lambda ic: "Is cell" if ic else "Not cell")
    )
    if not inplace:
        return mdata

