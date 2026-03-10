# %%
import re
import os
# import toolz
# from functools import partial
import pandas as pd
from natsort import natsort_keygen
from itertools import chain
from pathlib import Path

# %% REGEX Patterns
SAMPLE_NUMBER = "S[0-9]+"
LANE_NUMBER = "L[0-9]+"
READ_NUMBER = "[IR][0-9]"
READ_NUMBER_NO_INDEX = "[R][0-9]"
SET_NUMBER = "001"
EXT = "[.]fastq[.]gz"

LANE_SPLIT_FASTQ = f"_{SAMPLE_NUMBER}_{LANE_NUMBER}_{READ_NUMBER}_{SET_NUMBER}{EXT}"
NO_LANE_SPLIT_FASTQ = f"_{SAMPLE_NUMBER}_{READ_NUMBER}_{SET_NUMBER}{EXT}"

LANE_SPLIT_FASTQ_NO_INDEX = f"_({SAMPLE_NUMBER})_({LANE_NUMBER})_({READ_NUMBER_NO_INDEX})_({SET_NUMBER}){EXT}"
NO_LANE_SPLIT_FASTQ_NO_INDEX = f"_({SAMPLE_NUMBER})_({READ_NUMBER_NO_INDEX})_({SET_NUMBER}){EXT}"

# %%
def is_lane_split_fastq(fn: str) -> bool:
    if re.search(LANE_SPLIT_FASTQ, fn) is not None:
        return True
    else:
        if re.search(NO_LANE_SPLIT_FASTQ, fn) is not None:
            return False
        else:
            raise ValueError(f"filename doesn't matter either expected pattern: {fn}")

# %%
def filter_file_paths(fs: list[os.PathLike], matching_pattern=LANE_SPLIT_FASTQ_NO_INDEX):
    filt = []
    for f in fs:
        if re.search(matching_pattern, f) is not None:
            filt.append(f)
    return filt

# %%

def get_sample_id(f, illumina_suffix_pattern=LANE_SPLIT_FASTQ_NO_INDEX):
    return re.sub(illumina_suffix_pattern, "", f)

# %%
def split_filenames(fs_str: list[str], pat: str):
    ldict = []
    for f in fs_str:
        s = re.search(pat, f)
        if s is None:
            raise ValueError("re.search is None when it really shouldn't be able to be")
        sid = get_sample_id(f, pat)
        ldict.append({k:v for k,v in zip(["sid", "sn", "ln", "rn", "set", "path"], chain([sid], s.groups(), [f]))})
    return pd.DataFrame(ldict).sort_values(["sid", "sn", "ln", "rn"])

# %%
def combine_lane_split_fastqs(
    fastq_dir: os.PathLike,
    output_dir: os.PathLike,
    # keep_originals: bool = False,
    exclude_index_reads=True
):
    if exclude_index_reads:
        pat = LANE_SPLIT_FASTQ_NO_INDEX
    else:
        pat = LANE_SPLIT_FASTQ
    fs_str = [p.name for p in Path(fastq_dir).iterdir()]
    fs_str = filter_file_paths(fs_str, pat)
    df = split_filenames(fs_str, pat)
    # grouped = {}
    for _key,grp in df.groupby(["sid", "rn"]):
        sorted = grp.sort_values(["ln"], key=natsort_keygen())
        assert len(sorted["sid"].unique())==1
        assert len(sorted["sn"].unique())==1
        assert len(sorted["rn"].unique())==1
        sid = sorted["sid"].unique()[0]
        sn = sorted["sn"].unique()[0]
        rn = sorted["rn"].unique()[0]
        out_base = f"{sid}_{sn}_{rn}_001.fastq.gz"
        with open(os.path.join(output_dir, out_base), 'wb') as combined:
            for f in sorted["path"]:
                with open(os.path.join(fastq_dir, f), 'rb') as fastq:
                    combined.write(fastq.read())

# %%
if __name__ == "__main__":
    pa = "/home/amsesk/super2/jayme_shiv/data/cr_atac_mkfastq/outs/fastq_test/"
    combine_lane_split_fastqs(pa, output_dir="/home/amsesk/super2/jayme_shiv/data/cr_atac_mkfastq/outs/fastq_test/")

