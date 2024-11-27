import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pmbi.bio.dna import revcomp
from pmbi.plotting import Paneler

# %%
ss = pd.read_csv("/home/ubuntu/projmnt/betts-dl/287879594/PmbiSampleSheet.csv", sep=",")
tt_idx = pd.read_csv(
    "/home/ubuntu/dev/pmbi/data/10X_idx/Dual_Index_Kit_TT_Set_A.csv",
    sep=",",
    skiprows=3,
)
tn_idx = pd.read_csv(
    "/home/ubuntu/dev/pmbi/data/10X_idx/Dual_Index_Kit_TN_Set_A.csv",
    sep=",",
    skiprows=3,
)

idx = pd.concat([tt_idx, tn_idx])
idx["index2_workflow_b_rc"] = idx["index2_workflow_b(i5)"].apply(lambda i: revcomp(i))

# %%
idx = idx.set_index("index_name")
used_idx = idx.loc[ss["Index"]][
    ["index(i7)", "index2_workflow_b(i5)", "index2_workflow_b_rc"]
]
used_idx = used_idx.rename(
    columns={
        "index(i7)": "i7",
        "index2_workflow_b(i5)": "i5",
        "index2_workflow_b_rc": "i5_rc",
    }
)

# %% Make bcl2fastq samplesheet instead
pd.concat([used_idx, ss.set_index("Index")], axis=1)[
    ["Sample", "i7", "i5"]
].reset_index(drop=True).to_csv(
    "/home/ubuntu/projmnt/betts-dl/287879594/PmbiBcl2fastqSampleSheet.csv",
    sep=",",
    index=False,
)

# %%
idx_count_dir = "/home/ubuntu/projmnt/betts-dl/287879594/fastq_mux/"
idx_count_f = [
    os.path.join(idx_count_dir, x)
    for x in os.listdir(idx_count_dir)
    if x.endswith("_idxcounts")
]
i1_count_f = {
    re.search("_(L00[0-9])_", x).group(1): x for x in idx_count_f if "_I1_" in x
}
i2_count_f = {
    re.search("_(L00[0-9])_", x).group(1): x for x in idx_count_f if "_I2_" in x
}
# %%
all_counts = {}
for k, v in {"I1": i1_count_f, "I2": i2_count_f}.items():
    counts = []
    for l, f in v.items():
        these_counts = pd.read_csv(f, sep=",")
        these_counts = these_counts.set_index("barcode", drop=True)
        these_counts = these_counts.rename(columns={"count": l})
        counts.append(these_counts)
    comb = pd.concat(counts, axis=1, ignore_index=False)
    comb = comb.replace(np.nan, 0)
    comb = comb.assign(total=lambda r: r.L001 + r.L002 + r.L003 + r.L004)
    all_counts[k] = comb

# %%
all_counts

# %%
i1_top200 = i1_counts.nlargest(n=200, columns="count")
i2_top200 = i2_counts.nlargest(n=200, columns="count")

i1_counts
# %%
pres = pd.DataFrame(
    {
        "expected_i7_present": used_idx["i7"].isin(i1_top200["barcode"]),
        "expected_i5_present": used_idx["i5"].isin(i2_top200["barcode"]),
    }
)
ss[ss.Index.isin(used_idx[pres["expected_i7_present"].values].index)].to_csv(
    "/home/ubuntu/projmnt/betts-dl/287879594/PmbiSampleSheet_withPresentBarcodes.csv",
    sep=",",
    index=False,
)
np.where(pres["expected_i7_present"])[0].shape
# %%

i1_missing.to_csv("i1_missing.txt")
i1_missing = i1_top200[~i1_top200["barcode"].isin(used_idx["i7"])]["barcode"]
i1_missing.values
any([x in idx["index(i7)"].values for x in i1_missing.values])

i2_top200[~i2_top200["barcode"].isin(used_idx["i5"])]

not_pres = pd.DataFrame(
    {
        "expected_i7_present": used_idx["i7"].isin(i1_top200["barcode"]),
        "expected_i5_present": used_idx["i5"].isin(i2_top200["barcode"]),
    }
)
# %%
np.where([x in i1_counts["barcode"].values for x in used_idx["i7"].valu])


# %%
panel = Paneler(1, 1, figsize=(2, 6), format="pdf")
# panel.next_ax().hist(x=np.log(top50["count"]), bins=100)
panel.next_ax().imshow(pres, aspect="auto")
pres
panel.current_ax.set_xticks(
    np.arange(len(pres.columns)), labels=pres.columns, size=3, rotation=90
)
panel.current_ax.set_yticks(np.arange(len(pres.index)), labels=pres.index, size=3)
panel.fig.savefig("/srv/http/betts/coculture/expected_barcode_present.pdf")


# %%
