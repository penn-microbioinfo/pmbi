# %%
import scipy.sparse as sparse
import pandas as pd
import numpy as np
from pmbi.plotting import Paneler
import matplotlib.pyplot as plt
import os
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--reference", help="Path to genome reference, requires .fai index")
parser.add_argument("-d", "--depth", help="Path to `samtools depth` output file.")
parser.add_argument("-o", "--output_prefix", help="Path to `samtools depth` output file.")
args = parser.parse_args()

def get_ref_index(refpath: str):
    idxpath = Path(f"{refpath}.fai")
    if not os.path.isfile(idxpath):
        raise ValueError(f"Unable to locate reference fai index: {idxpath}")
    else:
        return idxpath


chr_sizes = pd.read_csv(get_ref_index(args.reference), sep="\t", header=None).iloc[:,0:2]
chr_sizes.columns=["chr","size"]
chr_sizes=chr_sizes.set_index("chr")


dep = pd.read_csv(args.depth, sep="\t", header=None, compression="gzip")
dep.columns = ["chr", "pos", "depth"]

# %%
for this_chr in dep["chr"].unique():
    print(this_chr)
    plt.close("all")
    plt.clf()
    this_chr_size = chr_sizes.loc[this_chr,"size"]
    sub = dep[dep["chr"]==this_chr]
    depth_arr = sub["depth"].to_numpy()
    cols_arr = np.full(len(sub["pos"]), 0)
    pos_arr = sub["pos"].to_numpy()
    D = sparse.coo_array((depth_arr, (pos_arr, cols_arr)), shape=(this_chr_size,1))
    D_arr = D.toarray()
    panel = Paneler(1,1,(8,8))
    panel.next_ax().scatter(sub["pos"].values, sub["depth"].values, s=0.75, marker=".", edgecolors="none")
    panel.fig.savefig(f"{args.output_prefix}/AES103_Pip1_depth_{this_chr}.png", dpi=1500)

sys.exit()

# %%
r = np.array([200e6, 250e6]).astype(np.int64)
panel = Paneler(1,1,(12,12))
# panel.next_ax().plot(range(0,10000), D_arr[0:10000,0], linewidth=0.1)
panel.next_ax().scatter(range(r[0], r[1]), D_arr[r[0]:r[1],0], s=0.5, marker=".")
panel.fig.savefig("/storage/anat/figures/AES103_Pip1_depth.png")
# %%
D.toarray()
len(sub["depth"])

# %%
np.full(5, 10)
rows = np.array([0,1,3])
cols = np.array([1,2,4])
data = np.array([10,20,30])
sparse.coo_matrix((data,(rows,cols)), shape=(4,5)).toarray()

# %%
loci_sam = pd.read_csv("/storage/anat/sam/loci_Pip1.sorted.sam", sep="\t", header=None, usecols=[0,2], names=["query","reference"])
loci_sam["count"] = 1
loci_mat = loci_sam.pivot_table(index="query", values="count", columns="reference", aggfunc=sum).fillna(0).drop(columns="*")
loci_mat.to_csv("/storage/anat/loci_Pip1_chrom_mat.csv", sep=",")
                       
# %%
loci_sam = pd.read_csv("/storage/anat/sam/loci_Pip1.sorted.sam", sep="\t", header=None, usecols=[0,2,3,9], names=["query","reference","pos","seq"])
loci_sam=loci_sam[loci_sam["reference"]!="*"]
loci_sam["seqlen"] = loci_sam["seq"].apply(len)
loci_sam["start"] = loci_sam["pos"]-loci_sam["seqlen"]
loci_sam["stop"] = loci_sam["pos"]+loci_sam["seqlen"]
loci_sam[["reference", "start", "stop", "query"]].to_csv("/storage/anat/sam/loci_Pip1_mapped.regions.bed", sep="\t", index=False, header=False)

