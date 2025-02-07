# %%
from collections import Counter

import numpy as np
import pandas as pd
import pysam

from pmbi.bio.aln import samflag_bin_repr
from pmbi.plotting import Paneler

# %%
# d=pd.read_csv("/home/amsesk/super1/cdiff_evo/bam/Ancestor.Day0_depth.txt", sep="\t", compression= "gzip", header=None)
# d.columns = ["chrom", "pos", "depth"]
# d.chrom.unique()
#
# d[d.chrom == "NZ_JACGTL010000002.1"].depth.mean()

valid_single_flags = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
# %%
rbam = pysam.AlignmentFile(
    "/home/amsesk/super1/cdiff_evo/bam/Ancestor.Day0_trimmed.bam", "rb"
)
fl = []
for _i, a in enumerate(rbam):
    fl.append(a.flag)

len(fl)
# %%
fls = pd.DataFrame(pd.Series(fl).value_counts()).reset_index()
fls.columns = ["flag", "count"]

# %%
fls["flag_repr"] = fls["flag"].apply(lambda val: samflag_bin_repr(val))
fls["aligned"] = fls["flag_repr"].apply(lambda val: int(val[12 - 3]) == 0)
len(fl) - fls[fls.aligned]["count"].sum()

(20000000 * 150) / 4400000

# %%

coords = pd.read_csv(
    "/home/amsesk/super1/cdiff_evo/bam/Ancestor.Day0_depth.txt",
    sep="\t",
    compression="gzip",
    header=None,
)
coords.columns = ["chrom", "pos", "depth"]

coords[coords.chrom == chrom]
# %%
from matplotlib.scale import FuncScale


# def log1p_forward(values):
#     return np.log(values+1)
#
# def log1p_reverse(values):
#     return np.e**(values+1)-1
#
# FuncScale(panel.current_ax, functions = (log1p_forward, log1p_reverse)).get_transform()


chroms = coords.chrom.unique()
panel = Paneler(len(chroms), 1, (10, 10))
for chrom in chroms:
    these_coords = coords[coords.chrom==chrom]
    x = these_coords.pos.to_numpy()
    y = these_coords.depth.to_numpy()
    y = np.log10(y+1)
    panel.next_ax().plot(x,y, linewidth = 0.2)

panel.fig.savefig("/home/amsesk/figures/cdiff_evo/depth.png")

# %%
y
np.sort(np.array([3,4,5,2,10,1]))

# %%

