# %%
import scipy.sparse as sparse
import pandas as pd
import numpy as np
from pathlib import Path
import pmbi.plotting as pmbip
import importlib
import matplotlib.pyplot as plt
from os import PathLike
from typing import Iterable
import pysam

PREFIX=Path("/home/amsesk/super1")
FIGS=Path("/home/amsesk/figures/anat")

# %%
def fai_into_chr_sizes(fai_path: PathLike):
    chr_sizes = pd.read_csv(fai_path, sep="\t", header=None).iloc[:,0:2]
    chr_sizes.columns=["chr","size"]
    chr_sizes=chr_sizes.set_index("chr")
    return chr_sizes

def stdepth_into_coo_arrays(path: PathLike, chr_sizes: pd.DataFrame) -> Iterable[sparse.sparray]:
    dep = pd.read_csv(path, sep="\t", header=None, names=["chr", "pos", "depth"])
    for this_chr in dep["chr"].unique():
        print(f"Processing: {this_chr}")
        this_chr_size = chr_sizes.loc[this_chr,"size"]
        sub = dep[dep["chr"]==this_chr]
        depth_arr = sub["depth"].to_numpy()
        cols_arr = np.full(len(sub["pos"]), 0)
        pos_arr = sub["pos"].to_numpy()
        yield sparse.coo_array((depth_arr, (pos_arr, cols_arr)), shape=(this_chr_size,1))

def read_bam(bam_path: PathLike) -> pysam.AlignmentFile:
    """Read a BAM file using pysam.

    Args:
        bam_path: Path to the BAM file

    Returns:
        pysam.AlignmentFile object for reading alignments
    """
    return pysam.AlignmentFile(bam_path, "rb")

# %%

chr_sizes = fai_into_chr_sizes(PREFIX.joinpath("anat/genome_ref/GCA_037367395.2_aLitPip1_p1.2_genomic.fna.fai"))
# D_arrs = stdepth_into_coo_arrays(PREFIX.joinpath("anat/sam/AES103_Pip1.sorted.bam.depth"), chr_sizes)


# %%
loci_bam = read_bam(PREFIX.joinpath("anat/sam/loci_Pip1.sorted.bam"))
ldict = []
for aln in loci_bam:
    row = {
        "query": aln.query_name, 
        "reference": aln.reference_name,
        "ref_start": aln.reference_start, 
        "ref_end": aln.reference_end
    }
    ldict.append(row)

aln_df = pd.DataFrame(ldict)

# %% Write BED of locus alignment regions
aln_df[["reference", "ref_start", "ref_end", "query"]].to_csv("/media/md0/webServer/data/amsesk/loci_Pip1_regions.bed", sep="\t", index=False, header=False)

# %% Print some statistics about loci mapping to reference genome
aln_df["query"].value_counts().reset_index(drop=False).sort_values(["count", "query"]).to_csv("/home/amsesk/figures/anat/locus_counts.csv", sep=",", index=True)
unmapped_loci = aln_df[aln_df["ref_start"]==-1]["query"].unique()
mapped_loci = aln_df[aln_df["ref_start"]!=-1]["query"].unique()
aln_df.iloc[aln_df[aln_df["query"].str.contains("Mhc")].index]["query"].value_counts()


print(unmapped_loci.shape)
print(mapped_loci.shape)


np.where(pd.Series(unmapped_loci).str.startswith("annotation-"))[0].shape
np.where(pd.Series(mapped_loci).str.startswith("annotation-"))[0].shape

aln_df[aln_df["ref_start"]!=-1]["query"].value_counts().value_counts()
aln_df[(aln_df["ref_start"]!=-1) & (aln_df["query"].str.startswith("annotation-"))]["query"].value_counts().value_counts()
aln_df[(aln_df["ref_start"]!=-1) & (~aln_df["query"].str.startswith("annotation-"))]["query"].value_counts().value_counts()

# %%

# %%
importlib.reload(pmbip)
singular_mapped_loci = aln_df.iloc[np.where(aln_df[aln_df["ref_start"]!=-1]["query"].value_counts()==1)[0],:]
fs=list(Path(PREFIX.joinpath("anat/sam")).iterdir())
fs.sort()
nucl_order = ["A", "C", "G", "T"]
nucl_colors = ["blue", "yellow", "black", "red"]
covs = {}

try:
    panel.close()
except:
    pass

for path in fs:
    if path.name.startswith("23007FL") and path.name.endswith(".sorted.bam"):
        panel = pmbip.Paneler(10,5,(16,16),format="pdf", output_prefix=FIGS.joinpath(f"{path.name}_cov"))
        plt.close("all")
        plt.clf()
        print(f"Reading BAM: {path}")
        sample_bam = read_bam(path)
        covs[path.name] = {}
        print(f"Counting coverages...")
        for _idx, query, ref, ref_start, ref_end in singular_mapped_loci.itertuples():
            ref_start, ref_end = (int(ref_start), int(ref_end))
            print(f"Working on query #{_idx}: {query}")
            counts_df = pd.DataFrame({n:c for n,c in zip(nucl_order, sample_bam.count_coverage(contig=ref, start=ref_start, stop=ref_end))}) 
            counts_df["pos"] = range(ref_start, ref_end)
            ax = panel.next_ax()
            for n,col in zip(nucl_order, nucl_colors):
                ax.plot(counts_df["pos"], counts_df[n], c=col, linewidth=0.1)
                ax.set_title(query, fontdict={"fontsize":4.0})
            covs[path.name][query] = counts_df
        panel.close()

# %% Average coverage in locus regions
mean_covs = []
for f,d in covs.items():
    these_means = []
    for q,df in d.items():
        these_means.append(df.drop(columns="pos").apply(sum, axis=1).mean())
    mean_covs.append(pd.Series(these_means, index=d.keys()))

# %%

D_arr = D.to_array()
panel = Paneler(1,1,(8,8))
panel.next_ax().scatter(sub["pos"].values, sub["depth"].values, s=0.75, marker=".", edgecolors="none")
panel.fig.savefig(FIGS.joinpath(f"AES103_Pip1_depth/AES103_Pip1_depth_{this_chr}.png"), dpi=1500)

# %%
r = np.array([200e6, 250e6]).astype(np.int64)
panel = Paneler(1,1,(12,12))
# panel.next_ax().plot(range(0,10000), D_arr[0:10000,0], linewidth=0.1)
panel.next_ax().scatter(range(r[0], r[1]), D_arr[r[0]:r[1],0], s=0.5, marker=".")
panel.fig.savefig(FIGS.joinpath("AES103_Pip1_depth.png"))
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
loci_sam = pd.read_csv(, sep="\t", header=None, usecols=[0,2], names=["query","reference"])
loci_sam["count"] = 1
loci_mat = loci_sam.pivot_table(index="query", values="count", columns="reference", aggfunc=sum).fillna(0).drop(columns="*")
loci_mat.to_csv(PREFIX.joinpath("anat/loci_Pip1_chrom_mat.csv"), sep=",")
                       
# %%
loci_sam = pd.read_csv("/storage/anat/sam/loci_Pip1.sorted.sam", sep="\t", header=None, usecols=[0,2,3,9], names=["query","reference","pos","seq"])
loci_sam=loci_sam[loci_sam["reference"]!="*"]
loci_sam["seqlen"] = loci_sam["seq"].apply(len)
loci_sam["start"] = loci_sam["pos"]-loci_sam["seqlen"]
loci_sam["stop"] = loci_sam["pos"]+loci_sam["seqlen"]
loci_sam[["reference", "start", "stop", "query"]].to_csv("/storage/anat/sam/loci_Pip1_mapped.regions.bed", sep="\t", index=False, header=False)

