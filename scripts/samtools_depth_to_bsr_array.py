# %%
import os
from pandas.core.computation import align
import scipy.sparse as sparse
import pandas as pd
import numpy as np
from pathlib import Path
import pmbi.plotting as pmbip
import importlib
import matplotlib.pyplot as plt
from os import PathLike
from typing import Iterable
from pmbi.bio.aln.paf import read_paf_to_df
from pmbi.bio.aln.bam import read_bam_to_df, aligned_segment_get_gaps, cigar_get_ops, cigar_ref_consuming_ops, cigar_spec
import pmbi.bio.gtf as pgtf
from Bio import SeqIO 
import pysam
from enum import Enum

PREFIX=Path("/home/amsesk/super1")
FIGS=Path("/home/amsesk/figures/anat")

# %%
def fai_into_chr_sizes(fai_path: PathLike):
    chr_sizes = pd.read_csv(fai_path, sep="\t", header=None).iloc[:,0:2]
    chr_sizes.columns=["chr","size"]
    chr_sizes=chr_sizes.set_index("chr")
    return chr_sizes

# %%

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

# %% GBFF Annotation
# for record in SeqIO.parse("/home/amsesk/super1/anat/genome_ref/GCA_037367395.2_aLitPip1_p1.2_genomic.gbff", "genbank"):
#     feat_types = list()
#     for feat in record.features:
#         feat_types.append(feat.type)
#     print(pd.Series(feat_types).value_counts())

# %%

chr_sizes = fai_into_chr_sizes(PREFIX.joinpath("anat/genome_ref/GCA_037367395.2_aLitPip1_p1.2_genomic.fna.fai"))
# D_arrs = stdepth_into_coo_arrays(PREFIX.joinpath("anat/sam/AES103_Pip1.sorted.bam.depth"), chr_sizes)


# %%
gtf_path="/home/amsesk/super1/anat/genome_ref/GCF_037356315.1_aLitPip1_a1.2_genomic.gtf"
with open(gtf_path,'r') as h:
    gtf = pgtf.FeatureFile.from_gtf(h)

gtf[gtf["feature"]=="gene"]
print(gtf)
# %%
loci_bam = read_bam_to_df(PREFIX.joinpath("anat/sam/loci_Pip1_lenient.sorted.bam"))

# %%
loci_bam_af = pysam.AlignmentFile(PREFIX.joinpath("anat/sam/loci_Pip1_lenient.sorted.bam"), "rb")
qad = {}
for a in loci_bam_af:
    if a.query_name in qad:
        qad[a.query_name].append(a)
    else:
        qad[a.query_name] = [a]

list(qad.keys())
# %%

a=qad["annotation-ENSACAG00000000778:k676-1076"]

bs
# %%
class BlockAnnot(Enum):
    COMPLETELY_WITHIN=0
    STARTS_WITHIN=1
    ENDS_WITHIN=2
    NEAREST=3

# %%
def aligned_segment_gtf_overlaps(a, gtf, feature_type="gene"):
    gtf_sub = gtf[(gtf["seqname"]==a.reference_name) & (gtf["feature"]==feature_type)].reset_index(drop=True)
    start_sort = gtf_sub.sort_values("start")
    # end_sort = gtf_sub.sort_values("end")
    annots = []
    for bs,be in a.get_blocks():
        s_idx = np.searchsorted(a=start_sort["start"].to_numpy(), v=bs, side="right")
        # e_idx = np.searchsorted(a=end_sort["end"].to_numpy(), v=be, side="left")
        # start_sort_r = start_sort.iloc[s_idx,:][["start","end"]].to_numpy()
        # end_sort_r = end_sort.iloc[e_idx,:][["start","end"]].to_numpy()
        # print(start_sort)
        # print(s_idx)
        # print(bs,be)
        cand = start_sort.iloc[s_idx-1,:]
        cs,ce = (cand["start"], cand["end"])
        if s_idx>start_sort.shape[0]:
            next_cand = start_sort.iloc[s_idx,:]
        else:
            next_cand = None
        if bs>=cs:
            if bs<=ce:
                if be<=ce:
                    ret = (BlockAnnot(0), cand) 
                else:
                    ret = (BlockAnnot(1), cand)
            else: # bs>ce
                if next_cand is None:
                    ret = (BlockAnnot(3), cand)
                else:
                    ncs,_nce = (next_cand["start"], next_cand["end"])
                    if (bs-ce)<(ncs-bs):
                        ret = (BlockAnnot(3), cand)
                    elif (bs-ce)>(ncs-bs):
                        ret = (BlockAnnot(3), next_cand)
                    else:
                        raise ValueError("eqidistant??")
        else: # bs<cs
            if s_idx!=0:
                raise ValueError("bs<cs AND s_idx!=0 despite right-sided searchsorted... should not happen")
            else:
                first = start_sort.iloc[0,:]
                if be>=first["start"]:
                    ret = (BlockAnnot(2), first)
                else:
                    ret = (BlockAnnot(3), first)
        annots.append(ret) 
    return annots

aligned_segment_gtf_overlaps(qad["Esculentin-PeptideArea-Rchir"][0],gtf)
# %%
qad_overlaps = {}
for q,segs in qad.items():
    print(q)
    qad_overlaps[q] = []
    for s in segs:
        qad_overlaps[q].append(aligned_segment_gtf_overlaps(s, gtf))

# %%
k1=list(qad_overlaps.keys())
k2=list(qad.keys())
assert all([x==y for x,y in zip(k1,k2)])

# %%
ldict=[]
for qn in k1:
    for seg,overlaps in zip(qad[qn], qad_overlaps[qn]):
        for bid,(bs,be) in enumerate(seg.get_blocks()):
            overlap=overlaps[bid]
            ldict.append({
                "query": seg.query_name,
                "target": seg.reference_name,
                "block_id": bid,
                "block_start": bs,
                "block_end": be,
                "annot_relation": overlap[0].name,
                "gene_id": pgtf.FeatureFile._parse_attributes_string(overlap[1]["attributes"])["gene_id"],
                "annot_seqname": overlap[1]["seqname"],
                "annot_start": overlap[1]["start"],
                "annot_end": overlap[1]["end"]
            })

# %%
ov_df=pd.DataFrame(ldict)

# Test relations
np.where(ov_df.apply(lambda r: r["block_start"]<r["annot_start"] and r["block_end"]<=r["annot_end"] and r["block_end"]>=r["annot_start"], axis=1))
ov_df.iloc[1938,:]
assert ov_df[ov_df["annot_relation"]=="STARTS_WITHIN"].apply(lambda r: r["block_start"]>=r["annot_start"] and r["block_end"]>r["annot_end"], axis=1).all()
assert ov_df[ov_df["annot_relation"]=="COMPLETEY_WITHIN"].apply(lambda r: r["block_start"]>=r["annot_start"] and r["block_end"]<=r["annot_end"], axis=1).all()
assert ov_df[ov_df["annot_relation"]=="NEAREST"].apply(lambda r: (r["block_start"]<r["annot_start"] and r["block_end"]<r["annot_end"]) or (r["block_start"]>r["annot_start"] and r["block_end"]>r["annot_end"]), axis=1).all()

ov_df.to_csv("/home/amsesk/figures/anat/aln_block_gtf_gene_annotations.tsv", sep="\t", index=False)
ov_df[ov_df["block_id"]==0]["query"].value_counts().reset_index()["count"].sum()

ov_df.iloc[872,:]

# %%
qad_overlaps[k]
blks=a[0].get_blocks()
ovs =aligned_segment_gtf_overlaps(a[0], gtf)
# print(blks[0][0], blks[0][1])
print(ovs)
# %%


a[0].get_blocks()
a[0].get_blocks()
aligned_segment_get_gaps(a[0])

# %%

loci_bam
loci_bam.cigar_s
loci_bam.cigar_t
# loci_paf = read_paf_to_df(open(PREFIX.joinpath("anat/sam/loci_Pip1_lenient.paf"), 'r'))

bam_mapped = loci_bam[loci_bam["ref_start"]!=-1]
bam_mapped_loci = loci_bam[loci_bam["ref_start"]!=-1]["query"].unique()
loci_bam[loci_bam["ref_start"]!=-1]["query"].shape
loci_bam[loci_bam["ref_start"]!=-1]["query"].unique().shape
loci_bam["query"].unique().shape


# %% Write BED of locus alignment regions
loci_bam.iloc[[0,1761],:]
bam_to_bed(loci_bam.iloc[[0,1761],:], "/media/md0/webServer/data/amsesk/loci_Pip1_lenient_regions.bed")
loci_bam.columns
loci_bam.iloc[0,].cigar_t
loci_bam.iloc[0,].query
loci_bam

aln_df.iloc[1,:]["tags"]

# %% Print some statistics about loci mapping to reference genome
query_map_counts = bam_mapped["query"].value_counts().reset_index(drop=False).sort_values(["count", "query"])
query_map_counts["count"].sum()
query_map_counts[query_map_counts["count"]==1]

loci_bam.columns

qmap2 = query_map_counts[query_map_counts["count"]==2].iloc[0,0]
qmap1 = query_map_counts[query_map_counts["count"]==1]
qmapMulti = query_map_counts[query_map_counts["count"]>1]
qmapMulti.to_csv("/home/amsesk/figures/anat/loci_mapped_multi.tsv", sep="\t", index=False)
qmap1.shape
qmap1
-np.log10(1e-100)
loci_bam_mapped = loci_bam[loci_bam["ref_start"]!=-1]
# %%
mapq_as = pd.DataFrame({"mq": loci_bam_mapped["mq"], "as": loci_bam_mapped["tags"].apply(lambda ts: ts["AS"])})
panel = pmbip.Paneler(1,1,(3,3))
panel.next_ax().scatter(x=mapq_as["mq"], y=mapq_as["as"], marker=".")
pan

panel.fig.savefig("/home/amsesk/figures/anat/mapq_as_scatter.pdf")
# %%
loci_bam[loci_bam["query"]==qmap2]
loci_bam[loci_bam["query"]==qmap2].tags.tolist()
loci_paf[loci_paf["qname"]==qmap2]

pd.set_option('display.max_rows', None)
loci_bam_map1 = loci_bam[(loci_bam["query"].isin(qmap1["query"]))]
loci_bam_map1_mapped = loci_bam_map1[loci_bam_map1["ref_start"]!=-1]
loci_bam_map1_mapped.shape
loci_bam_map1_mapped[~loci_bam_map1_mapped["query"].str.startswith("annotation")]["query"].unique().shape
loci_bam_map1_mapped[loci_bam_map1_mapped["query"].str.startswith("annotation")]["query"].unique().shape

# %%
loci_paf[loci_paf["qname"].str.startswith("Brevinin")]

query_map_counts
query_map_counts.to_csv("/home/amsesk/figures/anat/locus_counts.csv", sep=",", index=False)
unmapped_loci = aln_df[aln_df["ref_start"]==-1]["query"].unique()
mapped_loci = aln_df[aln_df["ref_start"]!=-1]["query"].unique()
aln_df.iloc[aln_df[aln_df["query"].str.contains("Mhc")].index]["query"].value_counts()

aln_df.shape
query_map_counts[query_map_counts["count"]==1].shape

# %%
for mm in query_map_counts[query_map_counts["count"]==2]["query"]:
    mm_alns = aln_df[aln_df["query"]==mm]
    print(mm_alns)
    print(mm_alns["tags"].tolist())
    break

# %%
loci_paf["qname"].value_counts().reset_index(drop=False).sort_values(["count","qname"])
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

#############################
#############################
#############################

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

