
import gc
import gzip
import importlib
import os
import re
import itertools
import functools
from collections import namedtuple
from termios import N_MOUSE

import matplotlib.pyplot as plt
import numpy as np
import palettable
import pandas as pd
from Bio import SeqIO
from toolz import *

import pmbi.bio.aln.psl as PSL
import pmbi.bio.dna
import pmbi.plotting as pmbip

# %%
importlib.reload(PSL)
importlib.reload(pmbip)

# %% CHUNK: Change to working directory 
os.chdir("/storage/anat/")

# %% CHUNK: Read in the reference sequences
# tar = "s0004895-augustus-gene-0.9-mRNA-1:1-401"
# tar_seq = None
tars = {}
for rec in SeqIO.parse("loci_ref/Unique_1388_nov2019_noDegenerateNucl.fna", "fasta"):
    tars[rec.id] = str(rec.seq)

# %% CHUNK: Read in the query sequences
# que = "A00572:174:HGGKNDRXX:1:2101:3215:1078"
# que_seq = None
ques = {
    "R1": {},
    "R2": {},
}
# with gzip.open(os.path.join(QUE_PREFIX, "AES537_R1.fa.gz"), 'rb') as gzfa:
with gzip.open("fasta/AES103_R1.fa.gz", "rt") as fa:
    for rec in SeqIO.parse(fa, "fasta"):
        ques["R1"][rec.id] = str(rec.seq)

with gzip.open("fasta/AES103_R2.fa.gz", "rt") as fa:
    for rec in SeqIO.parse(fa, "fasta"):
        ques["R2"][rec.id] = str(rec.seq)

# # %% Some unused functions {{{
# def tostr(series):
#     return "".join(series.tolist())
#
# def insert(df, r, c, s):
#     ins_d = {k:v for k,v in zip(r,s)}
#     df.loc[r,c] = ins_d
#
# }}}
# %%

def _tName_in(r, l):
    if r.tName in l:
        return True
    else:
        return False

def _qName_in(r, l):
    if r.qName in l:
        return True
    else:
        return False

# %% CHUNK: Read in the PSL files
importlib.reload(PSL)
fs = [
    "psl/AES103_R1_blat.psl.gz",
    "psl/AES103_R2_blat.psl.gz",
]
ldict = []
importlib.reload(PSL)
pac = {}
with gzip.open(fs[0], "rt") as stream:
    pac["R1"] = PSL.PslAlignmentCollection.from_file(stream, target_seqs=tars, query_seqs=ques["R1"])

# with gzip.open(fs[1], "rt") as stream:
#     pac["R2"] = PSL.PslAlignmentCollection.from_file(stream, target_seqs=tars, query_seqs=ques["R2"])


# %%
this_one = "annotation-ENSXETG00000005933:1-401"
pa = pac["R1"]
pa_filt = pa.filter(lambda r: r.tName == this_one)


# %%
# n_blocks = []
# all_mm = []
# for r in pa_filt.rows:
#     blks = r._blocks()
#     mm = []
#     for blk in blks:
#         mm.append([x.targetPos for x in blk.mismatches(r.qSeq, r.tSeq)])
#     all_mm.append(mm)

# %%
nucl_to_int = {
    "A": 0,
    "T": 1,
    "C": 2,
    "G": 3,
    "N": 4,
    "-": 5
}
# TODO: Incorporate gaps, if needed
# TODO: Test with multiple blocks - currently testing on an alignment that only has one block
def blk_target_idx(blk, qSeq, tSeq, qidx):
    qSeq = np.array(list(qSeq))
    tSeq = np.array(list(tSeq))
    t_pos = pd.Index(np.arange(blk.tS, blk.tE))
    q_aln_range = np.arange(blk.qS, blk.qE)
    return pd.DataFrame({"t_pos": t_pos, f"q{qidx}": [nucl_to_int[n] for n in qSeq[q_aln_range]]})

def leading_target_idx(r, qidx):
    qSeq = np.array(list(r.qSeq))
    t_pos = np.arange( r.tStart-r.n_leading_query(), r.tStart)  
    q_range = np.arange(0, r.n_leading_query())
    return pd.DataFrame({"t_pos": t_pos, f"q{qidx}": [nucl_to_int[n] for n in qSeq[q_range]]})

def trailing_target_idx(r, qidx):
    qSeq = np.array(list(r.qSeq))
    t_pos = np.arange( r.tEnd, r.tEnd+r.n_trailing_query())  
    q_range = np.arange(r.qEnd, r.qEnd+r.n_trailing_query())
    return pd.DataFrame({"t_pos": t_pos, f"q{qidx}": [nucl_to_int[n] for n in qSeq[q_range]]})

[x.n_trailing_query() for x in pa_filt.rows]
# %%
alns = []
for qidx,row in enumerate(pa_filt.rows):
    la = leading_target_idx(row, qidx)
    ba = blk_target_idx(row._blocks()[0], row.qSeq, row.tSeq, qidx)
    ta = trailing_target_idx(row, qidx)
    fa = pd.concat([la, ba, ta], axis=0)
    assert all(fa["t_pos"].to_numpy() == np.arange(fa["t_pos"].min(), fa["t_pos"].max()+1))
    alns.append(fa)

# %%
full_aln = functools.reduce(lambda left,right: pd.merge(left, right, how="outer", on="t_pos"), alns).transpose()
full_aln.columns = np.array(full_aln.loc["t_pos"])
full_aln = full_aln.drop("t_pos", axis=0)

# %%
assert all(np.array(full_aln.columns) == np.arange(xticks.min(), xticks.max()+1))
xticks = np.array(full_aln.columns)
xticks = np.arange(xticks.min(), xticks.max()+1, step=5)


panel = pmbip.Paneler(nrow=1, ncol=1, figsize=(12,12))
pmbip.heatmap(matrix=full_aln, ax = panel.next_ax(), xlab = "pos", ylab = "query") 
# panel.current_ax.set_xticks(ticks=xticks, labels=xticks)
panel.fig.savefig("/storage/anat/figures/full_aln_hm.pdf")

# %%
pd.Series(n_blocks).value_counts()
# %%
    print(r)
    break

# %%
pa.rows[0]._aln_trailing()
# %%

# # pac["R1"].rows
# padf = pa.as_dataframe()
# padf.shape
# these_ones = pa.filter(_tName_in, [this_one]).select(["qName"])["qName"].tolist()
# padf[padf["qName"].isin(these_ones)]["tStarts"].astype(str).value_counts()[1:50]
# pa.filter(_qName_in, these_ones).select(["tName"]).value_counts()
# pa.select(["qName"])["qName"].value_counts()
#
# pa.rows[1]._aln_repr_blocks()
# pa.rows[2]._aln_repr()
# pa.rows[2].leading
# pa.rows[2].trailing
# pa.rows[2].show()
# pa.rows[0]._aln_repr()
# pa.rows[1]._aln_leading()
# pa.rows[1]._aln_trailing()
# pa.rows[0].tSize
# pa.rows[0].show()
# pac["R1"].filter(lambda
# %%
importlib.reload(PSL)
taR1 = PSL.TargetAlignments(pac=pac["R1"], target="annotation-ENSACAG00000000778:676-1076")
taR2 = PSL.TargetAlignments(pac=pac["R2"], target="annotation-ENSACAG00000000778:676-1076")
taR1.target_coverage()
taR2.target_coverage()
taR1_qname = [x.qName for x in taR1.rows]
taR2_qname = [x.qName for x in taR2.rows]
for qn in taR1_qname:
    if qn in taR2_qname:
        print(qn)

print([x for x in taR1.rows if x.qName == "A00572:174:HGGKNDRXX:2:2258:7075:36699"][0])
[x for x in taR1.rows if x.qName == "A00572:174:HGGKNDRXX:2:2258:7075:36699"][0].show()
print("----")
print([x for x in taR2.rows if x.qName == "A00572:174:HGGKNDRXX:2:2258:7075:36699"][0])
[x for x in taR2.rows if x.qName == "A00572:174:HGGKNDRXX:2:2258:7075:36699"][0].show()

[x.qName for x in taR1.rows]
[x.qName for x in taR2.rows]
len(taR2.rows)
ta.rows[0].show()
ta.aln

# %%
ta.target_coverage()
# %%
n_mm=list()
all_mm = list()
ra = list()
nucl2int = {
    "N": 0,
    "A": 1,
    "T": 2,
    "C": 3,
 "G": 4,
}
for a in ta.rows:
    mm = pipe(map(lambda b: b.mismatches(a.qSeq, a.tSeq), a._blocks()), concat, list)
    all_mm.append(mm)
    # mm_qBase = [x.queryBase for x in mm]
    # tar_pos = [x.targetPos for x in mm]
    # tar_arr = np.array([np.nan]*len(a.tSeq))
    # for ti in range(a.tStart,a.tEnd):
    #     tar_arr[ti] = nucl2int[a.tSeq[ti]]
    # for m in mm:
    #     tar_arr[m.targetPos] = nucl2int[m.queryBase]
    # ra.append(tar_arr)


# %%
TARS = pd.Series([x.tName for x in pa.rows]).value_counts()[0:20].index
TARS = pd.Index(["annotation-ENSACAG00000000778:676-1076", "annotation-ENSACAG00000028334:1095-1495"])
TARS
# with gzip.open(fs[0], "rt") as stream:
#     pac = PSL.PslAlignmentCollection.from_file(stream, target_seqs=tars, query_seqs=ques, n=None, targetSeq_include = TARS.tolist())

# %%
panel = pmbip.Paneler(4,5,figsize=(12,16))
for TAR in TARS:
    ta = PSL.TargetAlignments(pac=pac["R1"], target=TAR)
    all_mm = list()
    for a in ta.rows:
        mm = pipe(map(lambda b: b.mismatches(a.qSeq, a.tSeq), a._blocks()), concat, list)
        all_mm.append(mm)
    print(f"{TAR} - tcov")
    tcov = ta.target_coverage()
    print(f"{TAR} - mm_vc")
    mm_vc = pd.DataFrame(pipe(all_mm, concat, list))["targetPos"].value_counts()
    print(f"{TAR} - mm_plt_p_mm")
    # plt_p_mm = pd.DataFrame(index=range(0,len(ta.tSeq))).join(mm_vc).fillna(0).join(pd.DataFrame({"cov":tcov})).assign(prop_mm = lambda r: r["count"]/r["cov"])
    queryBaseCounts = pd.DataFrame(pipe(all_mm, concat, list)).groupby(["targetPos", "queryBase"]).size().reset_index(name="Count").pivot_table(values="Count", columns="queryBase", index="targetPos").reindex(range(0,len(ta.tSeq)))
    print(queryBaseCounts)
    mm_queryBase_prop = pd.DataFrame({"cov":tcov}).join(queryBaseCounts).reset_index(drop=False).rename(columns={"index": "targetPos"}).fillna(0).melt(id_vars=["targetPos", "cov"]).sort_values("targetPos").assign(prop_mm = lambda r: r["value"]/r["cov"])
    mm_queryBase_prop["prop_mm"] = mm_queryBase_prop["prop_mm"].fillna(0)
    print(mm_queryBase_prop)
    qB = mm_queryBase_prop["variable"].unique()
    mm_queryBase_prop_wide = mm_queryBase_prop.pivot_table(values = "prop_mm", index="targetPos", columns="variable")
    print(mm_queryBase_prop_wide)
    print(f"{TAR} - bar")
    ax = panel.next_ax()
    bot = pd.Series([0]*len(ta.tSeq))
    for b in qB:
        ax.bar(x=mm_queryBase_prop_wide.index, height=mm_queryBase_prop_wide[b], bottom=bot, width=1)
        bot = bot+mm_queryBase_prop_wide[b]

panel.fig.savefig("/storage/anat/figures/targetPos_hm.pdf")
# %%
tn=[]
for a in pac["R1"].rows: 
    tn.append(a.tName)

len(set(tn))
# %%
ta.rows[156]._aln_repr_blocks()
panel = pmbip.Paneler(1,1,figsize=(6,6))
panel.fig.savefig("/storage/anat/figures/targetPos_hm.pdf")

# %%
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

Z = linkage(X, method='ward') 
Z

# %%
# np.apply_along_axis(lambda x: np.where(np.isnan(x))
X = np.stack(ra)
alnmat = pd.DataFrame(X)
panel = pmbip.Paneler(1,2,figsize=(6,6))
# dn = dendrogram(Z, ax=panel.next_ax())
# pmbip.heatmap(alnmat.iloc[map(int, dn["ivl"]),], ax=panel.next_ax(), xlab="targetPos", ylab="aln")
pmbip.heatmap(alnmat, ax=panel.next_ax(), xlab="targetPos", ylab="aln")
# panel.next_ax().bar(x=list(range(0,len(sums))), height=sums, width=1)
panel.fig.savefig("/storage/anat/figures/targetPos_hm.pdf")

# %%






pipe([mm.queryBase for mm in concat(all_mm) if mm.targetPos==12], pd.Series).value_counts()
pipe(map(lambda mms: [mm.targetPos for mm in mms], all_mm), concat, pd.Series).value_counts()
apply(len, all_mm)
[mm.targetPos for mm in all_mm[2]]
len(pa_sub)
len(all_mm)
all_mm[0]

# %%
for idx,pa in enumerate(pa_sub):
    blks = pa._blocks()
    for b in blks:
        if b.q_insert_bases > 0 and b.t_insert_bases > 0:
            print(idx)

print(pa_sub[121404])
pa_sub[121404].show()
# %%
for b in concat(map(lambda pa: pa._blocks(), pa_sub)):
    if b.q_insert_bases > 0 and b.t_insert_bases > 0:
        print(b)

# %%
all_mm[2]
pa_sub[2].show()
pa_sub[2]._aln_repr_blocks()
pa_sub[2]._blocks()[0].mismatches(qSeq, t)


# pa_sub[2].show()

# n_mm

# %%
type(pa.rows[1])
pa.rows[1w
].inner["qSeq"]
type(tars)

a = pa.rows[1]
print(a)
a.qStarts
a.tStarts
a.blockSizes
a.tSeq
a.misMatches
prinjt(a._blocks()[1].mismatches(a.qSeq, a.tSeq))
print(a._blocks()[0].q_gapped(a.qSeq))
print(a._blocks()[0].t_gapped(a.tSeq))
a.show()
a._aln_repr_blocks()
a._fields()
# %%
blks = a._blocks()
qss = list(map(lambda b: b.q_gapped(a.qSeq), a._blocks()))
list(itertools.accumulate(qss))
list(qss)
list(itertools.accumulate(a._blocks(), lambda l,c: l.q_gapped(a.qSeq) + c.q_gapped(a.qSeq)))

# %%
ss = ""
for b in blks:
    ss += b.q_gapped(a.qSeq)

print(ss)



# %%
importlib.reload(PSL)
with gzip.open(fs[0], "rt") as stream:
    while not stream.readline().startswith("-"):
        continue
    for line in stream:
        print(line)
        p=PSL.PslAlignment.from_str(line)
        break

p.strand
# %%
# %% CHUNK: Functions that should be moved into one of the PSL-related classes {{{
def leading_query_r(aln):
    l = aln.leading_query_pos()
    print(l)
    lqr = np.array([min(l), max(l) + 1]) + (-1 * min(l))
    print(lqr)
    return aln.r.qSeq[lqr[0] : lqr[1]][::-1]


def leading_seqs(pa: PSL.PslAln, ques, tars):
    leaders = []
    for idx, row in enumerate(pa):
        aln = PSL.PslRowAln(row, ques, tars)
        if len(aln.leading_query_pos()) > 0:
            # aln.show()
            leaders.append(leading_query_r(aln))
    return leaders


def trailing_query_f(aln):
    t = aln.trailing_query_pos()
    tqr = np.array([aln.r.qSize - ((max(t) + 1) - min(t)), aln.r.qSize])
    return aln.r.qSeq[tqr[0] : tqr[1]][::1]


def trailing_seqs(pa: PSL.PslAln, ques, tars):
    trailers = []
    for idx, row in enumerate(pa):
        aln = PSL.PslRowAln(row, ques, tars)
        if len(aln.trailing_query_pos()) > 0:
            trailers.append(trailing_query_f(aln))
    return trailers


def count_nt_at_pos(seqs):
    longest = max([len(s) for s in seqs])
    ldict = []
    for i in range(0, longest):
        pos = []
        for s in seqs:
            if i >= len(s):
                continue
            else:
                pos.append(s[i])
        pos_count = pd.Series(pos).value_counts().to_dict()
        ldict.append(pos_count)
    counts = pd.DataFrame(ldict)
    counts[np.isnan(counts)] = 0
    return counts


# %%
def sequence_composition_bars(
    weights, ax, flip_x=False, show_legend=True, letter_order=["A", "T", "C", "G", "N"], pal = None
):
    if pal is None:
        pal = {k:v for k,v in zip(letter_order, ['red', 'green', 'blue', 'orange', 'purple'])}
    lens = [len(v) for k, v in weights.items()]
    assert len(set(lens)) == 1
    bottom = np.zeros(lens[0])
    for letter in letter_order:
        if letter not in weights:
            continue
        else:
            wt = weights[letter]
            ax.bar(np.arange(0, lens[0]), wt, 1, label=letter, bottom=bottom, color=pal[letter])
            bottom += wt
    if show_legend:
        ax.legend()
    if flip_x:
        ax.invert_xaxis()
    t = pmbip.Theme()
    t.axislabels.fontsize = 6
    t.apply_to(ax)
    return None


# }}}

# %% {{{
    
seqname = "annotation-ENSXP-018414552:2449-2849"
sub = [x for x in pa.rows if x.tName == seqname]
len(sub)
leading_w = count_nt_at_pos(leading_seqs(sub, ques, tars))
leading_w_p = leading_w.apply(lambda row: row/sum(row), axis=1)
trailing_w = count_nt_at_pos(trailing_seqs(sub, ques, tars))
trailing_w_p = trailing_w.apply(lambda row: row/sum(row), axis=1)

# %% 
def which_nucl(row):
    where_max = np.where(row==max(row))[0]
    if len(where_max) > 1:
        raise ValueError
    else:
        return row.index[where_max[0]]

lead_cons = ''.join(leading_w_p.apply(which_nucl, axis=1).to_list())
tail_cons = ''.join(trailing_w_p.apply(which_nucl, axis=1).to_list())

comb = lead_cons[::-1] + tars[seqname] + tail_cons
"/stor/home/ka33933/work/Run3/blat/AES537_R1_blat.psl.gz"
with open("/stor/home/ka33933/work/Run3/AES537_R1_blat_round1.fasta", 'w') as out:
    out.write(f">{seqname}__round1\n{comb}\n")

# %%
tars = {}
for rec in SeqIO.parse(
    "/stor/home/ka33933/work/Run3/realigned_round1/AES537_R1_blat_round1.fa", "fasta"
):
    tars[rec.id] = str(rec.seq)


    # %% CHUNK: Working with a round_1 realign file...
ldict = []
with gzip.open("/stor/home/ka33933/work/Run3/realigned_round1/AES537_R1_round1.psl.gz", "rt") as stream:
    pa = PSL.PslAln.from_file(stream, target_seqs=tars, query_seqs=ques)


aln = PSL.PslRowAln(pa.rows[len(pa.rows)-1], ques, tars)
aln.leading_query_pos()
aln.trailing_query_pos()
aln._aln_repr_blocks() 
aln.leading_query_pos()
aln.n_leading_query()
aln._aln_repr()
ls = leading_seqs(pa.rows, ques, tars)

print(aln.r)
ques[aln.r.qName]
aln._blocks()[0].q_gapped(qSeq = ques[aln.r.qName])
aln._blocks()[0].t_gapped(tSeq = tars['annotation-ENSXP-018414552:2449-2849__round1'])
# %%
nb = []
for idx, row in enumerate(pa.rows):
    aln = PSL.PslRowAln(row, ques, tars)
    nb.append(len(aln._blocks()))
    if len(aln._blocks()) > 1:
        print(aln._aln_repr_blocks())

pd.Series(nb).value_counts()
pd.Series(nb).value_counts().sum()

# %%
# }}}
# %% {{{
target_counts = (
    pd.DataFrame({"tar": [r.tName for r in pa.rows], "que": [r.qName for r in pa.rows]})
    .iloc[:, 0]
    .value_counts()
)
target_counts_sort = target_counts.sort_values(ascending=False)
target_counts
# %%

pal = {k:v for k,v in zip(["A","T","C","G","N"], palettable.wesanderson.FantasticFox2_5.mpl_colors)}
for i, tar in enumerate(target_counts_sort.iloc[0:20].index):
    panel = pmbip.Paneler(2, 2, figsize=(9, 3), format="png")
    pa_sub = [x for x in pa.rows if x.tName == tar]
    def plottable_weights(weights):
        w_a = {
            k: [vv for kk, vv in weights[k].items()]
            for k, v in weights.to_dict().items()
        }
        w_p = weights.apply(lambda r: r / r.sum(), axis=1)
        w_p = {k: [vv for kk, vv in w_p[k].items()] for k, v in w_p.to_dict().items()}
        return (w_a, w_p)
    # Leading
    leading_w = count_nt_at_pos(leading_seqs(pa_sub, ques, tars))
    leading_w_a, leading_w_p = plottable_weights(leading_w)
    print(panel.axs)
    sequence_composition_bars(leading_w_a, ax=panel.axs[0][0], flip_x=True, pal=pal)
    sequence_composition_bars(
        leading_w_p, ax=panel.axs[1][0], flip_x=True, show_legend=False, pal=pal
    )
    panel.axs[0][0].set_xlabel("pos")
    panel.axs[0][0].set_ylabel("count")
    panel.axs[1][0].set_xlabel("pos")
    panel.axs[1][0].set_ylabel("freq")
    # Trailing
    trailing_w = count_nt_at_pos(trailing_seqs(pa_sub, ques, tars))
    trailing_w_a, trailing_w_p = plottable_weights(trailing_w)
    sequence_composition_bars(trailing_w_a, panel.axs[0][1], pal=pal)
    sequence_composition_bars(trailing_w_p, panel.axs[1][1], show_legend=False, pal=pal)
    panel.axs[0][1].set_xlabel("pos")
    panel.axs[0][1].set_ylabel("count")
    panel.axs[1][1].set_xlabel("pos")
    panel.axs[1][1].set_ylabel("freq")
    panel.fig.suptitle(tar)
    panel.fig.savefig(f"/stor/home/ka33933/figures/bar_test_{i}.png")

# }}}

# %% TODO: Need to identify bad alignments somehow... for the purposes of making sure that the base proportions makes sense in some of the weirder edge cases (e.g., #19)
# target_counts_sort.index[0]
# pa.rows[0]
# target_counts_sort.index[0]
# pa_test = [x for x in pa.rows if x.tName == target_counts_sort.index[5]]
# pa_test = [PSL.PslRowAln(r, ques, tars) for r in pa_test]
# [(rr.r.qNumInsert, rr.r.qBaseInsert, rr.r.tNumInsert, rr.r.tBaseInsert) for rr in pa_test]
# pa_test[0].r.qEnd
# print(trailing_seqs(pa_test, ques, tars))
# # %%
# trailing_panel = sequence_composition_bars(count_nt_at_pos(trailing_seqs(pa_test, ques, tars)))
# trailing_panel.fig.savefig("/stor/home/ka33933/figures/trailing_bar_test.png")

# %% Some more unused code
# pa_test = [x for x in pa.rows if x.tName == "annotation-ENSXP-018414552:2449-2849"]
# leaders = leader_seqs(pa_test, ques, tars)
# weights = count_nt_at_pos(leaders)
# weights


# %%
# rs = []
# seq = tars["annotation-ENSXP-018414552:2449-2849"]
# for idx,row in enumerate(pa_test):
#     aln = PSL.PslRowAln(row, ques, tars)
#     rs += list(aln.w
# leading_query_pos()) + list(aln.trailing_query_pos())
#     if len(aln.trailing_query_pos()) > 0:
#         print(idx)
#         # leaders.append(leading_query_r(aln))

# %%
# aln = PSL.PslRowAln(pa_test[113034], ques, tars)
# aln.show()
# aln.trailing_query_pos()
# max(aln.trailing_query_pos())
# aln.r.qSeq[(150-73):150]
# trailing_query_f(aln)
# aln = PSL.PslRowAln(pa_test[92109], ques, tars)
# leading_query_r(aln)
# %%

# pos
# # %%
# aln = PSL.PslRowAln(pa_test[113040], ques, tars)
# aln.show()
# leading_query_r(aln)
# aln = PSL.PslRowAln(pa_test[92109], ques, tars)
# leading_query_r(aln)
#
# # %%
# print(aln.r)
# aln.r.qSeq[0:117]
# aln.show()
# aln._block_coords()
# aln.trailing_query_pos()
# aln.show()
#
# # %%
# rs = np.array(rs)
# nbins=abs(min(rs))+max(rs)
# nbins
# hist = np.histogram(rs, bins=nbins)
# np.array(list(hist[1])+[521])
#
# # %%
# t = pmbip.Theme()
# t.axislabels.fontsize=6
# panel = pmbip.Paneler(1,1,figsize=(9,3), format = "png")
# panel.next_ax().bar(x=hist[1], height=np.array(list(hist[0])+[0]))
# t.apply_to(panel.current_ax)
# panel.fig.savefig("/stor/home/ka33933/figures/overhang_hists.png")
# # %%
#
# # a=PslRowAln(pa_test[113421], ques, tars)
# # print(a.r)
# # a.show()
#
# hist
#
# >>>>>>> 2e0029999e68121b368a4b34060db000b226a28e
#
# print(i)
# # %%
#
# aln=PslRowAln(pa.rows[137], ques, tars)
# print(aln.r)
# blks = aln._blocks()
# aln.qSize_gapped
# aln.tSize_gapped
# aln.show()
# <<<<<<< HEAD
# aln.trailing_query_pos()
# =======
#
# np.histogram(np.array(list(aln.trailing_query_pos())+list(aln.trailing_query_pos())))
#
# >>>>>>> 2e0029999e68121b368a4b34060db000b226a28e
# sum([b.t_insert_bases for b in aln._blocks()])
# sum([b.q_insert_bases for b in aln._blocks()])
#
# aln.
# aln=PslRowAln(pa.rows[1267], ques, tars)
# aln=PslRowAln(pa.rows[137], ques, tars)
# print(aln.r)
# len(re.sub("[^ATCG]", "", aln._aln_repr()[0]))
# len(re.sub("[^ATCG]", "", aln._aln_repr()[1]))
# al
# aln._aln_track()
# aln.show(chunksize=150)
# aln.r.tBaseInsert
# aln.leading_query()
# aln.trailing_query()
# al
# aln._blocks()[8].reaches_end_of_query()
# aln._aln_repr_blocks()
# mm=[b.mismatches(qSeq=ques[aln.r.qName], tSeq=tars[aln.r.tName]) for b in aln._blocks()]
# mm=[b.mismatches(qSeq=aln.r.qSeq, tSeq=aln.r.tSeq) for b in aln._blocks()]
# mm
# len(mm)
# # %%
#
# print([b.q_has_insert for b in blks])
# print([b.t_has_insert for b in blks])
# print([b.q_insert_bases for b in blks])
# print([b.t_insert_bases for b in blks])
# # %%
# # # Calulate number bases inserted manually
# #     if block.qbS != block.psl_bs:
# #         assert block.psl_bs < block.qbS
# #         qbi += block.qbS - block.psl_bs
# #     if block.tbS != block.psl_bs:
# #         assert block.psl_bs < block.tbS
# #         tbi += block.tbS - block.psl_bs
#
#
# print(aln.r)
# qbi
# tbi
# # %%
#
# print(aln.r)
# aln._query_insert_point()
# aln.show()
# print(pa.rows[0])
# # %%
# n_mis = 0
# for qs,ts,bs in blks:
#
#     print(qs,ts,bs)
#
#
# # %%
# blks
# block_cov = []
# block_cov = [(qs,qs+bs) for qs,ts,bs in blks]
# inbetween = []
# last_end = 0
# for s,e in block_cov:
#     inbetween.append((last_end,s))
#     last_end = e
# inbetween.append((last_end,aln.r.qSize))
#
# block_cov
# inbetween
# # %%
# for i in range(0,10):
#     print(aln.r.qSeq[inbetween[i][0]:inbetween[i][1]])
#     print(aln.r.qSeq[block_cov[i][0]:block_cov[i][1]])
#
# # %%
# inbetween
# block_cov
# # %%
# for block in blks:
#     qs,ts,bs = block
#     for i in list(range(qs,qs+bs)):
#         block_cov.append(i)
#
# [x for x in range(0,len(qstr)) if x not in block_cov]
# # %%
# block_cov
# pd.Series(aln.r.tStarts).diff(-1)
# d = pd.concat([pd.Series(aln.r.tStarts).diff(-1), pd.Series(aln.r.qStarts).diff(-1)], axis=1)
# d.columns = ["tsd", "qsd"]
# d["t_diff"] = d["tsd"]-d["qsd"]
# pd.concat([pd.DataFrame(blks), d], axis=1)
# d
# tstr =list(aln.r.tSeq)
# qstr =list(aln.r.qSeq)
# tEnd =  aln.r.tEnd
# qEnd =  aln.r.qEnd
# for bidx,block in enumerate(aln._block_coords()):
#     print(block)
#     qs,ts,bs = block
#     tb=aln.r.tSeq[ts:ts+bs]
#     qb=aln.r.qSeq[qs:qs+bs]
#     for i in range(0,bs):
#         if tb[i] != qb[i]:
#             n_mis+=1
#     t_diff= d.iloc[bidx,:]["t_diff"]
#     if not pd.isna(t_diff):
#         if t_diff > 0:
#             inc = int(abs(t_diff))
#             tstr = tstr[0:ts+bs] + ["-"]*inc + tstr[ts+bs:]
#             tEnd+=inc
#         elif t_diff < 0:
#             inc = int(abs(t_diff))
#             qstr = qstr[0:qs+bs] + ["-"]*inc + qstr[qs+bs:]
#             qEnd+=inc
#         else:
#             pass
#     # print(tb,qb,t_diff)
# print(f"{''.join(tstr[aln.r.tStart:tEnd])}\n{''.join(qstr[aln.r.qStart:qEnd])}")
#
# blks
# n_mis
# # %%
# qa = list(aln.r.qSeq[aln.r.qStart:aln.r.qEnd])
# ta = list(aln.r.tSeq[aln.r.tStart:aln.r.tEnd])
# qa
# ta
#
# new = []
# for i in range(0,len(qa)):
#     if qa[i] != ta[i]:
#         qa.insert(i, "-")
#     else:
#         pass
#
# print(qa)
#
# # %%
# aln.r.tSeq[280:378]
# pmbi.bio.dna.revcomp(aln.r.qSeq)
# ques[pa.rows[0].qName]
# aln.df.loc[range(0,400),:]
# # %%
# a = pd.DataFrame({
#     "aln": [" "]*1200,
#     "tseq": [" "]*1200,
#     "qseq": [" "]*1200
#     })
# a.index = range(-400,800)
#
#
# insert(a, range(0,r.tSize), "tseq", list(tars[r.tName]))
# insert(a, t_range_q, "qseq", list(ques[r.qName]))
#
#
#
# a.loc[t_range_q,"qseq"] = {k:v for k,v in zip(t_range_q, list(ques[r.qName]))}
# a.loc[t_range_q,:]
#
#
#
# # %%
#
# # %%
# def target_aln_repr(p):
#     s = list("."*(p.tSize))
#     for i,ts in enumerate(p.tStarts):
#         for x in range(ts, (ts)+p.blockSizes[i]):
#             s[x] = '|'
#     chunk_size = 100
#     breaks = list(range(0, p.tSize, chunk_size))
#     return [''.join(s[i:i+chunk_size:]) for i in breaks]
#
#
# print(pa.rows[1])
# 187-243
# 151-95
# pa.rows[0].tStarts
# pa.rows[0].tStart
# pa.rows[0].qStart
# print(dir(pa.rows[0]))
# def aln_repr(p):
#     s = list("."*(p.tSize))
#     q = list(" "*(p.tSize))
#     head_hang = 0
#     if p.tStart-p.qStart<0:
#         head_hang = abs(p.tStart-p.qStart)
#     t_pos = range(0-head_hang, p.tSize)
#     q_pos = range(p.qStart-p.tStart, p.tSize)
#     # print(t_pos)
#     # print(len(t_pos), len(q_pos))
#     # print(t_pos, s, q_pos, q)
#     print(pd.Series(t_pos), pd.Series(q_pos))
#     df = pd.DataFrame({
#         "t_pos": pd.Series(t_pos),
#         "t_repr": s,
#         "q_pos": pd.Series(q_pos),
#         "q_repr": q
#         })
#     # print(df)
#     df.loc[np.where(df.q_pos>p.qSize)[0], "q_pos"] = np.inf
#     print(df)
#     # if head_hang < 0:
#     #     head = pd.DataFrame({"t_repr": " "*abs(head_hang), "q_repr": "Q"*abs(head_hang)})
#     #     df = pd.concat([head, df], axis=1)
#     return(df)
#     s_head = []
#     s_tail = []
#     q_head = []
#     q_tail = []
#     chunk_size = 100
#     print(p.qStart, p.qEnd)
#     if p.qStart < p.qEnd:
#         for i,ts in enumerate(p.tStarts):
#             print(ts-p.qStarts[i], (ts)+p.blockSizes[i])
#             for x in range(ts-p.qStarts[i], (ts)+p.blockSizes[i]):
#                 if x<0:
#                     s_head = [" "]+s_head
#                     q_head = ["Q"]+q_head
#                 else:
#                     if x <= (p.tSize-1):
#                         if q[x] != ' ':
#                             add = 'q'
#                         else:
#                             add = 'Q'
#                         if x>=ts:
#                             s[x] = '|'
#                         q[x] = add
#                     else:
#                         s_tail = s_tail+[" "]
#                         q_tail = q_tail+["Q"]
#         s = s_head + s + s_tail
#         q = q_head + q + q_tail
#         breaks = list(range(0, len(s), chunk_size))
#         return ( [''.join(s[i:i+chunk_size:]) for i in breaks], [''.join(q[i:i+chunk_size:]) for i in breaks] )
#     else:
#         return (['N','N'], ['N','N'])
#
# # %%
# # p.tSize-p.tStarts[0]
# n_gt1_block = 0
# # for p in pa.rows[9991:9992]:
# # for p in pa.rows[9980:10000]:
# for p in pa.rows:
#     if p.blockCount>1:
#         n_gt1_block+=1
#     if isinstance(p.qBaseInsert, list):
#         print(p)
#         break
#     # print('X'*(p.qStart-1) + '|'*(p.qSize-p.qStart))
#     print(
#           # p.aln_len()["query"],
#           # p.aln_len()["target"],
#           p.qStart,
#           p.qEnd,
#           p.tStart,
#           p.tEnd,
#           p.blockSizes,
#           # p.blockCount,
#           [(t,q) for t,q in zip(p.qStarts,p.tStarts)],
#           # p.tStarts,
#           # p.qNumInsert,
#           # p.qBaseInsert,
#           # p.tNumInsert,
#           # p.tBaseInsert,
#           # p.strand,
#           p.dist_from_target_ends(),
#           p.qSize,
#           p.tSize
#           )
#     # t,q = aln_repr(p)
#     df = aln_repr(p)
#     # print(df)
#     # for i in range(0,len(t)):
#     #     print(t[i])
#     #     print(q[i])
#     #     print('')
#     # print("\n\\\\\\\\\\\\\\\\\\\\\\\ \n")
#     # break
#
# # %%
# n_gt1_block
# len(pa.rows)
#
# print(len(pa.rows))
# pa.rows[0]
# df_r1 = pd.DataFrame(ldict)
#
# # %%
# class PslRow(object):
#     def __init__(self,
#                  row,
#                 ):
#         self.target_id=row.t_id
#         self.target_len=row.t_len
#         self.target_start=row.t_start
#         self.target_end=row.t_end
#         self.query_id=row.q_id
#         self.query_start=row.q_start
#         self.query_end=row.q_end
#         self.mismatch=row.mismatch
#         self.length=row.length
#     def __repr__(self):
#         return f"{self.query_id}--{self.target_id}"
#     def aln_len(self):
#         return {
#                 "query": abs(self.query_end - self.query_start),
#                 "target": abs(self.target_end - self.target_start),
#                 }
#     def ngaps(self):
#         aln_len = self.aln_len()
#         return {
#                 "query": self.length - aln_len["query"],
#                 "target": self.length - aln_len["target"],
#                 }
#     def query_orient(self):
#         if self.query_start > self.query_end:
#             return "-"
#         elif self.query_start < self.query_end:
#             return "+"
#         else:
#             return "0"
#
# # %%
# df_r1.shape
# for i,r in enumerate(df_r1.itertuples()):
#     p = PslRow(r)
#     al = p.aln_len()
#     if al["query"] != al["target"]:
#         break
#     if i == df_r1.shape[0]-1:
#         print("END")
#
#
# # %%
# df_r1
# df_r1.iloc[684996,:]
# p.ngaps()
# p.target_id
# df_r1
# # %%
# df_r1.apply(lambda r: PslRow(r), axis=0)
# df_r1
# df_r1["q_aligned_len"] = abs(df_r1["q_end"] - df_r1["q_start"])
# df_r1["t_aligned_len"] = abs(df_r1["t_end"] - df_r1["t_start"])
# df_r1["q_gaps"] = df_r1["length"] - df_r1["q_aligned_len"]
# df_r1["t_gaps"] = df_r1["length"] - df_r1["t_aligned_len"]
# df_r1[df_r1["q_gaps"]==df_r1["t_gaps"]]
# df_r1[df_r1["q_aligned_len"]!=df_r1["t_aligned_len"]]
#
# # %%
# aln_rank = df.groupby("t_id").size().sort_values(ascending=False)
# aln_rank
# n_alignments_by_target = np.log10(aln_rank.to_numpy())
# most_mapped = df[df.t_id == aln_rank.index[0]]
#
# aln_rank_q = df.groupby("q_id").size().sort_values(ascending=False)
# aln_rank_q.to_numpy()
# n_alignments_by_query = np.log10(aln_rank_q.to_numpy())
# most_mapped_q = df[df.t_id == aln_rank_q.index[0]]
#
# # %%
# def label(ax, xlab, ylab):
#     ax.set_xlabel(xlab)
#     ax.set_ylabel(ylab)
#     return ax
#
# # %%
# # def overhang_reads(df, cutoff, strand = "+"):
# #     if strand == "+":
# #         correct_orientation = df[df.q_start<df.q_end]
# #         within = correct_orientation[correct_orientation.t_start>=(correct_orientation.t_len-cutoff)]
# #         return within
# #     elif strand == "-":
# #         correct_orientation = df[df.q_start>df.q_end]
# #         within = correct_orientation[correct_orientation.t_end<=cutoff]
# #         return within
#
# # %%
# def overhang_reads_unstranded(df, cutoff):
#     left = df[(df.t_end<=cutoff) & (df.q_start>df.q_end)]
#     right = df[(df.t_start>=(df.t_len-cutoff)) & (df.q_start<df.q_end)]
#     return right,left
#
# # %%
# right,left = overhang_reads_unstranded(df[df.t_id==examp], 150)
# left["q_max_extent"] = left.t_end - (left.q_len-(left.q_len-left.q_start)) - left.q_gaps
# right["q_max_extent"] = right.t_start + (right.q_len-right.q_start) + right.q_gaps
# right["q_max_extent"].min()
#
# left["q_max_extent"].max()
# left[left.q_max_extent>100]
#
# # %%
# colors = palettable.cartocolors.qualitative.Antique_10.mpl_colors
# colors
# # %%
# t = pmbip.Theme()
# t.inner
# t.axislabels.fontsize=6
# panel = pmbip.Paneler(2,3,figsize=(9,3), format = "png")
# panel.next_ax().hist(n_alignments_by_target, 100)
# panel.current_ax = label(panel.current_ax, "log10 n_alignments", "count")
# t.apply_to(panel.current_ax)
#
# t=pmbip.Theme()
# t.title.text = aln_rank.index[0]
# t.title.fontdict = {"fontsize":8}
# t.axislabels.fontsize=6
#
# panel.next_ax().hist(most_mapped["t_start"], 50, color=colors[0])
# panel.current_ax = label(panel.current_ax, "position", "count")
# t.apply_to(panel.current_ax)
#
# panel.next_ax().hist(most_mapped["t_end"], 50, color=colors[1])
# panel.current_ax = label(panel.current_ax, "position", "count")
# t.apply_to(panel.current_ax)
#
# t = pmbip.Theme()
# t.inner
# t.axislabels.fontsize=6
# panel.next_ax().hist(n_alignments_by_query, 100)
# panel.current_ax = label(panel.current_ax, "log10 n_alignments", "count")
# t.apply_to(panel.current_ax)
#
#
# t = pmbip.Theme()
#
# panel.next_ax().hist(left["q_max_extent"], 50, color=colors[3])
# panel.current_ax.hist(right["q_max_extent"], 50, color=colors[4])
# # panel.current_ax.vlines(x=[0,400], ymin=0, ymax=1500)
# panel.current_ax = label(panel.current_ax, "position", "count")
# t.apply_to(panel.current_ax)
#
# panel.fig.savefig("/stor/home/ka33933/figures/blat_hists.png")
# # %%
#
