import importlib
import os
import re
import itertools
from collections import namedtuple
from termios import N_MOUSE

import matplotlib.pyplot as plt
import numpy as np
import palettable
import pandas as pd
from Bio.Align.psl import AlignmentIterator
from Bio import SeqIO

import pmbi.plotting as pmbip
import pmbi.bio.dna
import pmbi.bio.aln.psl as PSL

importlib.reload(PSL)
importlib.reload(pmbip)

# %% PSL fields {{{
# psl_fields = ["matches",
#               "misMatches",
#               "repMatches",
#               "nCount",
#               "qNumInsert",
#               "qBaseInsert",
#               "tNumInsert",
#               "tBaseInsert",
#               "strand",
#               "qName",
#               "qSize",
#               "qStart",
#               "qEnd",
#               "tName",
#               "tSize",
#               "tStart",
#               "tEnd",
#               "blockCount",
#               "blockSizes",
#               "qStarts",
#               "tStarts"]
# }}}

# %% 
REF_PREFIX="/home/amsesk/penn-microbioinfo/anat/ref/"
QUE_PREFIX="/home/amsesk/penn-microbioinfo/anat/psl"

tar = "s0004895-augustus-gene-0.9-mRNA-1:1-401"
tar_seq = None
tars = {}
for rec in SeqIO.parse(os.path.join(REF_PREFIX, "Unique_1388_nov2019_noDegenerateNucl.fna"), "fasta"):
    tars[rec.id]=str(rec.seq)

que = "A00572:174:HGGKNDRXX:1:2101:3215:1078"
que_seq = None
ques = {}
for rec in SeqIO.parse(os.path.join(QUE_PREFIX, "AES537_R1.fa"), "fasta"):
    ques[rec.id]=str(rec.seq)

# # %% {{{
# def tostr(series):
#     return "".join(series.tolist())
#
# def insert(df, r, c, s):
#     ins_d = {k:v for k,v in zip(r,s)}
#     df.loc[r,c] = ins_d
#
# }}}



# %%
fs = ["/stor/home/ka33933/work/Run3/blat/AES537_R1_blat.psl", "/stor/home/ka33933/work/Run3/blat/AES537_R2_blat.psl"]
ldict = []
with open(fs[0], 'r') as stream:
    pa = PSL.PslAln.from_file(stream, target_seqs=tars, query_seqs=ques)

# %%
def leading_query_r(aln):
    l = aln.leading_query_pos()
    lqr=np.array([min(l), max(l)+1])+(-1*min(l))
    return aln.r.qSeq[lqr[0]:lqr[1]][::-1]

def leading_seqs(pa: PslAln, ques, tars):
    leaders = []
    for idx,row in enumerate(pa):
        aln = PSL.PslRowAln(row, ques, tars)
        if len(aln.leading_query_pos()) > 0:
            leaders.append(leading_query_r(aln))
    return leaders

def trailing_query_f(aln):
    t = aln.trailing_query_pos()
    tqr = np.array([aln.r.qSize-((max(t)+1)-min(t)), aln.r.qSize])
    return aln.r.qSeq[tqr[0]:tqr[1]][::1]

def trailing_seqs(pa: PslAln, ques, tars):
    trailers = []
    for idx,row in enumerate(pa):
        aln = PSL.PslRowAln(row, ques, tars)
        if len(aln.trailing_query_pos()) > 0:
            trailers.append(trailing_query_f(aln))
    return trailers

def count_nt_at_pos(seqs):
    longest = max([len(s) for s in seqs])
    ldict = []
    for i in range(0,longest):
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

def sequence_composition_bars(weights):
    xrange = weights.index.to_list()
    weights_a = {k: [vv for kk,vv in weights[k].items()] for k,v in weights.to_dict().items()}
    weights_p = weights.apply(lambda r: r/r.sum(), axis=1)
    weights_p = {k: [vv for kk,vv in weights_p[k].items()] for k,v in weights_p.to_dict().items()}
    t = pmbip.Theme()
    t.axislabels.fontsize=6
    panel = pmbip.Paneler(2,1,figsize=(9,3), format = "png")
    ax = panel.next_ax()
    bottom = np.zeros(len(xrange))
    for nt,c in weights_a.items():
        ax.bar(xrange, c , 1, label = nt, bottom = bottom)
        bottom += c
    t.apply_to(panel.current_ax)
    ax = panel.next_ax()
    bottom = np.zeros(len(xrange))
    for nt,c in weights_p.items():
        ax.bar(xrange, c , 1, label = nt, bottom = bottom)
        bottom += c
    t.apply_to(panel.current_ax)
    return panel

# %%
pa_test = [x for x in pa.rows if x.tName == "annotation-ENSXP-018414552:2449-2849"]
leaders = leader_seqs(pa_test, ques, tars)
weights = count_nt_at_pos(leaders)
weights


# %% 
rs = []
seq = tars["annotation-ENSXP-018414552:2449-2849"]
for idx,row in enumerate(pa_test):
    aln = PSL.PslRowAln(row, ques, tars)
    rs += list(aln.leading_query_pos()) + list(aln.trailing_query_pos())
    if len(aln.trailing_query_pos()) > 0:
        print(idx)
        # leaders.append(leading_query_r(aln))

# %%
aln = PSL.PslRowAln(pa_test[113034], ques, tars)
aln.show()
aln.trailing_query_pos()
max(aln.trailing_query_pos())
aln.r.qSeq[(150-73):150]
trailing_query_f(aln)
aln = PSL.PslRowAln(pa_test[92109], ques, tars)
leading_query_r(aln)
# %%
target_counts = pd.DataFrame({"tar": [r.tName for r in pa.rows], "que": [r.qName for r in pa.rows]}).iloc[:,0].value_counts()
target_counts_sort = target_counts.sort_values(ascending=False)
# %%
for i,tar in enumerate(target_counts_sort.iloc[0:20].index):
    pa_sub = [x for x in pa.rows if x.tName == tar]
    leading_panel = sequence_composition_bars(count_nt_at_pos(leading_seqs(pa_sub, ques, tars)))
    trailing_panel = sequence_composition_bars(count_nt_at_pos(trailing_seqs(pa_sub, ques, tars)))
    leading_panel.fig.savefig(f"/stor/home/ka33933/figures/leading_bar_test_{i}.png")
    trailing_panel.fig.savefig(f"/stor/home/ka33933/figures/trailing_bar_test_{i}.png")
    print(tar)


# %% TODO: Need to identify bad alignments somehow... for the purposes of making sure that the base proportions makes sense in some of the weirder edge cases (e.g., #19)
target_counts_sort.index[0]
pa.rows[0]
target_counts_sort.index[0]
pa_test = [x for x in pa.rows if x.tName == target_counts_sort.index[5]]
pa_test = [PSL.PslRowAln(r, ques, tars) for r in pa_test]
[(rr.r.qNumInsert, rr.r.qBaseInsert, rr.r.tNumInsert, rr.r.tBaseInsert) for rr in pa_test]
pa_test[0].r.qEnd
print(trailing_seqs(pa_test, ques, tars))
# %%
trailing_panel = sequence_composition_bars(count_nt_at_pos(trailing_seqs(pa_test, ques, tars)))
trailing_panel.fig.savefig("/stor/home/ka33933/figures/trailing_bar_test.png")

# %% Stacked bar

# %%

pos
# %%
aln = PSL.PslRowAln(pa_test[113040], ques, tars)
aln.show()
leading_query_r(aln)
aln = PSL.PslRowAln(pa_test[92109], ques, tars)
leading_query_r(aln)

# %%
print(aln.r)
aln.r.qSeq[0:117]
aln.show()
aln._block_coords()
aln.trailing_query_pos()
aln.show()

# %%
rs = np.array(rs)
nbins=abs(min(rs))+max(rs)
nbins
hist = np.histogram(rs, bins=nbins)
np.array(list(hist[1])+[521])

# %%
t = pmbip.Theme()
t.axislabels.fontsize=6
panel = pmbip.Paneler(1,1,figsize=(9,3), format = "png")
panel.next_ax().bar(x=hist[1], height=np.array(list(hist[0])+[0]))
t.apply_to(panel.current_ax)
panel.fig.savefig("/stor/home/ka33933/figures/overhang_hists.png")
# %%

# a=PslRowAln(pa_test[113421], ques, tars)
# print(a.r)
# a.show()

hist

>>>>>>> 2e0029999e68121b368a4b34060db000b226a28e

print(i)
# %%

aln=PslRowAln(pa.rows[137], ques, tars)
print(aln.r)
blks = aln._blocks()
aln.qSize_gapped
aln.tSize_gapped
aln.show()
<<<<<<< HEAD
aln.trailing_query_pos()
=======

np.histogram(np.array(list(aln.trailing_query_pos())+list(aln.trailing_query_pos())))

>>>>>>> 2e0029999e68121b368a4b34060db000b226a28e
sum([b.t_insert_bases for b in aln._blocks()])
sum([b.q_insert_bases for b in aln._blocks()])

aln.
aln=PslRowAln(pa.rows[1267], ques, tars)
aln=PslRowAln(pa.rows[137], ques, tars)
print(aln.r)
len(re.sub("[^ATCG]", "", aln._aln_repr()[0]))
len(re.sub("[^ATCG]", "", aln._aln_repr()[1]))
al
aln._aln_track()
aln.show(chunksize=150)
aln.r.tBaseInsert
aln.leading_query()
aln.trailing_query()
al
aln._blocks()[8].reaches_end_of_query()
aln._aln_repr_blocks()
mm=[b.mismatches(qSeq=ques[aln.r.qName], tSeq=tars[aln.r.tName]) for b in aln._blocks()]
mm=[b.mismatches(qSeq=aln.r.qSeq, tSeq=aln.r.tSeq) for b in aln._blocks()]
mm
len(mm)
# %%

print([b.q_has_insert for b in blks])
print([b.t_has_insert for b in blks])
print([b.q_insert_bases for b in blks])
print([b.t_insert_bases for b in blks])
# %%
# # Calulate number bases inserted manually
#     if block.qbS != block.psl_bs:
#         assert block.psl_bs < block.qbS
#         qbi += block.qbS - block.psl_bs
#     if block.tbS != block.psl_bs:
#         assert block.psl_bs < block.tbS
#         tbi += block.tbS - block.psl_bs


print(aln.r)
qbi
tbi
# %%

print(aln.r)
aln._query_insert_point()
aln.show()
print(pa.rows[0])
# %%
n_mis = 0
for qs,ts,bs in blks:

    print(qs,ts,bs)


# %%
blks
block_cov = []
block_cov = [(qs,qs+bs) for qs,ts,bs in blks]
inbetween = []
last_end = 0
for s,e in block_cov:
    inbetween.append((last_end,s))
    last_end = e
inbetween.append((last_end,aln.r.qSize))

block_cov
inbetween
# %%
for i in range(0,10):
    print(aln.r.qSeq[inbetween[i][0]:inbetween[i][1]])
    print(aln.r.qSeq[block_cov[i][0]:block_cov[i][1]])

# %%
inbetween
block_cov
# %%
for block in blks:
    qs,ts,bs = block
    for i in list(range(qs,qs+bs)):
        block_cov.append(i)

[x for x in range(0,len(qstr)) if x not in block_cov]
# %%
block_cov
pd.Series(aln.r.tStarts).diff(-1)
d = pd.concat([pd.Series(aln.r.tStarts).diff(-1), pd.Series(aln.r.qStarts).diff(-1)], axis=1)
d.columns = ["tsd", "qsd"]
d["t_diff"] = d["tsd"]-d["qsd"]
pd.concat([pd.DataFrame(blks), d], axis=1)
d
tstr =list(aln.r.tSeq) 
qstr =list(aln.r.qSeq)
tEnd =  aln.r.tEnd
qEnd =  aln.r.qEnd
for bidx,block in enumerate(aln._block_coords()):
    print(block)
    qs,ts,bs = block
    tb=aln.r.tSeq[ts:ts+bs]
    qb=aln.r.qSeq[qs:qs+bs]
    for i in range(0,bs):
        if tb[i] != qb[i]:
            n_mis+=1
    t_diff= d.iloc[bidx,:]["t_diff"]
    if not pd.isna(t_diff):
        if t_diff > 0:
            inc = int(abs(t_diff))
            tstr = tstr[0:ts+bs] + ["-"]*inc + tstr[ts+bs:]
            tEnd+=inc
        elif t_diff < 0:
            inc = int(abs(t_diff))
            qstr = qstr[0:qs+bs] + ["-"]*inc + qstr[qs+bs:]
            qEnd+=inc
        else:
            pass
    # print(tb,qb,t_diff)
print(f"{''.join(tstr[aln.r.tStart:tEnd])}\n{''.join(qstr[aln.r.qStart:qEnd])}")

blks
n_mis
# %%
qa = list(aln.r.qSeq[aln.r.qStart:aln.r.qEnd])
ta = list(aln.r.tSeq[aln.r.tStart:aln.r.tEnd])
qa
ta

new = []
for i in range(0,len(qa)):
    if qa[i] != ta[i]:
        qa.insert(i, "-")
    else:
        pass

print(qa)

# %%
aln.r.tSeq[280:378]
pmbi.bio.dna.revcomp(aln.r.qSeq)
ques[pa.rows[0].qName]
aln.df.loc[range(0,400),:]
# %%
a = pd.DataFrame({
    "aln": [" "]*1200,
    "tseq": [" "]*1200, 
    "qseq": [" "]*1200
    })
a.index = range(-400,800)


insert(a, range(0,r.tSize), "tseq", list(tars[r.tName]))
insert(a, t_range_q, "qseq", list(ques[r.qName]))



a.loc[t_range_q,"qseq"] = {k:v for k,v in zip(t_range_q, list(ques[r.qName]))}
a.loc[t_range_q,:]



# %%

# %%
def target_aln_repr(p):
    s = list("."*(p.tSize))
    for i,ts in enumerate(p.tStarts):
        for x in range(ts, (ts)+p.blockSizes[i]):
            s[x] = '|'
    chunk_size = 100
    breaks = list(range(0, p.tSize, chunk_size))
    return [''.join(s[i:i+chunk_size:]) for i in breaks]


print(pa.rows[1])
187-243
151-95
pa.rows[0].tStarts
pa.rows[0].tStart
pa.rows[0].qStart
print(dir(pa.rows[0]))
def aln_repr(p):
    s = list("."*(p.tSize))
    q = list(" "*(p.tSize))
    head_hang = 0
    if p.tStart-p.qStart<0:
        head_hang = abs(p.tStart-p.qStart)
    t_pos = range(0-head_hang, p.tSize)
    q_pos = range(p.qStart-p.tStart, p.tSize)
    # print(t_pos)
    # print(len(t_pos), len(q_pos))
    # print(t_pos, s, q_pos, q)
    print(pd.Series(t_pos), pd.Series(q_pos))
    df = pd.DataFrame({
        "t_pos": pd.Series(t_pos),
        "t_repr": s, 
        "q_pos": pd.Series(q_pos),
        "q_repr": q
        })
    # print(df)
    df.loc[np.where(df.q_pos>p.qSize)[0], "q_pos"] = np.inf
    print(df)
    # if head_hang < 0:
    #     head = pd.DataFrame({"t_repr": " "*abs(head_hang), "q_repr": "Q"*abs(head_hang)})
    #     df = pd.concat([head, df], axis=1)
    return(df)
    s_head = []
    s_tail = []
    q_head = []
    q_tail = []
    chunk_size = 100
    print(p.qStart, p.qEnd)
    if p.qStart < p.qEnd:
        for i,ts in enumerate(p.tStarts):
            print(ts-p.qStarts[i], (ts)+p.blockSizes[i])
            for x in range(ts-p.qStarts[i], (ts)+p.blockSizes[i]):
                if x<0:
                    s_head = [" "]+s_head
                    q_head = ["Q"]+q_head
                else:
                    if x <= (p.tSize-1):
                        if q[x] != ' ':
                            add = 'q'
                        else:
                            add = 'Q'
                        if x>=ts:
                            s[x] = '|'
                        q[x] = add
                    else:
                        s_tail = s_tail+[" "]
                        q_tail = q_tail+["Q"]
        s = s_head + s + s_tail
        q = q_head + q + q_tail
        breaks = list(range(0, len(s), chunk_size))
        return ( [''.join(s[i:i+chunk_size:]) for i in breaks], [''.join(q[i:i+chunk_size:]) for i in breaks] )
    else:
        return (['N','N'], ['N','N'])

# %%
# p.tSize-p.tStarts[0]
n_gt1_block = 0
# for p in pa.rows[9991:9992]:
# for p in pa.rows[9980:10000]:
for p in pa.rows:
    if p.blockCount>1:
        n_gt1_block+=1
    if isinstance(p.qBaseInsert, list):
        print(p)
        break
    # print('X'*(p.qStart-1) + '|'*(p.qSize-p.qStart))
    print(
          # p.aln_len()["query"], 
          # p.aln_len()["target"],
          p.qStart,
          p.qEnd,
          p.tStart,
          p.tEnd,
          p.blockSizes,
          # p.blockCount,
          [(t,q) for t,q in zip(p.qStarts,p.tStarts)],
          # p.tStarts,
          # p.qNumInsert,
          # p.qBaseInsert,
          # p.tNumInsert, 
          # p.tBaseInsert,
          # p.strand,
          p.dist_from_target_ends(),
          p.qSize,
          p.tSize
          )
    # t,q = aln_repr(p)
    df = aln_repr(p)
    # print(df)
    # for i in range(0,len(t)):
    #     print(t[i])
    #     print(q[i])
    #     print('')
    # print("\n\\\\\\\\\\\\\\\\\\\\\\\ \n")
    # break

# %%
n_gt1_block
len(pa.rows)

print(len(pa.rows))
pa.rows[0]
df_r1 = pd.DataFrame(ldict) 

# %%
class PslRow(object):
    def __init__(self,
                 row,
                ):
        self.target_id=row.t_id
        self.target_len=row.t_len
        self.target_start=row.t_start
        self.target_end=row.t_end
        self.query_id=row.q_id
        self.query_start=row.q_start
        self.query_end=row.q_end
        self.mismatch=row.mismatch
        self.length=row.length
    def __repr__(self):
        return f"{self.query_id}--{self.target_id}"
    def aln_len(self):
        return {
                "query": abs(self.query_end - self.query_start),
                "target": abs(self.target_end - self.target_start),
                }
    def ngaps(self):
        aln_len = self.aln_len()
        return {
                "query": self.length - aln_len["query"],
                "target": self.length - aln_len["target"],
                }
    def query_orient(self):
        if self.query_start > self.query_end:
            return "-"
        elif self.query_start < self.query_end:
            return "+"
        else:
            return "0"

# %%
df_r1.shape
for i,r in enumerate(df_r1.itertuples()):
    p = PslRow(r)
    al = p.aln_len()
    if al["query"] != al["target"]:
        break
    if i == df_r1.shape[0]-1:
        print("END")


# %%
df_r1
df_r1.iloc[684996,:]
p.ngaps()
p.target_id
df_r1
# %%
df_r1.apply(lambda r: PslRow(r), axis=0)
df_r1
df_r1["q_aligned_len"] = abs(df_r1["q_end"] - df_r1["q_start"])
df_r1["t_aligned_len"] = abs(df_r1["t_end"] - df_r1["t_start"])
df_r1["q_gaps"] = df_r1["length"] - df_r1["q_aligned_len"] 
df_r1["t_gaps"] = df_r1["length"] - df_r1["t_aligned_len"] 
df_r1[df_r1["q_gaps"]==df_r1["t_gaps"]]
df_r1[df_r1["q_aligned_len"]!=df_r1["t_aligned_len"]]

# %%
aln_rank = df.groupby("t_id").size().sort_values(ascending=False)
aln_rank
n_alignments_by_target = np.log10(aln_rank.to_numpy())
most_mapped = df[df.t_id == aln_rank.index[0]]

aln_rank_q = df.groupby("q_id").size().sort_values(ascending=False)
aln_rank_q.to_numpy()
n_alignments_by_query = np.log10(aln_rank_q.to_numpy())
most_mapped_q = df[df.t_id == aln_rank_q.index[0]]

# %%
def label(ax, xlab, ylab):
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ax

# %%
# def overhang_reads(df, cutoff, strand = "+"):
#     if strand == "+":
#         correct_orientation = df[df.q_start<df.q_end]
#         within = correct_orientation[correct_orientation.t_start>=(correct_orientation.t_len-cutoff)]
#         return within
#     elif strand == "-":
#         correct_orientation = df[df.q_start>df.q_end]
#         within = correct_orientation[correct_orientation.t_end<=cutoff]
#         return within

# %%
def overhang_reads_unstranded(df, cutoff):
    left = df[(df.t_end<=cutoff) & (df.q_start>df.q_end)]
    right = df[(df.t_start>=(df.t_len-cutoff)) & (df.q_start<df.q_end)]
    return right,left

# %%
right,left = overhang_reads_unstranded(df[df.t_id==examp], 150)
left["q_max_extent"] = left.t_end - (left.q_len-(left.q_len-left.q_start)) - left.q_gaps
right["q_max_extent"] = right.t_start + (right.q_len-right.q_start) + right.q_gaps
right["q_max_extent"].min()

left["q_max_extent"].max()
left[left.q_max_extent>100] 

# %%
colors = palettable.cartocolors.qualitative.Antique_10.mpl_colors
colors
# %%
t = pmbip.Theme()
t.inner
t.axislabels.fontsize=6
panel = pmbip.Paneler(2,3,figsize=(9,3), format = "png")
panel.next_ax().hist(n_alignments_by_target, 100)
panel.current_ax = label(panel.current_ax, "log10 n_alignments", "count")
t.apply_to(panel.current_ax)

t=pmbip.Theme()
t.title.text = aln_rank.index[0]
t.title.fontdict = {"fontsize":8}
t.axislabels.fontsize=6

panel.next_ax().hist(most_mapped["t_start"], 50, color=colors[0])
panel.current_ax = label(panel.current_ax, "position", "count")
t.apply_to(panel.current_ax)

panel.next_ax().hist(most_mapped["t_end"], 50, color=colors[1])
panel.current_ax = label(panel.current_ax, "position", "count")
t.apply_to(panel.current_ax)

t = pmbip.Theme()
t.inner
t.axislabels.fontsize=6
panel.next_ax().hist(n_alignments_by_query, 100)
panel.current_ax = label(panel.current_ax, "log10 n_alignments", "count")
t.apply_to(panel.current_ax)


t = pmbip.Theme()

panel.next_ax().hist(left["q_max_extent"], 50, color=colors[3])
panel.current_ax.hist(right["q_max_extent"], 50, color=colors[4])
# panel.current_ax.vlines(x=[0,400], ymin=0, ymax=1500)
panel.current_ax = label(panel.current_ax, "position", "count")
t.apply_to(panel.current_ax)

panel.fig.savefig("/stor/home/ka33933/figures/blat_hists.png")
# %%

