# %%

import shelve
import copy
import pmbi.plotting as pmbip
import pandas as pd
import argparse
import numpy as np
import palettable
from matplotlib.markers import MarkerStyle
import re

import vcf

# %%
import importlib
importlib.reload(pmbip)


# %%
parser = argparse.ArgumentParser()
parser.add_argument("gvcf", action="store", help="Path to the GVCF file.")
parser.add_argument("--include", action="store", help="Comma-separated list of contigs to include.")
parser.add_argument(
    "--skip_na_mqrs",

    action="store_true",
    required=False,
    help="Skips contigs that don't have a value in the MQRankSum field.",
)
args = parser.parse_args()
include = args.include
skip_na_mqrs = args.skip_na_mqrs
if include is not None:
    select_contigs = args.include.split(',')

# %%


vcf_reader = vcf.Reader(open(args.gvcf))
vcf_reader = vcf.Reader(open("/home/amsesk/super1/cdiff_evo/combined/combined.g.vcf", 'r'))
# next(vcf_reader.reader)
# for r in vcf_reader:
#     r =r 
#     break

# %%
from Bio import SeqIO
ref = SeqIO.parse("/home/amsesk/super1/cdiff_evo/ref/GCF_021378415.1_ASM2137841v1_genomic.fna", "fasta")

seqnames = []
for r in ref:
    print(r.name)
    # if r.name == "NZ_JACGTL010000001.1":
    #     core = r
    #     chrom_len = len(r.seq)

# %%
vcf_reader = vcf.Reader(open("/home/amsesk/super1/cdiff_evo/combined/combined.g.vcf", 'r'))
include = "NZ_CP059592.1"
skip_na_mqrs = False
select_contigs = include.split(',')
ldict = []
n_pos = 0
allele_counts = {}
i=1
for record in vcf_reader:
    if len(record.ALT) > 1:
        n_alt = len(record.ALT)-1
        if len([s["AD"] for s in record.samples if hasattr(s.data, "AD")]) == 0:
            continue
        for sample in record.samples:
            # print(sample.data)
            if sample["AD"] is None:
                ad = None
                # ad = [0]*(1+len(record.ALT)-1)
                # ref_count = 0
                # alt_count = [0]*n_alt
            else:
                # ad = sample["AD"][0:len(sample["AD"])]
                ad = sample["AD"][0:len(sample["AD"])]
                # ref_count = sample["AD"][0]
                # alt_count = sample["AD"][1:(1+(n_alt+1))]
                if ad[0] != 0 and all([x==0 for x in ad[1::]]):
                    print(record)
                    print(sample)
                # [0,1,2,3,4]
                # l[1:5-]
            dp = sample["DP"]
            if "MQRankSum" not in record.INFO:
                if skip_na_mqrs:
                    continue
                else:
                    mqrs = np.nan
            else:
                mqrs = record.INFO["MQRankSum"]
            newrow = {
                "sample": sample.sample,
                "CHROM": record.CHROM, 
                "POS": record.POS,
                "REF": record.REF,
                "ALT": record.ALT,
                "AD": ad,
                # "REF_count": ref_count,
                # "ALT_count": alt_count,
                "DP": dp,
                "MQRS": mqrs,
                }
            # print(newrow)
            ldict.append(newrow)
        i+=1
    # if i>=10:
    #     break
#     if include is not None:
#         if record.CHROM not in select_contigs:
#             continue
#     if "END" in record.INFO:
#         block_range = range(record.POS-1, record.INFO["END"])
#     else:
#         block_range = range(record.start, record.end)
#     n_pos += len(block_range)
#     if n_pos > max(block_range)+1:
#         n_pos -= (n_pos-(max(block_range)+1))
#     alleles = record.alleles
#     if len(alleles) in allele_counts:
#         allele_counts[len(alleles)] += 1
#     else:
#         allele_counts[len(alleles)] = 1
#     if len(alleles) == 2:
#         for sample in [s for s in record.samples if re.search("(^Sample|^Ancestor)", s.sample) is not None]:
#             ad = sample["AD"]
#             dp = sample["DP"]
#             if "MQRankSum" not in record.INFO:
#                 if skip_na_mqrs:
#                     continue
#                 else:
#                     mqrs = np.nan
#             else:
#                 mqrs = record.INFO["MQRankSum"]
#             if ad is None:
#                 continue
#             s = sum(ad)
#             if s == 0:
#                 continue
#                 # p = np.nan
#                 # q = np.nan
#             else:
#                 p = ad[0]/s
#                 q = ad[1]/s
#             ldict.append({
#                 "sample": sample.sample,
#                 "CHROM": record.CHROM, 
#                 "POS": record.POS,
#                 "AD": ad,
#                 "p": p,
#                 "q": q,
#                 "DP": dp,
#                 "MQRS": mqrs,
#                 })
#

# allele_counts
# %%
def variant_in(row, sample):
    if row[sample] is not None:
        return True
    else:
        return False

def variant_only_in(row, samples):
    variant_not_none_in = []
    # print(row)
    # return None
    for sample in row.index:
        if row[sample] is not None:
            variant_not_none_in.append(sample)
    if all([x in samples for x in variant_not_none_in]):
        return True
    else:
        return False
    

# %% Convert to DataFrame and pivot to POS x sample_AD 
df = pd.DataFrame(ldict)
df_ad = df.pivot(index="POS", columns="sample", values="AD")
n_var_total = df_ad.shape[0]
df_ad.columns


# %% Drop variants that are only variant in the control samples
control_samples = ["DNAfreewater1.20241030", "Extractblankswab1.20241030", "Extractemptywell1.20241030" , "mockdna1.20241030"]
df_ad = df_ad[~df_ad.apply(variant_only_in, axis=1, args=(control_samples,))]
df_ad = df_ad.drop(control_samples, axis=1)
n_var_noControl = df_ad.shape[0]

# %% Drop variants that are variant in the Ancestor
ancestor = "Ancestor.Day0"
df_ad = df_ad[~df_ad.apply(variant_in, axis=1, args=(ancestor,))]
n_var_noControl_noAncestor = df_ad.shape[0]

# %% Print some stats
print(f"""
      Total variant: {n_var_total}
      sans variant control-only variants: {n_var_noControl}
      sans variant in Ancestor: {n_var_noControl_noAncestor}
      """)

# %% Pull filtered variants from main Dataframe and filter out control and ancestor samples
df_filt = df[df["POS"].isin(df_ad.index)]
df_filt = df_filt[~df_filt["sample"].isin(control_samples+[ancestor])]
df_filt[["REF", "ALT", "AD"]]
df_filt["alleles"] = df_filt.apply(lambda row: ([row.REF]+list(row.ALT)), axis=1)

# %% Get rid of NONREF values
df_filt["alleles"] = df_filt["allelw
                             es"].apply(lambda row: row[0:-1:])
df_filt["AD"] = [row[0:-1:] if row is not None else None for row in df_filt["AD"]]

# %%
df_filt = df_filt[["sample", "CHROM", "POS", "REF", "alleles", "AD", "DP", "MQRS"]]

# %% Just to make sure that all of the non-None allele lists are the same length
assert df_filt.pivot(index="POS", columns="sample", values="AD").apply(lambda row: [len(x) for x in row if x is not None], axis=1).apply(lambda row: all([x==row[0] for x in row])).all()

# %% Add count of alleles per variant position
n_alleles = pd.DataFrame(df_filt.pivot(index="POS", columns = "sample", values="alleles").apply(lambda row: [len(x) for x in row if x is not None][0], axis=1)).reset_index().rename(columns={0: "n_alleles"})
df_filt = pd.merge(left=df_filt, right=n_alleles, how="left")
df_filt

# %% Convert AD==None into lists with counts of REF allele = DP, followed by 0's for the other alleles
df_filt["AD"] = df_filt.apply(lambda row: row["AD"] if row["AD"] is not None else [row["DP"]]+[0]*(row["n_alleles"]-1), axis=1)

# %%
def ad_to_af(ads, round_to=3):
    total = sum(ads)
    return [round(np.true_divide(x, total),4) for x in ads]

na_af_at = np.where([sum(x)==0 for x in df_filt["AD"]])

# NOTE: np.nan here means that the depth at that variant position was also 0
df_filt["AF"] = df_filt["AD"].apply(ad_to_af)

# %%

df_filt["POS"].unique()
df_filt_posBySample = df_filt.pivot(index="POS", columns="sample", values=["AD","DP"])
df_filt_posBySample

# %%
plt_n_alleles = df_filt[["POS", "n_alleles"]].drop_duplicates()
plt_dp = df_filt_posBySample["DP"].apply(lambda x: sum(x)/len(x), axis=1)
df_filt
panel = pmbip.Paneler(nrow=2, ncol=2, figsize=(4,4))
panel.next_ax().hist(plt_n_alleles["n_alleles"])
panel.current_ax.set_xlabel("n alleles")
panel.current_ax.set_ylabel("frequency")
panel.next_ax().hist(plt_dp)
panel.current_ax.set_xlabel("mean depth at site")
panel.current_ax.set_ylabel("frequency")
pal = palettable.scientific.diverging.Roma_10.mpl_colors
ax = panel.next_ax()
for i,nall in enumerate(df_filt["n_alleles"].unique()):
    sub = df_filt[df_filt["n_alleles"] == nall]
    ax.scatter(np.log10(sub["DP"]), [x[0] for x in sub["AF"]], s=0.5, marker=MarkerStyle(marker="o", fillstyle="full"), c=pal[i])
panel.current_ax.set_xlabel("depth")
panel.current_ax.set_ylabel("freq(REF)")
pal = palettable.scientific.diverging.Tofino_10.mpl_colors
ax = panel.next_ax()
for i,nall in enumerate(df_filt["n_alleles"].unique()):
    sub = df_filt[df_filt["n_alleles"] == nall]
    ax.scatter(np.log10(sub["DP"]), [x[1] for x in sub["AF"]], s=0.5, marker=MarkerStyle(marker="o", fillstyle="full"), c=pal[i])
panel.current_ax.set_xlabel("depth")
panel.current_ax.set_ylabel("freq(ALT_major)")
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/filtered_snps_stats.pdf")
# %%
# df.POS.unique().shape

# %%
df
samps = df["sample"].unique()
ancest_pq = df[df["sample"]==samps[0]][["p","q"]]
df[(df["sample"]==samps[0])]
df[(df["sample"]==samps[0]) & (df["q"] >= 0.05)]
pos_not_fixed_in_anc = df[~df["POS"].isin(pos_anc)]
pos_not_fixed_in_anc

ancest_pq = df[df["sample"]==samps[3]][["p","q"]]
ancest_pq

len(samps)
# %%
panel = pmbip.Paneler(4,4,(9,9))
panel_d = pmbip.Paneler(4,4,(9,9))
for samp in samps:
    # tp = pos_not_fixed_in_anc[pos_not_fixed_in_anc["sample"]==samp]
    tp = df[df["sample"]==samp]
    [x for x in tp.iterrows() if np.isnan(x[1]["p"])]
    list(tp.iterrows())[0]
    kde = gaussian_kde(tp["p"])
    x = np.arange(0.0,1.0,0.01)
    y = kde(x)
    print(x,y)
    panel.next_ax().scatter(tp["p"], tp["DP"], s = 0.1)
    panel_d.next_ax().scatter(x,y,s=0.1)
    panel.current_ax.set_title(f"{samp} - n={tp.shape[0]}")

panel.fig.savefig("/home/amsesk/figures/cdiff_evo/ancestor_p_depth.png")
panel_d.fig.savefig("/home/amsesk/figures/cdiff_evo/ancestor_p_density.png")

# %%
poly = pos_not_fixed_in_anc[["sample", "POS", "p"]].pivot(index="sample", columns = "POS", values="p").fillna(1.0)
poly
panel = pmbip.Paneler(1,1,(12,4), dpi=1400)
pmbip.heatmap(poly.iloc[:], panel.next_ax(), xlab="pos", ylab="sample")
panel.current_ax.set_xticklabels(labels=[])
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/snp_hm_ancestorMinus.png", dpi=1400)

# %%
print(df)
print(df.shape)
df_finite = df[ (~np.isnan(df["MQRS"])) & (~np.isnan(df["DP"]))]
print(df_finite)

panel = pmbip.Paneler(1,2,(6,3))
panel.next_ax().hist(df_finite["MQRS"], bins=100)
panel.next_ax().hist2d(x=df_finite["MQRS"], y=df_finite["DP"], bins = (25,25), )
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/mqrs_hist_filt.png")

# %%
def sliding_window(full_size, chunk_size, one_based = True):
    if one_based:
        raise NotImplementedError
        # for start in np.arange(1, full_size+1, 1):
        #     stop = start+chunk_size+1
        #     if stop>(full_size+1):
        #         break
        #     else:
        #         # chunks.append(np.arange(start, stop, 1))
        #         yield np.arange(start, stop, 1)
    else:
        for start in np.arange(0, full_size, 1):
            stop = start+chunk_size
            if stop==full_size:
                break
            else:
                yield np.arange(start, stop, 1)

# %%

# %s
from scipy.sparse import lil_array
sp_arr = csr_array(snp_binary.values)
sp_arr[0:1000].nnz

# %%

# r = range(0, chrom_len)
# snp_pos = sub["POS"].values
# snp_binary = lil_array( (len(r), 1) )
# snp_binary.shape
# snp_binary[:] = 0
# snp_binary[snp_pos] = 1
# %%
def _compute(sample, df, n_pos, window_size):
    sub = df[df["sample"] == sample]
    snp_pos = sub["POS"].values
    snp_binary = lil_array( (len(range(0,n_pos)), 1) )
    snp_binary[:] = 0
    snp_binary[snp_pos] = 1
    n_pos = snp_binary.shape[0]
    dens = lil_array( (n_pos-window_size,1) )
    idx = 0
    for chunk in sliding_window(full_size=n_pos, chunk_size=window_size, one_based=False):
        n_snp = snp_binary[chunk].nnz
        if n_snp == 0:
            d = 0.0
        else:
            d = np.true_divide(n_snp, window_size)
        dens[idx] = d
        idx+=1
    return dens.toarray().reshape(-1)

# %%
from joblib import Parallel, delayed
all_samples = df_finite["sample"].unique()
ret = Parallel(n_jobs=15)(
    delayed(_compute)(s, df, chrom_len, 1000) for s in samples
)

# %%
named_ret = {k:v for k,v in zip(all_samples, ret)}
d=shelve.open("/home/amsesk/super1/cdiff_evo/snp_density_arrays.shelve")
d["1kb_sliding_window"] = copy.deepcopy(named_ret)
d.close()

# %%
mehmeh=shelve.open("/home/amsesk/super1/cdiff_evo/snp_density_arrays.shelve")
mehmeh["1kb_sliding_window"]



# %% CHUNK: Plot 2 {{{
x = np.arange(0,named_ret["Ancestor.Day0"].shape[0])/1e6
samples = ("Ancestor.Day0", "Sample3.Day1", "Sample3.Day120")
y_arrs = [named_ret[s] for s in samples]
pal = palettable.cmocean.sequential.Dense_20.mpl_colors[4:17:6]
y_max = np.concatenate(y_arrs).max()
y_max_ax = y_max*1.10
panel = pmbip.Paneler(3,1,(4,3))
chrom_len
i = 0
for s,color in zip(samples, pal):
    y = named_ret[s]
    ax = panel.next_ax()
    ax.plot(x, y, linewidth = 1.0, c=color)
    ax.text(4.25,0.055,s, fontsize=6, c = color, ha="right")
    ax.set_xlabel("", fontsize=0)
    ax.set_ylim(0,y_max_ax)
    ax.set_xlim(0,4.15)
    ax.set_xticks(ticks=np.arange(0,4.5, 0.125), labels='', minor = True)
    ax.set_xticks(ticks=np.arange(0, 4.5, 0.5), labels='', minor = False)
    # ax.set_title(s)
    i = i + 1
    if i == len(samples):
        ax.set_xlabel("chromosome position (Mbp)", fontsize=5)
        ax.set_xticks(ticks=np.arange(0, 4.5, 0.5), labels=np.arange(0, 4.5, 0.5), minor = False)
    if i == 2:
        ax.set_ylabel("SNP density (SNP/kb)", fontsize=5)
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/snp_density_Ancestor_Sample3.png")
# }}}

# %% CHUNK: Plot 2 {{{
named_ret_subtract = copy.deepcopy(named_ret)
x = np.arange(0,named_ret["Ancestor.Day0"].shape[0])/1e6
y_arrs = [named_ret[s] for s in samples]
samples = ("Ancestor.Day0", "Sample3.Day1", "Sample3.Day120")
named_ret_subtract["Ancestor_mAncestor"] = named_ret_subtract["Ancestor.Day0"]-named_ret_subtract["Ancestor.Day0"]
named_ret_subtract["Sample3.Day1_mAncestor"] = named_ret_subtract["Sample3.Day1"]-named_ret_subtract["Ancestor.Day0"]
named_ret_subtract["Sample3.Day120_mDay1"] = named_ret_subtract["Sample3.Day120"]-named_ret_subtract["Sample3.Day1"]
samples = (("Ancestor.Day0", "Ancestor_mAncestor"), ("Sample3.Day1", "Sample3.Day1_mAncestor"), ("Sample3.Day120", "Sample3.Day120_mDay1"))
samples
pal = palettable.cartocolors.qualitative.Prism_10.mpl_colors[3:10:3]
y_max = np.concatenate(y_arrs).max()
y_max_ax = y_max*1.10
panel = pmbip.Paneler(3,1,(4,3))
chrom_len
i = 0
for s,color in zip(samples, pal):
    y = named_ret_subtract[s[0]]
    y_subtracted = named_ret_subtract[s[1]]
    ax = panel.next_ax()
    ax.plot(x, y, linewidth = 1.0, c=color, alpha=0.2)
    ax.plot(x, y_subtracted, linewidth = 1.0, c=color)
    ax.text(4.25,0.055,s[0], fontsize=6, c = color, ha="right")
    ax.set_xlabel("", fontsize=0)
    ax.set_ylim(0,y_max_ax)
    ax.set_xlim(0,4.15)
    ax.set_xticks(ticks=np.arange(0,4.5, 0.125), labels='', minor = True)
    ax.set_xticks(ticks=np.arange(0, 4.5, 0.5), labels='', minor = False)
    i = i + 1
    if i == len(samples):
        ax.set_xlabel("chromosome position (Mbp)", fontsize=5)
        ax.set_xticks(ticks=np.arange(0, 4.5, 0.5), labels=np.arange(0, 4.5, 0.5), minor = False)
    if i == 2:
        ax.set_ylabel("SNP density (SNP/kb)", fontsize=5)
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/snp_density_Ancestor_Sample3_subtracted.png")
# }}}

# %% CHUNK: Plot 3 {{{
# day120_samps = np.array([s for s in all_samples if re.search("Day120", s) is not None])[[0,3,6,8]]
day120_samps = np.array(["Sample3.Day120", "Sample15.Day120"])
samples = np.concatenate( (np.array(["Ancestor.Day0"]), day120_samps) )
arr = {f"{s}_mAncestor": named_ret[s]-named_ret["Ancestor.Day0"] for s in samples}
sub_names = list(arr.keys())
orig_names = [str(s) for s in list(np.concatenate((np.array(["Ancestor.Day0"]), day120_samps)))]
arr.update({s:named_ret[s] for s in orig_names})
arr

pal = palettable.cartocolors.qualitative.Vivid_10.mpl_colors[0:10:2]
list(arr.values())
y_max = np.concatenate(list(arr.values())).max()
y_max_ax = y_max*1.10

i = 0
titles = orig_names.copy()
titles[1] = f"{titles[1]} (Het)"
titles[2] = f"{titles[2]} (KO)"
panel = pmbip.Paneler(3,1,(4,3))
for orig,sub,color,tit in zip(orig_names, sub_names, pal, titles):
    y = arr[orig]
    y_subtracted = arr[sub]
    ax = panel.next_ax()
    ax.plot(x, y, linewidth = 1.0, c=color, alpha=0.2)
    ax.plot(x, y_subtracted, linewidth = 1.0, c=color)
    ax.text(4.4,0.053, tit, fontsize=5, c = color, ha="right")
    ax.set_xlabel("", fontsize=0)
    ax.set_ylim(0,y_max_ax)
    ax.set_xlim(0,4.15)
    ax.set_yticks(ticks=np.linspace(0,0.05,num=3,), labels=np.linspace(0,0.05,num=3), minor = False, fontsize= 4.5)
    ax.set_xticks(ticks=np.linspace(0,4.5, num=37), labels='', minor = True)
    ax.set_xticks(ticks=np.linspace(0, 4.5, num=10), labels='', minor = False)
    i = i + 1
    if i == len(samples):
        ax.set_xlabel("chromosome position (Mbp)", fontsize=5.5)
        ax.set_xticks(ticks=np.linspace(0,4.5,num=10), labels=np.linspace(0, 4.5, num=10), minor = False, fontsize=4.5)
    if i == 2:
        ax.set_ylabel("SNP density (SNP/kb)", fontsize=5.5)
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/snp_density_Ancestor_Samples3_15_subtracted.pdf")
np.linspace(0,4.5,num=10)

# }}}


# %%

