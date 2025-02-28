# %%
import shelve
import copy
import pmbi.plotting as pmbip
import pandas as pd
import argparse
import numpy as np
import palettable
import re

import vcf

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
vcf_reader = vcf.Reader(open("/home/amsesk/super1/cdiff_evo/combined/combined.vcf", 'r'))
# next(vcf_reader.reader)
# for r in vcf_reader:
#     r =r 
#     break

# %%
from Bio import SeqIO
ref = SeqIO.parse("/home/amsesk/super1/cdiff_evo/ref/GCF_019880625.1_ASM1988062v1_genomic.fna", "fasta")
for r in ref:
    print(dir(r))
    if r.name == "NZ_JACGTL010000001.1":
        core = r
        chrom_len = len(r.seq)

# %%
include = "bbb"
skip_na_mqrs = False
select_contigs = ["NZ_JACGTL010000001.1"]
ldict = []
n_pos = 0
for record in vcf_reader:
    if include is not None:
        if record.CHROM not in select_contigs:
            continue
    if "END" in record.INFO:
        block_range = range(record.POS-1, record.INFO["END"])
    else:
        block_range = range(record.start, record.end)
    n_pos += len(block_range)
    if n_pos > max(block_range)+1:
        n_pos -= (n_pos-(max(block_range)+1))
    if len(record.get_hets()) == 1:
        print(record.samples)
        for sample in [s for s in record.samples if re.search("(^Sample|^Ancestor)", s.sample) is not None]:
            ad = sample["AD"]
            dp = sample["DP"]
            if "MQRankSum" not in record.INFO:
                if skip_na_mqrs:
                    continue
                else:
                    mqrs = np.nan
            else:
                mqrs = record.INFO["MQRankSum"]
            print(sample)
            if ad is None:
                continue
            s = sum(ad)
            if s == 0:
                p = np.nan
                q = np.nan
            else:
                p = ad[0]/s
                q = ad[1]/s
            ldict.append({
                "sample": sample.sample,
                "CHROM": record.CHROM, 
                "POS": record.POS,
                "AD": ad,
                "p": p,
                "q": q,
                "DP": dp,
                "MQRS": mqrs,
                })

# %%
df = pd.DataFrame(ldict)
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

