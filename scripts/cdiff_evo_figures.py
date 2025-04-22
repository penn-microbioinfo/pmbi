
# %%
pd.options.display.max_rows=100
pd.set_option("display.max_rows", 100)
df_filt_long[(df_filt_long["POS"]== 83429) & (df_filt_long["snpEff__Annotation_Impact"]=="LOW")][["sample","snpEff__Annotation_Impact"]]
df_filt_long.columns

# %% SNP density
# %% FUNC: Return lists of intergers constituting sliding windows {{{
def sliding_window(full_size, chunk_size, step_size=1, one_based=True):
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
        for start in np.arange(0, full_size, step_size):
            stop = start + chunk_size
            if stop == full_size:
                break
            else:
                # yield np.arange(start, stop, 1)
                yield (start, stop)

# %% 
def is_all_ref(afs):
    is_ref = list()
    for af in afs:
        if af[0]==1.0:
            is_ref.append(True)
        else:
            is_ref.append(False)
    return pd.Series(is_ref)
                
def pos_in_range(pos, r): 
    return pos>=min(r) and pos<max(r)

# %%
df_filt[["sample_num", "sample_day"]] = df_filt["sample"].str.split(".", expand=True)
colors = pd.read_csv("/Users/amsesk/penn-microbioinfo/cdiff_evo/colors.tsv", sep="\t")
colors["sample"]=colors["sample"].str.replace("SM004_", "Sample")
colors = colors.set_index("sample")
allsamp = df_filt["sample"].unique()


# %%
max_snp_densities = []
for s in allsamp:
    sn = s.split(".")[0]
    print(s)
    sample_has = df_filt[(df_filt["sample"] == s) & (~is_all_ref(df_filt["AF"]))].POS.reset_index(drop=True)
    window_size=10000
    dens = []
    for i,(start,stop) in enumerate(sliding_window(full_size=ref_len, chunk_size=window_size, step_size=10, one_based=False)):
        if i%10000==0:
            print(i)
        snp_in = sample_has.apply(lambda pos: pos>=start and pos<stop )
        dens.append(np.where(snp_in)[0].shape[0]/window_size)
    dens_arr = np.array(dens)
    dens = pd.DataFrame({"pos": range(0,len(dens_arr)), "snp_density": dens_arr})
    dens["pos"] = dens["pos"]*10
    max_snp_densities.append(dens["snp_density"].max())
    panel = pmbip.Paneler(1,1,(5,2), dpi=800)
    ax = panel.next_ax()
    ax.plot(dens["pos"], dens["snp_density"], linewidth=0.75, c=colors.loc[sn,"color"])
    # ax.set_ylim((0,0.0085))
    ax.set_ylabel("SNP Density")
    ax.set_xlabel("Position")
    ax.set_title(s)
    panel.fig.savefig(f"/Users/amsesk/penn-microbioinfo/cdiff_evo/snp_density_figs/diff_y/{s}_density_plts.pdf")



# sdf = df_filt[df_filt["sample"]=="Sample1.Day1"]
# is_all_ref(sdf["AF"])
# max(max_snp_densities)
##################################################################
##################################################################
# %% STOP  {{{
##################################################################
##################################################################

# %%
index_cols = ["CHROM", "POS", "REF", "allele", "n_alleles"] + [
    f"snpEff__{xx}" for xx in ann_format_columns
]
df_filt_long.pivot(index=index_cols, columns="sample", values="AF").reset_index()

# %% CHUNK: Assign gene names based on POS
gtf._features


df_filt["gene_name_or_id"] = df_filt["POS"].apply(
    pos_to_gene_name, args=(gtf._features.attributes(),)
)
df_filt.columns

# gtf_df[(gtf_df["start"] <= 11000) & (gtf_df["end"] >= 12000)]
# %% CHUNK: Rearrange columns and print
df_filt = df_filt[
    [
        "sample",
        "CHROM",
        "POS",
        "gene_name_or_id",
        "REF",
        "alleles",
        "AD",
        "AF",
        "DP",
        "MQRS",
        "ANN",
    ]
]
df_filt.to_csv(
    "/mnt/raid1/webServer/export/amsesk/cdiff_evo_long.csv", sep=",", index=False
)
df_filt.columns
df_filt.pivot(
    columns="sample",
    index=["CHROM", "POS", "gene_name_or_id", "REF", "MQRS"],
    values="AF",
).to_csv("/mnt/raid1/webServer/export/amsesk/cdiff_evo_wide.csv", sep=",", index=False)


#################################
#################################

# %% CHUNK:
df_filt[(df_filt.sample == "Sample11.Day120") & (df_filt.POS == 93516)]
df_filt
df_filt[(df_filt["sample"] == "Sample12.Day120") & (df_filt["POS"] == 93516)]
df_filt[(df_filt["POS"] == 93516)]
# %% Convert to long format
df_filt_long = df_filt.explode(["alleles", "AD", "AF"]).rename(
    columns={"alleles": "allele"}
)
df_filt_long = df_filt_long.reset_index(drop=True)

gtf_df.columns


# %% Add gene names from GTF
df_filt_long["snpEff__Gene_Name"].apply(
    lambda x: gtf_df[gtf_df["gene_id"] == x]["gene"][0] if x is not None else None
)
df_filt_long = df_filt_long[~pd.isnull(df_filt_long["snpEff__Allele"])]
gtf_df.drop_duplicates()

df_filt_long["snpEff__Gene_Name"].value_counts()

# %% Add gene names from GTF
gene_names = []
gene_rename = gtf_df[["gene_id", "gene"]]
gene_rename.columns = ["snpEff__Gene_ID", "snpEff__Gene_Name"]
df_filt_long = df_filt_long.drop(columns="snpEff__Gene_Name").merge(
    gene_rename, how="left", on="snpEff__Gene_ID"
)

df_filt_long.to_csv(
    "/home/amsesk/figures/cdiff_evo/snp_table_noAncestor.tsv", sep="\t", index=False
)


# %% Drop modifier-class variant annotations for subset table of putatively more interesting variants
df_filt_long_noMod = df_filt_long[
    df_filt_long["snpEff__Annotation_Impact"] != "MODIFIER"
][
    [
        "sample",
        "CHROM",
        "POS",
        "REF",
        "allele",
        "AF",
        "AD",
        "DP",
        "MQRS",
        "snpEff__Gene_Name",
    ]
].reset_index(
    drop=True
)
df_filt_long_noMod
assert (
    df_filt_long_noMod[["CHROM", "POS", "REF", "allele", "snpEff__Gene_Name"]]
    .drop_duplicates()
    .value_counts()
    == 1
).all()

df_filt_long_noMod.pivot_table(
    columns="sample", index=["POS", "REF", "allele", "snpEff__Gene_Name"], values="AF"
).reset_index()
# .to_csv("/home/amsesk/figures/cdiff_evo/snp_table_noAncestor_wide.tsv", index=False, sep=",")
# df_filt_long_noMod["sample"].unique().shape
# df_filt_long[df_filt_long["POS"]==200598]["snpEff__Gene_Name"]
# 2800/15


# %%
["CHROM", "POS", "REF", "allele"] + list(df_filt_long.columns[snpEff_column_idx])
(
    df_filt_long[df_filt_long["snpEff__Annotation_Impact"] != "MODIFIER"][
        ["CHROM", "POS", "REF", "allele"]
        + list(df_filt_long.columns[snpEff_column_idx])
    ].drop_duplicates()
)
snpEff_column_idx = np.where(
    [x.startswith("snpEff__") for x in df_filt_long.columns.values]
)[0]
df_filt_long[df_filt_long["snpEff__Annotation_Impact"] != "MODIFIER"][
    ["CHROM", "POS", "REF", "allele"] + list(df_filt_long.columns[snpEff_column_idx])
]
d
# %%

df_filt_long.apply(
    lambda row: row["ANN"][row["allele"]] if row["ANN"] is not None else None, axis=1
)
df_filt
df_filt_long
df_filt[df_filt["n_alleles"] > 2]
df_filt.loc[120]
df_filt["ANN"].iloc[0]
df_filt.iloc[50:60, :]


df_filt["POS"].unique()
df_filt[df_filt.POS == 3414401]
df_filt_posBySample = df_filt.pivot(index="POS", columns="sample", values=["AD", "DP"])
df_filt_posBySample


# %%
plt_n_alleles = df_filt[["POS", "n_alleles"]].drop_duplicates()
plt_dp = df_filt_posBySample["DP"].apply(lambda x: sum(x) / len(x), axis=1)
df_filt
panel = pmbip.Paneler(nrow=2, ncol=2, figsize=(4, 4))
panel.next_ax().hist(plt_n_alleles["n_alleles"])
panel.current_ax.set_xlabel("n alleles")
panel.current_ax.set_ylabel("frequency")
panel.next_ax().hist(plt_dp)
panel.current_ax.set_xlabel("mean depth at site")
panel.current_ax.set_ylabel("frequency")
pal = palettable.scientific.diverging.Roma_10.mpl_colors
ax = panel.next_ax()
for i, nall in enumerate(df_filt["n_alleles"].unique()):
    sub = df_filt[df_filt["n_alleles"] == nall]
    ax.scatter(
        np.log10(sub["DP"]),
        [x[0] for x in sub["AF"]],
        s=0.5,
        marker=MarkerStyle(marker="o", fillstyle="full"),
        c=pal[i],
    )
panel.current_ax.set_xlabel("depth")
panel.current_ax.set_ylabel("freq(REF)")
pal = palettable.scientific.diverging.Tofino_10.mpl_colors
ax = panel.next_ax()
for i, nall in enumerate(df_filt["n_alleles"].unique()):
    sub = df_filt[df_filt["n_alleles"] == nall]
    ax.scatter(
        np.log10(sub["DP"]),
        [x[1] for x in sub["AF"]],
        s=0.5,
        marker=MarkerStyle(marker="o", fillstyle="full"),
        c=pal[i],
    )
panel.current_ax.set_xlabel("depth")
panel.current_ax.set_ylabel("freq(ALT_major)")
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/filtered_snps_stats.pdf")
# %%
# df.POS.unique().shape

# %%
df
samps = df["sample"].unique()
ancest_pq = df[df["sample"] == samps[0]][["p", "q"]]
df[(df["sample"] == samps[0])]
df[(df["sample"] == samps[0]) & (df["q"] >= 0.05)]
pos_not_fixed_in_anc = df[~df["POS"].isin(pos_anc)]
pos_not_fixed_in_anc

ancest_pq = df[df["sample"] == samps[3]][["p", "q"]]
ancest_pq

len(samps)
# %%
panel = pmbip.Paneler(4, 4, (9, 9))
panel_d = pmbip.Paneler(4, 4, (9, 9))
for samp in samps:
    # tp = pos_not_fixed_in_anc[pos_not_fixed_in_anc["sample"]==samp]
    tp = df[df["sample"] == samp]
    [x for x in tp.iterrows() if np.isnan(x[1]["p"])]
    list(tp.iterrows())[0]
    kde = gaussian_kde(tp["p"])
    x = np.arange(0.0, 1.0, 0.01)
    y = kde(x)
    print(x, y)
    panel.next_ax().scatter(tp["p"], tp["DP"], s=0.1)
    panel_d.next_ax().scatter(x, y, s=0.1)
    panel.current_ax.set_title(f"{samp} - n={tp.shape[0]}")

panel.fig.savefig("/home/amsesk/figures/cdiff_evo/ancestor_p_depth.png")
panel_d.fig.savefig("/home/amsesk/figures/cdiff_evo/ancestor_p_density.png")

# %%
poly = (
    pos_not_fixed_in_anc[["sample", "POS", "p"]]
    .pivot(index="sample", columns="POS", values="p")
    .fillna(1.0)
)
poly
panel = pmbip.Paneler(1, 1, (12, 4), dpi=1400)
pmbip.heatmap(poly.iloc[:], panel.next_ax(), xlab="pos", ylab="sample")
panel.current_ax.set_xticklabels(labels=[])
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/snp_hm_ancestorMinus.png", dpi=1400)

# %%
print(df)
print(df.shape)
df_finite = df[(~np.isnan(df["MQRS"])) & (~np.isnan(df["DP"]))]
print(df_finite)

panel = pmbip.Paneler(1, 2, (6, 3))
panel.next_ax().hist(df_finite["MQRS"], bins=100)
panel.next_ax().hist2d(
    x=df_finite["MQRS"],
    y=df_finite["DP"],
    bins=(25, 25),
)
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/mqrs_hist_filt.png")


# %%
def sliding_window(full_size, chunk_size, one_based=True):
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
            stop = start + chunk_size
            if stop == full_size:
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
    snp_binary = lil_array((len(range(0, n_pos)), 1))
    snp_binary[:] = 0
    snp_binary[snp_pos] = 1
    n_pos = snp_binary.shape[0]
    dens = lil_array((n_pos - window_size, 1))
    idx = 0
    for chunk in sliding_window(
        full_size=n_pos, chunk_size=window_size, one_based=False
    ):
        n_snp = snp_binary[chunk].nnz
        if n_snp == 0:
            d = 0.0
        else:
            d = np.true_divide(n_snp, window_size)
        dens[idx] = d
        idx += 1
    return dens.toarray().reshape(-1)


# %%
from joblib import Parallel, delayed

all_samples = df_finite["sample"].unique()
ret = Parallel(n_jobs=15)(delayed(_compute)(s, df, chrom_len, 1000) for s in samples)

# %%
named_ret = {k: v for k, v in zip(all_samples, ret)}
d = shelve.open("/home/amsesk/super1/cdiff_evo/snp_density_arrays.shelve")
d["1kb_sliding_window"] = copy.deepcopy(named_ret)
d.close()

# %%
mehmeh = shelve.open("/home/amsesk/super1/cdiff_evo/snp_density_arrays.shelve")
mehmeh["1kb_sliding_window"]


# %% CHUNK: Plot 2 {{{
x = np.arange(0, named_ret["Ancestor.Day0"].shape[0]) / 1e6
samples = ("Ancestor.Day0", "Sample3.Day1", "Sample3.Day120")
y_arrs = [named_ret[s] for s in samples]
pal = palettable.cmocean.sequential.Dense_20.mpl_colors[4:17:6]
y_max = np.concatenate(y_arrs).max()
y_max_ax = y_max * 1.10
panel = pmbip.Paneler(3, 1, (4, 3))
chrom_len
i = 0
for s, color in zip(samples, pal):
    y = named_ret[s]
    ax = panel.next_ax()
    ax.plot(x, y, linewidth=1.0, c=color)
    ax.text(4.25, 0.055, s, fontsize=6, c=color, ha="right")
    ax.set_xlabel("", fontsize=0)
    ax.set_ylim(0, y_max_ax)
    ax.set_xlim(0, 4.15)
    ax.set_xticks(ticks=np.arange(0, 4.5, 0.125), labels="", minor=True)
    ax.set_xticks(ticks=np.arange(0, 4.5, 0.5), labels="", minor=False)
    # ax.set_title(s)
    i = i + 1
    if i == len(samples):
        ax.set_xlabel("chromosome position (Mbp)", fontsize=5)
        ax.set_xticks(
            ticks=np.arange(0, 4.5, 0.5), labels=np.arange(0, 4.5, 0.5), minor=False
        )
    if i == 2:
        ax.set_ylabel("SNP density (SNP/kb)", fontsize=5)
panel.fig.savefig("/home/amsesk/figures/cdiff_evo/snp_density_Ancestor_Sample3.png")
# }}}

# %% CHUNK: Plot 2 {{{
named_ret_subtract = copy.deepcopy(named_ret)
x = np.arange(0, named_ret["Ancestor.Day0"].shape[0]) / 1e6
y_arrs = [named_ret[s] for s in samples]
samples = ("Ancestor.Day0", "Sample3.Day1", "Sample3.Day120")
named_ret_subtract["Ancestor_mAncestor"] = (
    named_ret_subtract["Ancestor.Day0"] - named_ret_subtract["Ancestor.Day0"]
)
named_ret_subtract["Sample3.Day1_mAncestor"] = (
    named_ret_subtract["Sample3.Day1"] - named_ret_subtract["Ancestor.Day0"]
)
named_ret_subtract["Sample3.Day120_mDay1"] = (
    named_ret_subtract["Sample3.Day120"] - named_ret_subtract["Sample3.Day1"]
)
samples = (
    ("Ancestor.Day0", "Ancestor_mAncestor"),
    ("Sample3.Day1", "Sample3.Day1_mAncestor"),
    ("Sample3.Day120", "Sample3.Day120_mDay1"),
)
samples
pal = palettable.cartocolors.qualitative.Prism_10.mpl_colors[3:10:3]
y_max = np.concatenate(y_arrs).max()
y_max_ax = y_max * 1.10
panel = pmbip.Paneler(3, 1, (4, 3))
chrom_len
i = 0
for s, color in zip(samples, pal):
    y = named_ret_subtract[s[0]]
    y_subtracted = named_ret_subtract[s[1]]
    ax = panel.next_ax()
    ax.plot(x, y, linewidth=1.0, c=color, alpha=0.2)
    ax.plot(x, y_subtracted, linewidth=1.0, c=color)
    ax.text(4.25, 0.055, s[0], fontsize=6, c=color, ha="right")
    ax.set_xlabel("", fontsize=0)
    ax.set_ylim(0, y_max_ax)
    ax.set_xlim(0, 4.15)
    ax.set_xticks(ticks=np.arange(0, 4.5, 0.125), labels="", minor=True)
    ax.set_xticks(ticks=np.arange(0, 4.5, 0.5), labels="", minor=False)
    i = i + 1
    if i == len(samples):
        ax.set_xlabel("chromosome position (Mbp)", fontsize=5)
        ax.set_xticks(
            ticks=np.arange(0, 4.5, 0.5), labels=np.arange(0, 4.5, 0.5), minor=False
        )
    if i == 2:
        ax.set_ylabel("SNP density (SNP/kb)", fontsize=5)
panel.fig.savefig(
    "/home/amsesk/figures/cdiff_evo/snp_density_Ancestor_Sample3_subtracted.png"
)
# }}}

# %% CHUNK: Plot 3 {{{
# day120_samps = np.array([s for s in all_samples if re.search("Day120", s) is not None])[[0,3,6,8]]
day120_samps = np.array(["Sample3.Day120", "Sample15.Day120"])
samples = np.concatenate((np.array(["Ancestor.Day0"]), day120_samps))
arr = {f"{s}_mAncestor": named_ret[s] - named_ret["Ancestor.Day0"] for s in samples}
sub_names = list(arr.keys())
orig_names = [
    str(s) for s in list(np.concatenate((np.array(["Ancestor.Day0"]), day120_samps)))
]
arr.update({s: named_ret[s] for s in orig_names})
arr

pal = palettable.cartocolors.qualitative.Vivid_10.mpl_colors[0:10:2]
list(arr.values())
y_max = np.concatenate(list(arr.values())).max()
y_max_ax = y_max * 1.10

i = 0
titles = orig_names.copy()
titles[1] = f"{titles[1]} (Het)"
titles[2] = f"{titles[2]} (KO)"
panel = pmbip.Paneler(3, 1, (4, 3))
for orig, sub, color, tit in zip(orig_names, sub_names, pal, titles):
    y = arr[orig]
    y_subtracted = arr[sub]
    ax = panel.next_ax()
    ax.plot(x, y, linewidth=1.0, c=color, alpha=0.2)
    ax.plot(x, y_subtracted, linewidth=1.0, c=color)
    ax.text(4.4, 0.053, tit, fontsize=5, c=color, ha="right")
    ax.set_xlabel("", fontsize=0)
    ax.set_ylim(0, y_max_ax)
    ax.set_xlim(0, 4.15)
    ax.set_yticks(
        ticks=np.linspace(
            0,
            0.05,
            num=3,
        ),
        labels=np.linspace(0, 0.05, num=3),
        minor=False,
        fontsize=4.5,
    )
    ax.set_xticks(ticks=np.linspace(0, 4.5, num=37), labels="", minor=True)
    ax.set_xticks(ticks=np.linspace(0, 4.5, num=10), labels="", minor=False)
    i = i + 1
    if i == len(samples):
        ax.set_xlabel("chromosome position (Mbp)", fontsize=5.5)
        ax.set_xticks(
            ticks=np.linspace(0, 4.5, num=10),
            labels=np.linspace(0, 4.5, num=10),
            minor=False,
            fontsize=4.5,
        )
    if i == 2:
        ax.set_ylabel("SNP density (SNP/kb)", fontsize=5.5)
panel.fig.savefig(
    "/home/amsesk/figures/cdiff_evo/snp_density_Ancestor_Samples3_15_subtracted.pdf"
)
np.linspace(0, 4.5, num=10)

# }}}


# %%


# %% TOSS {{{

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
# }}}

# }}}
