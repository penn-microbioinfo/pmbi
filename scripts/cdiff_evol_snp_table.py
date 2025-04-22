# %%
from Bio import SeqIO
import argparse
import copy
import gzip
import importlib
import re
import shelve

import numpy as np
import palettable
import pandas as pd
import vcf
from pathlib import Path
from matplotlib.markers import MarkerStyle

import pmbi.bio.gtf as pgtf
import pmbi.bio.vcf as pvcf
import pmbi.plotting as pmbip
from pmbi.util import is_ipython

pd.set_option("display.max_rows", 100)
importlib.reload(pmbip)
importlib.reload(pgtf)
importlib.reload(pvcf)

# %% FUNC: Class for checking and handling snpEff reference paths {{{
class snpEffRef(object):
    def __init__(self, refdir: Path):
        self.refdir = refdir
        self.sequences, self.genes, self.cds, self.proteins = self._parse_refdir(refdir)
    @staticmethod
    def _parse_refdir(refdir):
        expected_files = ["sequences.fa.gz",
                          "genes.gtf.gz",
                          "cds.fa.gz",
                          "protein.fa.gz"
                          ]
        fs = [f.name for f in refdir.iterdir()]
        if any([e not in fs for e in expected_files]):
            raise OSError("Missing files in snpEff reference directory")
        else:
            return tuple([refdir.joinpath(e) for e in expected_files])

# }}}

# %% FUNC: Get gene ID from GTF dataframe based on position {{{
def gene_id_from_pos(pos, gtf_df):
    after_start = gtf_df["start"] <= pos
    before_end = gtf_df["end"] >= pos
    entries = gtf_df[(after_start) & (before_end)]
    if entries.shape[0] == 0:
        return "intergenic"
    elif entries.shape[0] == 1:
        entry = entries.iloc[0, :]
        return entry["gene_id"]
    else:
        return ",".join(entries["gene_id"].tolist())

# }}}

# %% Parse arguments if not ipython or directly assign variables {{{
if is_ipython():
    skip_na_mqrs = False
    gvcf = Path("/home/amsesk/super1/cdiff_evo/combined/combined.snpeff.g.vcf")
    subtract_ancestor = False
    ref = snpEffRef(Path("/home/amsesk/pkgs/snpEff/data/cdiff_CD196/"))
    outdir = Path("/home/amsesk/figures/cdiff_evo")
else:
    parser = argparse.ArgumentParser()
    parser.add_argument("gvcf", action="store", help="Path to the GVCF file.")
    parser.add_argument(
        "-r", "--snpeff", action="store", help="Path to sneEff reference directory for reading reference fastas, gtfs, etc."
    )
    parser.add_argument(
        "--outdir", action="store", required=False, default="./", help="Path to write output files to."
    )
    parser.add_argument(
        "--subtract_ancestor", action="store", required=False, default="./", help="Whether or not to remove variants present in the ancestal sample. Defualt: False"
    )
    parser.add_argument(
        "--skip_na_mqrs",
        action="store_true",
        required=False,
        help="Skips contigs that don't have a value in the MQRankSum field. Default: False",
    )
    args = parser.parse_args()
    skip_na_mqrs = args.skip_na_mqrs
    gvcf = args.gvcf
    subtract_ancestor = args.subtract_ancestor
    outdir = Path(args.outdir)
    ref = snpEffRef(Path(args.snpeff))

# }}} 

# %% CHUNK: Define some sample pools
control_samples = [
    "DNAfreewater1.20241030",
    "Extractblankswab1.20241030",
    "Extractemptywell1.20241030",
    "mockdna1.20241030",
]
ancestor_sample = "Ancestor.Day0"

if subtract_ancestor:
    filter_out_samples = control_samples + [ancestor_sample]
    snp_table_out_base = "cdiff_evo_snp_table_long_subtractAnc.tsv"
else:
    filter_out_samples = control_samples
    snp_table_out_base = "cdiff_evo_snp_table_long.tsv"

# %% CHUNK: Read nucleotide reference
with gzip.open(ref.sequences, "rt") as fa_handle:
    nucl_ref = {s.name: str(s.seq) for s in SeqIO.parse(fa_handle, "fasta")}

# %% CHUNK: Parse the reference protein fasta
with gzip.open(ref.proteins, "rt") as fa_handle:
    prot_ref = {s.name: str(s.seq) for s in SeqIO.parse(fa_handle, "fasta")}

# %% CHUNK: Parse the referenece GTF, add protein sequences based on protein ids, and print attributes to TSV
with gzip.open(
    ref.genes, "rt"
) as gtf_handle:
    gtf = pgtf.FeatureFile.from_gtf(gtf_handle)

gtf_attr = gtf.attributes()
gtf_attr = gtf_attr.replace(np.nan, pd.NA)
gtf_attr["protein_sequence"] = gtf_attr["protein_id"].apply(
    lambda pid: prot_ref[pid] if not pd.isnull(pid) else pd.NA
)
gtf_attr = gtf_attr.drop_duplicates()
gtf_attr.to_csv(
    outdir.joinpath("cdiff_evo_gene_attributes.tsv"),
    sep="\t",
    index=False,
)

# %% CHUNK: Do initial reading and filtering of GCVF {{{
vcf_reader = vcf.Reader(open(gvcf, "r"))
ldict = []
n_pos = 0
allele_counts = {}
i = 1
for record in vcf_reader:
    if len(record.ALT) > 1:
        n_alt = len(record.ALT) - 1
        if len([s["AD"] for s in record.samples if hasattr(s.data, "AD")]) == 0:
            continue
        for sample in record.samples:
            if sample["AD"] is None:
                ad = None
            else:
                ad = sample["AD"][0 : len(sample["AD"])]
                if ad[0] != 0 and all([x == 0 for x in ad[1::]]):
                    print(record)
                    print(sample)
            dp = sample["DP"]
            if "MQRankSum" not in record.INFO:
                if skip_na_mqrs:
                    continue
                else:
                    mqrs = np.nan
            else:
                mqrs = record.INFO["MQRankSum"]
            if "ANN" in record.INFO:
                annot = record.INFO["ANN"]
            else:
                annot = None
            newrow = {
                "sample": sample.sample,
                "CHROM": record.CHROM,
                "POS": record.POS,
                "REF": record.REF,
                "ALT": record.ALT,
                "AD": ad,
                "DP": dp,
                "MQRS": mqrs,
                "ANN": annot,
            }
            ldict.append(newrow)
        i += 1

# }}}

# %% CHUNK: Do the rest

# %% CHUNK: %% Convert to DataFrame and pivot to POS x sample_AD
df = pd.DataFrame(ldict)
df["ANN"] = df["ANN"].apply(pvcf.split_snpeff_annots)
df_ad = df.pivot(index="POS", columns="sample", values="AD")
n_var_total = df_ad.shape[0]

# %% CHUNK: Drop variants that are only variant in the control samples
df_ad = df_ad[~df_ad.apply(pvcf.variant_only_in, axis=1, args=(control_samples,))]
df_ad = df_ad.drop(control_samples, axis=1)
n_var_noControl = df_ad.shape[0]

# %% CHUNK: Drop variants that are variant in the Ancestor, if subtract_ancestor is True
df_ad_subtract_Ancestral_variants = df_ad[~df_ad.apply(pvcf.variant_in, axis=1, args=(ancestor_sample,))]
n_var_noControl_noAncestor = df_ad_subtract_Ancestral_variants.shape[0]
if subtract_ancestor:
    df_ad = df_ad_subtract_Ancestral_variants

# %% Print some stats
print(
    f"""
      Total variant: {n_var_total}
      sans variant control-only variants: {n_var_noControl}
      sans variant in Ancestor: {n_var_noControl_noAncestor}

      Variants remaining now: {df_ad.shape[0]}
      """
)

# %% CHUNK: Pull filtered variants from main Dataframe and filter out control and ancestor samples
df_filt = df[df["POS"].isin(df_ad.index)]
df_filt = df_filt[~df_filt["sample"].isin(filter_out_samples)]
df_filt = df_filt.reset_index(drop=True)

# %% CHUNK: Remove NONREF 'allele' from ALT and AD columns, also convert ALT to list of str
df_filt["ALT"] = df_filt["ALT"].apply(lambda alts: [str(x) for x in alts[0:-1:]])
df_filt["AD"] = [row[0:-1:] if row is not None else None for row in df_filt["AD"]]

assert (
    df_filt["ALT"].apply(lambda x: [isinstance(xx, str) for xx in x]).all()
), "Not all lists in ALT are of type str"
assert (
    df_filt["ALT"].apply(lambda x: [len(xx) > 0 for xx in x]).all()
), "Not all lists in ALT are length > 0"

# %% CHUNK: Combine REF and ALT into a new 'alleles' column that contains a list of all alleles for a position
df_filt["alleles"] = df_filt.apply(lambda row: ([row.REF] + list(row.ALT)), axis=1)

assert (
    df_filt["alleles"].apply(lambda x: [isinstance(xx, str) for xx in x]).all()
), "Not all lists in ALT are of type str"
assert (
    df_filt["alleles"].apply(lambda x: [len(xx) > 0 for xx in x]).all()
), "Not all lists in ALT are length > 0"

# %% CHUNK: Make sure that all of the non-None allele lists are the same length for each position.
assert (
    df_filt.pivot(index="POS", columns="sample", values="AD")
    .apply(lambda row: [len(x) for x in row if x is not None], axis=1)
    .apply(lambda row: all([x == row[0] for x in row]))
    .all()
)

# %% CHUNK: Add a column that contains the number of alleles for each variant position
n_alleles = (
    pd.DataFrame(
        df_filt.pivot(index="POS", columns="sample", values="alleles").apply(
            lambda row: [len(x) for x in row if x is not None][0], axis=1
        )
    )
    .reset_index()
    .rename(columns={0: "n_alleles"})
)
df_filt = pd.merge(left=df_filt, right=n_alleles, how="left")

# %% CHUNK: %% Convert AD==None into lists with counts of REF allele = DP, followed by 0's for the other alleles
df_filt["AD"] = df_filt.apply(
    lambda row: (
        row["AD"]
        if row["AD"] is not None
        else [row["DP"]] + [0.0] * (row["n_alleles"] - 1)
    ),
    axis=1,
)

# %%CHUNK: Compute allele frequency AF from AD
# NOTE: np.nan here means that the depth at that variant position was also 0
df_filt["AF"] = df_filt["AD"].apply(pvcf.ad_to_af)
na_af_at = np.where([sum(x) == 0 for x in df_filt["AD"]])

# %% CHUNK: Convert to long format
df_filt_long = (
    df_filt.drop(columns=["ALT"])
    .explode(["alleles", "AD", "AF"])
    .rename(columns={"alleles": "allele"})
)

# %% CHUNK: Assign snfEff annotations to long format based on allele column
ann_format_columns = [
    re.sub(" ", "", x.strip())
    for x in re.search("['](.+)[']", vcf_reader.infos["ANN"].desc).group(1).split("|")
]
df_filt_long["ANN"] = df_filt_long.apply(
    lambda row: (
        [[pd.NA for x in range(0, len(ann_format_columns))]]
        if (row["allele"] == row["REF"] or row["allele"] in ["*"])
        else row["ANN"][row["allele"]]
    ),
    axis=1,
)
df_filt_long["ANN_colnames"] = df_filt_long["ANN"].apply(
    # lambda x: [f"snpEff__{xx}" for xx in ann_format_columns] if x is not None else None
    lambda x: [f"snpEff__{xx}" for xx in ann_format_columns]
)
df_filt_long["ANN_idx"] = df_filt_long["ANN"].apply(
    lambda x: list(range(0, len(x))) if len(x) > 0 else list(range(0, 1))
)

df_filt_long = df_filt_long.explode(["ANN", "ANN_idx"])
df_filt_long = df_filt_long.explode(["ANN", "ANN_colnames"])

# %%
df_filt_long = df_filt_long.pivot(
    columns="ANN_colnames",
    index=[
        "sample",
        "CHROM",
        "POS",
        "REF",
        "allele",
        "AD",
        "DP",
        "MQRS",
        "n_alleles",
        "AF",
        "ANN_idx",
    ],
    values="ANN",
).reset_index()

# %%
gene_gtf_df = gtf.update(gtf[gtf["feature"] == "gene"].reset_index(drop=True))
gene_gtf_df._features["gene_id"] = gene_gtf_df.get_attribute("gene_id")

df_filt_long["resident_gene_id"] = df_filt_long["POS"].apply(
    gene_id_from_pos, args=(gene_gtf_df._features,)
)

df_filt_long = df_filt_long[
    [
        "sample",
        "CHROM",
        "POS",
        "resident_gene_id",
        "REF",
        "allele",
        "AF",
        "AD",
        "DP",
        "MQRS",
        "n_alleles",
    ]
    + [f"snpEff__{xx}" for xx in ann_format_columns]
]

df_filt_long = df_filt_long.replace(np.nan, pd.NA)
df_filt_long = df_filt_long.sort_values(by=["POS", "sample"])
df_filt_long = df_filt_long.rename(columns={"snpEff_Gene_ID": "snpEff_gene_id"})

# %% CHUNK: Write out long annotated snp table
df_filt_long.to_csv(
    outdir.joinpath(snp_table_out_base),
    sep="\t",
    index=False,
)

# %%

