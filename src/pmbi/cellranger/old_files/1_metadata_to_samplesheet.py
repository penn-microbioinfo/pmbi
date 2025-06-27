import argparse
import sys

import pandas as pd

from pmbi.config import import_config

parser = argparse.ArgumentParser()
parser.add_argument(
    "-m", "--metadata", action="store", required=True, help="Path to metadata sheet."
)
parser.add_argument(
    "-o", "--output", action="store", required=True, help="Output filename."
)
parser.add_argument(
    "-c",
    "--config",
    action="store",
    required=False,
    default="./config.toml",
    help="Path to pmbi.cellranger config file.",
)
sample_filter_group = parser.add_mutually_exclusive_group()
sample_filter_group.add_argument(
    "-d",
    "--donors",
    action="store",
    required=False,
    help="Donors to include in output. Default is to include all donors.",
)
sample_filter_group.add_argument(
    "-s",
    "--samples",
    action="store",
    required=False,
    help="Samples to include in output. Default is to include all samples.",
)
args = parser.parse_args()

# %% CHUNK: Read and parse config {{{
config = import_config(args.config)
print(config.metadata.colnames.sample_id)
# }}}

md = pd.read_csv(args.metadata, dtype=str)
if args.donors is not None or args.samples is not None:
    if args.donors is not None:
        donors = args.donors.strip().split(",")
        md = md[md[config.metadata.colnames.donor_id].isin(donors)]
    if args.samples is not None:
        samples = args.samples.strip().split(",")
        md = md[md[config.metadata.colnames.sample_id].isin(samples)]

md["Lane"] = "*"
md["Sample"] = (
    md[config.metadata.colnames.sample_id]
    + "_"
    + md[config.metadata.colnames.modality]
    + "_"
    + md[config.metadata.colnames.chip_well]
)
md["Index"] = (
    "SI"
    + "-"
    + md[config.metadata.colnames.index_plate]
    + "-"
    + md[config.metadata.colnames.well]
)
ss = md[["Lane", "Sample", "Index"]]

ss.to_csv(args.output, index=False, header=True)
