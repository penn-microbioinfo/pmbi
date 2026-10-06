import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.style import library
import natsort

from pmbi.cellranger.util import (
    read_10x_chromium_dual_index_sheets,
    read_10x_chromium_single_index_sheets,
)
from pmbi.config import import_config
from pmbi.file_handlers import SheetHandler

if __name__ == "__main__":
    # %% CHUNK: Argument parser {{{
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "tenx_index_sheets",
        nargs="*",
        help="Paths to 10X index sheets including all indices used.",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        help="Your metadata sheet with atleast two columns, one with library name and one with 10X index name.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="./PmbiBcl2fastqSampleSheet.csv",
        help="Path to write the sample sheet out to.",
    )
    parser.add_argument(
        "-i",
        "--index_design",
        choices=["single", "dual"],
        default="dual",
        help="Design of the illumina indexing: `single` or `dual`",
    )
    dual_index_args = parser.add_argument_group(title="Dual Index Design")
    dual_index_args.add_argument(
        "-w",
        "--workflow",
        default=None,
        help="Workflow that you used - for choosing with i5 index to use.",
    )
    parser.add_argument(
        "-c",
        "--config",
        action="store",
        required=True,
        default="./config.toml",
        help="Path to pmbi.cellranger config file.",
    )
    parser.add_argument(
        "--filter_by",
        action="store",
        required=False,
        help="Column name of user-supplied metadata to apply filter by.",
    )
    parser.add_argument(
        "--filter_vals",
        action="store",
        required=False,
        help="Comma separated list of values allowed in filter_by column",
    )
    args = parser.parse_args()

    # Handle issues with passed values
    if args.index_design == "single" and args.workflow is not None:
        raise ValueError("`Single index design and `--workflow` are incompatible`")

    if args.index_design == "dual" and args.workflow is None:
        raise ValueError("`--workflow` must be specified with dual index designs")

    config_path = Path(args.config)
    if not config_path.exists():
        raise ValueError(f"Config file not found: {config_path}")

    sheet_path = Path(args.metadata)
    if not sheet_path.exists():
        raise ValueError(f"Metadata sheet not found: {sheet_path}")

    # }}}
    config = import_config(args.config)
    metadata = SheetHandler(sheet_path).sheet

    library_name_col, index_name_col = (
        config.metadata.colnames.library_name,
        config.metadata.colnames.index_name,
    )
    req_colnames = pd.Series([library_name_col, index_name_col])
    if not (req_colnames.isin(metadata.columns)).all():
        raise ValueError(
            f"Based on the config file, expected supplied metadata sheet to contain the following columns: {', '.join(req_colnames)}"
        )

    # %% Filter rows of metadata based on args
    if args.filter_by is not None:
        if args.filter_vals is None:
            raise argparse.ValueError(
                message="If filter_by is specified, so must filter_vals"
            )
        else:
            filter_vals = args.filter_vals.split(",")
            print(filter_vals)
            metadata = metadata[metadata[args.filter_by].astype(str).isin(filter_vals)]

    metadata = metadata[req_colnames]

    # Process separately based on index design
    # TODO: Functionalize and impl tests
    if args.index_design == "single":
        tenx_idx = read_10x_chromium_single_index_sheets(*args.tenx_index_sheets)

        merged = pd.merge(
            left=metadata,
            right=tenx_idx,
            how="left",
            left_on=index_name_col,
            right_index=True,
        )

        indices_idx = np.where(~merged.columns.isin(req_colnames))[0]
        samplesheet = (
            merged.melt(
                id_vars=req_colnames,
                value_vars=merged.columns[indices_idx],
                value_name="index",
            )
            .drop(columns=["variable"])
            .dropna()
            .sort_values(library_name_col, key=natsort.natsort_keygen())
            .reset_index(drop=True)
            .filter(items=[library_name_col, "index"])
        )
        libs_without_indices = ~merged[library_name_col].isin(samplesheet[library_name_col])
        if libs_without_indices.any():
            raise ValueError(f"library has no associated indices: {','.join(merged[library_name_col][libs_without_indices])}")

        samplesheet.columns = ["Sample_ID", "index"]

    # TODO: Functionalize and impl tests
    else:  # dual indices
        tenx_idx = read_10x_chromium_dual_index_sheets(
            *args.tenx_index_sheets, workflow=args.workflow
        )
        merged = pd.merge(
            left=metadata,
            right=tenx_idx,
            how="left",
            left_on=index_name_col,
            right_index=True,
        )
        samplesheet = merged[
            [
                library_name_col,
                "i7",
                f"i5_workflow_{args.workflow}",
            ]
        ]
        samplesheet.columns = ["Sample_ID", "index", "index2"]

    with open(args.output, "w") as ss:
        ss.write("[Data]\n")
        ss.write(samplesheet.to_csv(sep=",", index=False))

# }}}
