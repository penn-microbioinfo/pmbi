import argparse
from pathlib import Path

import pandas as pd

from pmbi.config import import_config
from pmbi.file_handlers import SheetHandler

if __name__ == "__main__":
# %% CHUNK: Argument parser {{{
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "tenx_index_sheets",
        nargs="*",
        help="Paths to 10X dual index sheets including all those used.",
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
        "-w",
        "--workflow",
        default="b",
        help="Workflow that you used - for choosing with i5 index to use.",
    )
    parser.add_argument(
        "-c",
        "--config",
        action="store",
        required=False,
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

    # row_filter_group = parser.add_mutually_exclusive_group()
    # row_filter_group.add_argument(
    #     "-d",
    #     "--donors",
    #     action="store",
    #     required=False,
    #     help="Donors to include in output. Default is to include all donors.",
    # )
    # row_filter_group.add_argument(
    #     "-s",
    #     "--samples",
    #     action="store",
    #     required=False,
    #     help="Samples to include in output. Default is to include all samples.",
    # )
    args = parser.parse_args()
# }}}

# %% CHUNK: Read in 10X index sheets {{{
    column_renamer = {
        "index_name": "index_name",
        "index(i7)": "i7",
        "index2_workflow_a(i5)": "i5_workflow_a",
        "index2_workflow_b(i5)": "i5_workflow_b",
    }
    expected_colnames = pd.Index(column_renamer.keys())

    tenx_index_sheets = []
    for path in args.tenx_index_sheets:
        sheet = pd.read_csv(path, sep=",", comment="#")
        if not expected_colnames.isin(sheet.columns).all():
            raise ValueError(f"Input sheet missing expected column names: {path}")
        else:
            sheet = sheet.rename(columns=column_renamer)[pd.Index(column_renamer.values())]
            tenx_index_sheets.append(sheet)

    tenx_idx = pd.concat(tenx_index_sheets, axis=0)[
        ["index_name", "i7", f"i5_workflow_{args.workflow}"]
    ]
     # }}}

# %% CHUNK: Read in user metadata, filter {{{
    config = import_config(args.config)
    metadata = SheetHandler(Path(args.metadata)).sheet
# %% Filter rows of metadata based on args
    if args.filter_by is not None:
        if args.filter_vals is None:
            raise argparse.ArgumentError(argument="filter", message="If filter_by is specified, so must filter_vals")
        else:
            filter_vals = args.filter_vals.split(',')
            print(filter_vals)
            metadata = metadata[metadata[args.filter_by].astype(str).isin(filter_vals)]

# }}}

# %% CHUNK: Merge user-provided 10X index names with their actual index sequences {{{
    metadata = metadata[
        [config.metadata.colnames.library_name, config.metadata.colnames.index_name]
    ]

    merged = pd.merge(
        left=metadata,
        right=tenx_idx,
        how="inner",
        left_on=config.metadata.colnames.index_name,
        right_on="index_name",
    )
# }}}

# %% CHUNK: Write the sample sheet {{{
    samplesheet = merged[
        [config.metadata.colnames.library_name, "i7", f"i5_workflow_{args.workflow}"]
    ]
    # TODO: Put function somewhere to validate these column names prior to demux
    samplesheet.columns = ["Sample_ID", "index", "index2"]

    with open(args.output, "w") as ss:
        ss.write("[Data]\n")
        ss.write(samplesheet.to_csv(sep=",", index=False))

# }}}
