import pandas as pd

# %% FUNC: This function tries to pull an accepted modality from a an input string based on an input pattern {{{
def get_modality_from_string(
    string: str,
    pattern: str,
    accepted_modalities: dict[str, str],
) -> str:
    modality = get_substring(string, pattern)
    if modality in accepted_modalities:
        return accepted_modalities[modality]
    else:
        raise ValueError(f"Unexpected modality: `{modality}`")


# %% FUNC: Read in 10X Chromium dual index sheets {{{
def read_10x_chromium_dual_index_sheets(*args, workflow):
    column_renamer = {
        "index_name": "index_name",
        "index(i7)": "i7",
        "index2_workflow_a(i5)": "i5_workflow_a",
        "index2_workflow_b(i5)": "i5_workflow_b",
    }
    expected_colnames = pd.Index(column_renamer.keys())

    if workflow not in ["a", "b"]:
        raise ValueError(f"Invalid workflow: {workflow}")

    tenx_index_sheets = []
    for path in args:
        sheet = pd.read_csv(path, sep=",", comment="#")
        if not expected_colnames.isin(sheet.columns).all():
            raise ValueError(f"Input sheet missing expected column names: {path}")
        else:
            sheet = sheet.rename(columns=column_renamer)[pd.Index(column_renamer.values())]
            tenx_index_sheets.append(sheet)

    tenx_idx = pd.concat(tenx_index_sheets, axis=0)[
        ["index_name", "i7", f"i5_workflow_{workflow}"]
    ].set_indx("index_name", keep=False)
    return tenx_idx

# %% FUNC: Read in 10X Chromium single index sheets {{{
# Files must be in CSV format
# Files must have no column names
# Rows should have the name of the 10x index in the first column
# Index names must correspond to values in the `index_name` column referenced in pmbi/cellranger/1_make_bcl2fastq_samplesheet.py
# Additional columns (n>=1) should contain all index sequences to associate with the index name in the first column
def read_10x_chromium_single_index_sheets(*sheet_paths):
    sheets = []
    for sp in sheet_paths:
        s = pd.read_csv(sp, sep=",", header=None).set_index(0)
        sheets.append(s)

    return pd.concat(sheets)

