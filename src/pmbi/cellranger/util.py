import pandas as pd

# %% CHUNK: Read in 10X index sheets {{{
def read_10x_index_sheets(*args, workflow):
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
    ]
    return tenx_idx

     # }}}
