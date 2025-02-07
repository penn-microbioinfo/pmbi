import os

import pandas as pd
import pytest

from pmbi.config import import_config
from pmbi.file_handlers import CellrangerHandler

CONTENT = "@a_read"


# %% CHUNK: Fixtures {{{
@pytest.fixture(scope="session")
def fastq_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("fastq")
    for f in range(0, 8):
        for m in CellrangerHandler.accepted_modalities():
            for r in ["R1", "R2"]:
                ff = d / f"HPAP-173_CC_{m}_{f}_S1_L001_{r}_001.fastq.gz"
                ff.write_text(CONTENT, encoding="utf-8")
    extra_file = d / "extra_something.txt"
    extra_file.write_text("extraneous", encoding="utf-8")
    return d


@pytest.fixture
def cellranger_handler(fastq_dir):
    config = import_config("/home/ubuntu/dev/pmbi/src/pmbi/cellranger/config.toml")
    return CellrangerHandler(path=fastq_dir, config=config)


@pytest.fixture
def expected_columns():
    return pd.Index(
        [
            "sample",
            "modality",
            "technical_rep",
            "read_number",
            "read_filename",
            "read_path",
            "read_dir",
            "backend",
            "sample_rep",
        ]
    )


@pytest.fixture
def n_fastq(fastq_dir):
    return len([f for f in fastq_dir.iterdir() if f.name.endswith("fastq.gz")])


# }}}

# %% CHUNK: Tests {{{


# %% TEST: Test that the CellrangerHandler detects any files at all in the fastq dir
def test_handler_detects_files(cellranger_handler):
    assert len(cellranger_handler.filelist) > 0


# %% TEST: Test that the CellrangerHandler correctly filters out files that don't match the supplied pattern
def test_handler_filters_files(cellranger_handler):
    assert len(cellranger_handler.filelist) > 0
    assert "extra_something.txt" not in [
        os.path.basename(f) for f in cellranger_handler.filelist
    ]

# %% TEST: Test the fildelity of the CellrangerHandler table
def test_handler_table(cellranger_handler, expected_columns, n_fastq):
    assert expected_columns.isin(cellranger_handler.table.columns).all()
    assert len(expected_columns) == len(cellranger_handler.table.columns)
    assert cellranger_handler.table.shape == (n_fastq, len(expected_columns))

# %% TEST: Test that the handler can subset itself
def test_handler_subset(cellranger_handler, expected_columns, n_fastq):
    sub = cellranger_handler.subset(f=lambda r: r["read_number"] == "R1")
    assert sub.table.shape == (n_fastq / 2, len(expected_columns))
    assert len(sub.filelist) == n_fastq / 2


# }}}
