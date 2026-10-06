from _pytest import monkeypatch
from pmbi.cellranger.collections import CellrangerCollection
from pmbi.cellranger.runners import CellrangerAtacRunner
from pmbi.subproc import run_and_log
from pmbi.logging import streamLogger
from pmbi.illumina.fastq import LANE_SPLIT_FASTQ_NO_INDEX, filter_file_paths
from pathlib import Path
import pytest
import shutil
from munch import Munch
import os
import re

@pytest.fixture
def fastq_dir(tmp_path, runopts):
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()

    for fq in runopts["fastqs"].iterdir():
        if re.search(LANE_SPLIT_FASTQ_NO_INDEX, fq.name) is not None:
            link_path = fastq_dir / fq.name.replace("testfastq", "testfastq_ATAC")
            link_path.symlink_to(fq)

    return fastq_dir 

# @pytest.fixture
# def handle_cr_wd(tmp_path, runopts):
#     cr_wd = tmp_path / "cr_wd"
#     cr_wd.mkdir()
#
#     yield cr_wd
#
#     # if cr_wd.exists():
#     #     shutil.rmtree(cr_wd)

@pytest.fixture
def config(runopts, tmp_path):
    config_d = {
        "modalities": [
            {
                "name": "ATAC",
                "fastqs": (tmp_path / "fastqs"),
                "reference": runopts["reference"]
            }
        ],
        "filename_patterns": {
            "sample": "^([^_]+)[_]",
            "modality": "[_](ATAC)[_]",
            "illumina_bits": LANE_SPLIT_FASTQ_NO_INDEX,
            "read_number": "_(R[0-9])_",
            "include": "fastq.gz",
        },
        "run": {
            "wd": tmp_path,
            "cellranger_flavor": "atac",
        },
        "command_line": {
                "localcores": 1,
                "localmem": 16
        }
    }
    return Munch.fromDict(config_d)


def test_cellranger_atac_run(tmp_path, runopts, monkeypatch, config, fastq_dir):
    coll = CellrangerCollection(path=fastq_dir, config=config)
    units = coll.get_units()
    assert len(units)==1
    runner = coll.RunnerType(unit=units[0], wd=config.run.wd, **config.command_line.__dict__)
    monkeypatch.chdir(tmp_path)
    runner.run()
    # cmd = [
    #     runopts["binary"],
    #     "count",
    #     "--id",
    #     runopts["id"],
    #     "--reference",
    #     runopts["reference"],
    #     "--fastqs",
    #     fastq_dir,
    #     "--sample",
    #     runopts["sample"]
    # ]
    # run_and_log(cmd, streamLogger("cellranger-atac-test"))

