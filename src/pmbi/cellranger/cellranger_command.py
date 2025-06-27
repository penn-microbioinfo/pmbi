from __future__ import annotations

import copy
import io
import os
import re
import shutil
from pathlib import Path
import subprocess
from typing import Generator, Iterator

import pandas as pd
from munch import Munch
import toolz

from pmbi.file_handlers import Backend, LocalBackend
from pmbi.util.misc import get_substring
from pmbi.logging import streamLogger
import pmbi.config as pmbiconf
import pmbi.subproc as pmbiproc


class CellrangerCollection:
    """Handles the initial collection of all files and metadata for cellranger multi operations"""
    
    def __init__(
        self, 
        path: Path, 
        config: Munch, 
        pattern: str = "[_]R[0-9][_].+[.]fastq[.]gz$", 
        backend: Backend = LocalBackend()
    ):
        self.path = path
        self.config = config
        self.backend = backend
        self.pattern = pattern
        self.filelist = self._gather_files()
        self.table = self._make_table()
        
    def _gather_files(self) -> list[Path]:
        """Gather and filter files based on pattern"""
        files = self.backend._list_files(src=self.path)
        filtered = []
        for f in files:
            fb = os.path.basename(f)
            if re.search(self.pattern, fb) is not None:
                filtered.append(f)
        return filtered

    def _feature_type_converter(self):
        """Convert modalities to cellranger feature types"""
        return pmbiconf.table_to_dataframe(self.config.modalities).set_index("name")[
            ["cellranger_multi_feature_type"]
        ]

    def _make_table(self) -> pd.DataFrame:
        """Create a table of file metadata"""
        ldict = []
        for f in self.filelist:
            fb = os.path.basename(f)
            ldict.append(
                {
                    "sample": get_substring(
                        string=fb, pattern=self.config.filename_patterns.sample
                    ),
                    "modality": get_substring(
                        string=fb, pattern=self.config.filename_patterns.modality
                    ),
                    "read_number": get_substring(
                        string=fb, pattern=self.config.filename_patterns.read_number
                    ),
                    "read_filename": fb,
                    "read_path": f,
                    "read_dir": os.path.split(f)[0],
                }
            )
        df = pd.DataFrame(ldict)
        df = df.sort_values(["sample", "modality", "read_number"])
        df["modality"] = pd.Categorical(
            values=df["modality"], categories=df["modality"].unique()
        )
        df["cellranger_multi_feature_type"] = df["modality"].apply(
            lambda m: self._feature_type_converter().loc[m].item()
        )
        return df

    def get_units(self) -> list[CellrangerUnit]:
        """Split collection into individual sample units"""
        units = []
        for sample in self.table["sample"].unique():
            unit_files = self.table[self.table["sample"] == sample].copy()
            units.append(
                CellrangerUnit(
                    table=unit_files,
                    config=self.config,
                    files=[Path(p) for p in unit_files["read_path"]]
                )
            )
        return units

class CellrangerUnit:
    
    def __init__(self, table: pd.DataFrame, config: Munch, files: list[Path]):
        self.table = table
        self.config = config
        self.files = files
        self.sample = self._validate_single_sample()
        
    def _validate_single_sample(self):
        """Ensure unit contains exactly one sample"""
        unique_samples = self.table["sample"].unique()
        if len(unique_samples) != 1:
            raise ValueError(
                f"CellrangerUnit must contain exactly one sample. Found: {unique_samples}"
            )
        return unique_samples[0]
    
    @property
    def config_csv(self) -> CellrangerConfigCsv:
        """Generate config CSV for this unit"""
        return CellrangerConfigCsv(self)
    
    def create_runner(self, wd: Path = Path("."), **kwargs) -> CellrangerRunner:
        """Create a runner for this unit"""
        return CellrangerRunner(self, wd, **kwargs)

class CellrangerConfigCsv:
    """Handles CSV generation for a single cellranger multi unit"""
    
    def __init__(self, unit: CellrangerUnit):
        self.unit = unit
        self.libraries = self._libraries_section()
        
    def _section_header_converter(self):
        """Get section headers from modality metadata"""
        return pmbiconf.table_to_dataframe(self.unit.config.modalities).set_index("name")[
            ["cellranger_multi_csv_section_header"]
        ]

    def _libraries_section(self) -> pd.DataFrame:
        """Generate libraries section of config CSV"""
        libraries = pd.DataFrame(
            {
                "fastq_id": self.unit.table["sample"].str.cat(
                    others=[
                        self.unit.table["modality"],
                    ],
                    sep="_",
                ),
                "fastqs": self.unit.table["read_dir"],
                "feature_types": self.unit.table["cellranger_multi_feature_type"],
            }
        )
        assert (
            libraries["fastq_id"].unique().shape
            == libraries["feature_types"].unique().shape
        ), "Number of unique fastq_ids should be equal to number of modalities"
        return libraries

    def _feature_type_sections(self) -> str:
        """Generate feature type sections of config CSV"""
        sections = (
            pmbiconf.table_to_dataframe(self.unit.config.modalities)
            .set_index("cellranger_multi_feature_type")
            .loc[
                self.unit.table["cellranger_multi_feature_type"]
                .astype(str)
                .unique(),
                :,
            ]
            .set_index("cellranger_multi_csv_section_header", drop=True)
            .drop(columns=["name"])
            .drop_duplicates()
        )
        buffer = io.StringIO()
        for sec in sections.iterrows():
            print(sec)
            buffer.write(f"[{sec[0]}]\n")
            buffer.write(
                "\n".join(
                    [f"{k},{v}" for k, v in sec[1].to_dict().items() if not pd.isna(v)]
                )
            )
            buffer.write("\n\n")
        return buffer.getvalue()

    def to_stringio(self):
        csv = io.StringIO()
        csv.write("[libraries]\n")
        self.libraries.drop_duplicates().to_csv(csv, index=False)
        csv.write("\n")
        csv.write(self._feature_type_sections())
        csv.seek(0)
        return csv

    def write(self, path: Path):
        with open(path, "w") as out:
            out.write(self.to_stringio().read())


class CellrangerRunner:
    
    def __init__(self, unit: CellrangerUnit, wd: Path = Path("."), **kwargs):
        self.unit = unit
        self.wd = wd
        self.id = self.unit.sample
        self.csv_base = Path(f"{self.id}__config.csv")
        self.outs_dest = self.wd / f"{self.id}__outs"
        self.logger = streamLogger("CellrangerRunner")
        self.kwargs = kwargs
        
    def more_args(self) -> list[str]:
        """Convert kwargs to command line arguments"""
        return list(toolz.concat([[f"--{k}", v] for k,v in self.kwargs.items()]))
        
    def cmd(self) -> list[str]:
        return list(map(str, (
            ["cellranger", "multi", "--id", self.id, "--csv", self.csv_base] 
                + self.more_args()
        )))

    def _check_output_exists(self):
        if self.outs_dest.exists():
            mes = f"Output directory exists. Skipping myself: {self.id}"
            self.logger.critical(mes)
            raise OSError(mes)

    def cleanup(self):
        outs_src = self.wd / self.id / "outs"
        outs_dest = self.outs_dest
        shutil.copytree(outs_src, outs_dest)
        shutil.rmtree(self.wd / self.id)
        self.csv_base.unlink()
        
    def run(self):
        logger = streamLogger("cellranger_multi")
        logger.info(f"Changing directory: {self.wd}")
        os.chdir(self.wd)
        logger.info(f"Writing config csv to: {self.wd / self.csv_base}")
        self.unit.config_csv.write(self.csv_base)
        pmbiproc.run_and_log(cmd=self.cmd(), logger=logger)
        self.cleanup()
