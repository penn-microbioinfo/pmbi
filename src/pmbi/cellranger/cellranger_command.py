from __future__ import annotations

import copy
import io
import os
import re
from pathlib import Path
import subprocess
from typing import Generator, Iterator

import pandas as pd
from munch import Munch
import toolz

from pmbi.file_handlers import Backend, LocalBackend
from pmbi.util import get_substring

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
        return self._modality_metadata().set_index("name")[
            ["cellranger_multi_feature_type"]
        ]

    def _modality_metadata(self):
        """Get modality metadata from config"""
        return pd.DataFrame(map(lambda x: x.__dict__, self.config.modalities))

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
                    "technical_replicate": get_substring(
                        string=fb,
                        pattern=self.config.filename_patterns.technical_replicate,
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
        df["sample_rep"] = df["sample"] + "_" + df["technical_replicate"]
        df = df.sort_values(["sample_rep", "modality", "read_number"])
        df["modality"] = pd.Categorical(
            values=df["modality"], categories=df["modality"].unique()
        )
        df["cellranger_multi_feature_type"] = df["modality"].apply(
            lambda m: self._feature_type_converter().loc[m].item()
        )
        return df

    def get_units(self) -> list[CellrangerUnit]:
        """Split collection into individual sample_rep units"""
        units = []
        for sample_rep in self.table["sample_rep"].unique():
            unit_files = self.table[self.table["sample_rep"] == sample_rep].copy()
            units.append(
                CellrangerUnit(
                    table=unit_files,
                    config=self.config,
                    files=[Path(p) for p in unit_files["read_path"]]
                )
            )
        return units

class CellrangerUnit:
    """Represents a single cellranger multi operation"""
    
    def __init__(self, table: pd.DataFrame, config: Munch, files: list[Path]):
        self.table = table
        self.config = config
        self.files = files
        self._validate_single_sample_rep()
        
    def _validate_single_sample_rep(self):
        """Ensure unit contains exactly one sample_rep"""
        unique_sample_reps = self.table["sample_rep"].unique()
        if len(unique_sample_reps) != 1:
            raise ValueError(
                f"CellrangerUnit must contain exactly one sample_rep. Found: {unique_sample_reps}"
            )
        self.sample_rep = unique_sample_reps[0]
    
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
        return self.unit._modality_metadata().set_index("name")[
            ["cellranger_multi_csv_section_header"]
        ]

    def _libraries_section(self) -> pd.DataFrame:
        """Generate libraries section of config CSV"""
        libraries = pd.DataFrame(
            {
                "fastq_id": self.unit.table["sample"].str.cat(
                    others=[
                        self.unit.table["modality"],
                        self.unit.table["technical_replicate"],
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
            self.unit._modality_metadata()
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
            .iterrows()
        )
        buffer = io.StringIO()
        for sec in sections:
            buffer.write(f"[{sec[0]}]\n")
            buffer.write(
                "\n".join(
                    [f"{k},{v}" for k, v in sec[1].to_dict().items() if not pd.isna(v)]
                )
            )
            buffer.write("\n\n")
        return buffer.getvalue()

    def write(self, path: Path):
        """Write config CSV to file"""
        with open(path, "w") as out:
            out.write("[libraries]\n")
            self.libraries.drop_duplicates().to_csv(out, index=False)
            out.write("\n")
            out.write(self._feature_type_sections())

class CellrangerRunner:
    """Handles command generation and execution for a single cellranger multi unit"""
    
    def __init__(self, unit: CellrangerUnit, wd: Path = Path("."), **kwargs):
        self.unit = unit
        self.wd = wd
        self.id = self.unit.sample_rep
        self.csv_path = self.wd.joinpath(f"{self.id}__cr_multi_config.csv")
        self.kwargs = kwargs
        
    @property
    def more_args(self) -> list[str]:
        """Convert kwargs to command line arguments"""
        return list(toolz.concat([[f"--{k}", v] for k,v in self.kwargs.items()]))
        
    def cmd(self) -> list[str]:
        """Generate cellranger multi command"""
        return list(map(str, (
            ["cellranger", "multi", "--id", self.id, "--csv", self.csv_path] 
            + self.more_args
        )))
        
    def run(self) -> subprocess.CompletedProcess:
        """Execute cellranger multi command"""
        # First write the config CSV
        self.unit.config_csv.write(self.csv_path)
        # Then run the command
        return subprocess.run(self.cmd(), capture_output=True)
