from __future__ import annotations

import io
import os
from itertools import chain
from pathlib import Path
from typing import Optional
import re
import yaml

import natsort
import pandas as pd
from munch import Munch

from pmbi.cellranger.runners import (
    CellrangerArcRunner,
    CellrangerAtacRunner,
    CellrangerMultiRunner,
)
from pmbi.cellranger.units import (
    CellrangerArcUnit,
    CellrangerAtacUnit,
    CellrangerMultiUnit,
)
from pmbi.file_handlers import Backend, LocalBackend
from pmbi.util.misc import get_substring, substring_detected


class CellrangerCollection:
    """Handles the initial collection of all files and metadata for cellranger multi operations"""

    def __init__(
        self,
        path: Path | list[Path],
        config: Munch,
        # include: str = "[_]R[0-9][_][0-9]+[.]fastq[.]gz$",
        # exclude: Optional[str] = None,
        backend: Backend = LocalBackend(),
    ):
        self.path = path
        self.config = config
        self.backend = backend
        # self.include = config.filename_patterns.include
        # self.exclude = config.filename_patterns.exclude

        if isinstance(path, Path):
            self.filelist = self._gather_files(self.path)
        elif isinstance(path, list):
            self.filelist = self._gather_files_multiple_paths(self.path)
        else:
            raise ValueError(f"invalid type for `path` arg: {type(path)}")

        if len(self.filelist) == 0:
            raise ValueError(
                "File list is of length 0. Check your include/exclude patterns."
            )

        if hasattr(self.config.filename_patterns, "exclude") and not hasattr(self.config.filename_patterns, "include"):
            raise ValueError("Cannot configure an `exclude` pattern without an `include` pattern")

        match self.config.run.cellranger_flavor:
            case "arc":
                self.UnitType = CellrangerArcUnit
                self.RunnerType = CellrangerArcRunner
            case "atac":
                self.UnitType = CellrangerAtacUnit
                self.RunnerType = CellrangerAtacRunner
            case "multi":
                self.UnitType = CellrangerMultiUnit
                self.RunnerType = CellrangerMultiRunner
            case _:
                raise ValueError(
                    f"Invalid `run.cellranger_flavor` in config: {self.config.run.cellranger_flavor}"
                )

        self.table = self._make_table()

    def _gather_files_multiple_paths(self, path: list[Path]) -> list[Path]:
        """Gather and filter files from a list[Path] based on include pattern"""
        files = []
        for p in path:
            files.append(self._gather_files(p))
        return list(chain.from_iterable(files))

    def _gather_files(self, path: Path) -> list[Path]:
        """Gather and filter files from a Path based on pattern"""
        files = self.backend._list_files(src=path)
        filtered = []
        for f in files:
            add_file = not hasattr(self.config.filename_patterns, "include")
            fb = os.path.basename(f)
            if hasattr(self.config.filename_patterns, "include"):
                add_file = substring_detected(self.config.filename_patterns.include, fb)
                if hasattr(self.config.filename_patterns, "exclude") and add_file:
                    add_file = not substring_detected(
                        self.config.filename_patterns.exclude, fb
                    )
            if add_file:
                filtered.append(f)
        return filtered

    # def _feature_type_converter(self):
    #     """Convert modalities to cellranger feature types"""
    #     return pmbiconf.table_to_dataframe(self.config.modalities).set_index("name")[
    #         ["cellranger_multi_feature_type"]
    #     ]

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
        df = df.sort_values(
            ["sample", "modality", "read_number"], key=natsort.natsort_keygen()
        ).reset_index(drop=True)
        df["modality"] = pd.Categorical(
            values=df["modality"], categories=df["modality"].unique()
        )
        # df["cellranger_multi_feature_type"] = df["modality"].apply(
        #     lambda m: self._feature_type_converter().loc[m].item()
        # )
        return df

    def get_units(self) -> list[CellrangerUnit]:
        """Split collection into individual sample units"""
        units = []
        for sample in self.samples():
            unit_files = self.table[self.table["sample"] == sample].copy()
            units.append(
                self.UnitType(
                    table=unit_files,
                    config=self.config,
                    files=[Path(p) for p in unit_files["read_path"]],
                )
            )
        return units

    def samples(self):
        return self.table["sample"].unique()

    def modalities_per_sample(self) -> pd.Series:
        return self.table.groupby("sample")["modality"].nunique()

    def to_index_hopping_filter_csv(self) -> io.StringIO:
        buf = io.StringIO()
        csv_df = (
            self.table[["read_number", "sample", "read_path"]]
            .pivot(columns="read_number", values="read_path", index="sample")
            .reset_index()
            .sort_values("sample", key=natsort.natsort_keygen())
            .rename(columns={"sample": "SampleId"})
        )
        read_columns = natsort.natsorted([c for c in csv_df.columns if c != "SampleId"])
        csv_df = csv_df[read_columns + ["SampleId"]]
        csv_df.to_csv(buf, sep=",", index=False)
        buf.seek(0)
        return buf

    def to_scc_proc_sample_config_yaml(self, modalities:list[str]=["ADT", "HTO"], sample_regex:str|None=None):
        if sample_regex is not None:
            self.table["sample"] = self.table["sample"].apply(lambda s: re.search(sample_regex, s).group(0))
        table = self.table[self.table["modality"].isin(modalities)]
        table["modality"] = table["modality"].astype(str)
        samples_list = []
        for s in table["sample"].unique():
            stable = table[table["sample"]==s]
            modal_grouped = stable.groupby("modality").apply(lambda g: [os.fspath(p) for p in g["read_filename"].to_list()])
            modal_grouped.index = pd.Index([f"{m.lower()}_fastqs" for m in modal_grouped.index])
            d = {"name": s}
            d.update(modal_grouped.to_dict())
            samples_list.append(d)

        return yaml.dump({"samples": samples_list}, sort_keys=False)
# %%

