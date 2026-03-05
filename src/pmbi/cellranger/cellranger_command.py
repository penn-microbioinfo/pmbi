from __future__ import annotations

import copy
import io
import os
import re
import shutil
import subprocess
import tempfile
from itertools import chain
from pathlib import Path
from typing import TYPE_CHECKING, Generator, Iterator, Optional

import natsort
import pandas as pd
import toolz
from munch import Munch

import pmbi.config as pmbiconf
import pmbi.subproc as pmbiproc
from pmbi.file_handlers import Backend, LocalBackend
from pmbi.logging import streamLogger
from pmbi.util.misc import get_substring, substring_detected

# if TYPE_CHECKING:
#     from pmbi.cellranger.runners import CellrangerRunner, CellrangerArcRunner, CellrangerMultiRunner


class CellrangerCollection:
    """Handles the initial collection of all files and metadata for cellranger multi operations"""

    def __init__(
        self,
        path: Path | list[Path],
        config: Munch,
        include: str = "[_]R[0-9][_][0-9]+[.]fastq[.]gz$",
        exclude: Optional[str] = None,
        backend: Backend = LocalBackend(),
    ):
        self.path = path
        self.config = config
        self.backend = backend
        self.include = include
        self.exclude = exclude

        if isinstance(path, Path):
            self.filelist = self._gather_files(self.path)
        elif isinstance(path, list):
            self.filelist = self._gather_files_multiple_paths(self.path)
        else:
            raise ValueError("invalid type for `path` arg")

        if len(self.filelist) == 0:
            raise ValueError(
                "File list is of length 0. Check your include/exclude patterns."
            )

        if self.config.run.cellranger_flavor == "arc":
            self.UnitType = CellrangerArcUnit
        elif self.config.run.cellranger_flavor == "atac":
            self.UnitType = CellrangerAtacUnit
        elif self.config.run.cellranger_flavor == "multi":
            self.UnitType = CellrangerMultiUnit
        else:
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
            fb = os.path.basename(f)
            if substring_detected(self.include, fb):
                if self.exclude is not None:
                    if not substring_detected(self.exclude, fb):
                        filtered.append(f)
                else:
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
        print(df)
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
        csv_df = self.table[["read_number", "sample", "read_path"]].pivot(
            columns="read_number", values="read_path", index="sample"
        ).reset_index().sort_values("sample", key=natsort.natsort_keygen()).rename(
            columns={"sample": "SampleId"}
        )
        read_columns = natsort.natsorted([c for c in csv_df.columns if c!="SampleId"])
        csv_df=csv_df[read_columns+["SampleId"]]
        csv_df.to_csv(buf, sep=",", index=False)
        buf.seek(0)
        return buf


class CellrangerUnit:
    def __init__(self, table: pd.DataFrame, config: Munch, files: list[Path]):
        self.table = table
        self.config = config
        self.files = files
        self.sample = self._validate_single_sample()
        self.modalities = {m.name: m for m in self.config.modalities}
        self.samples = None
        self.SUPPORTED_MODALITIES = []
        self.REQUIRED_MODALITIES = []
        self._tempfile_paths = {}

    def _validate_single_sample(self):
        """Ensure unit contains exactly one sample"""
        unique_samples = self.table["sample"].unique()
        if len(unique_samples) != 1:
            raise ValueError(
                f"CellrangerUnit must contain exactly one sample. Found: {unique_samples}"
            )
        return unique_samples[0]

    @property
    def config_csv(self) -> io.StringIO:
        return io.StringIO()

    def create_runner(self, wd: Path = Path("."), **kwargs):
        return CellrangerRunner(self, wd, **kwargs)

    # TODO: Functions like this should return bool as welll as info about the results of specefic tests
    def _verify_modalities(self):
        all_modalities_in_supported = (
            self.table["modality"].isin(self.SUPPORTED_MODALITIES).all().item()
        )
        all_req_represented_in_table = [
            m in self.table["modality"].values for m in self.REQUIRED_MODALITIES
        ]
        all_req_represented_in_config = all(
            [m in self.modalities for m in self.REQUIRED_MODALITIES]
        )
        config_and_table_agree = all(
            [m in self.table["modality"].values for m in self.modalities]
        )
        if (
            all_modalities_in_supported
            and all_req_represented_in_table
            and all_req_represented_in_config
            and config_and_table_agree
        ):
            return True
        else:
            print(
                all_modalities_in_supported,
                all_req_represented_in_table,
                all_req_represented_in_config,
                config_and_table_agree,
            )
            return False


class CellrangerMultiUnit(CellrangerUnit):

    def __init__(self, table: pd.DataFrame, config: Munch, files: list[Path]):
        super().__init__(table, config, files)
        self.CSV_LIBRARY_TYPES = {
            "GEX": "Gene Expression",
            "ADT": "Antibody Capture",
            "HTO": "Antibody Capture",
            "VDJ-T": "VDJ-T",
            "VDJ-B": "VDJ-B",
        }
        self.CSV_SECTION_HEADERS = {
            "GEX": "gene-expression",
            "ADT": "feature",
            "HTO": "feature",
            "VDJ-T": "vdj",
            "VDJ-B": "vdj",
        }
        self.SUPPORTED_MODALITIES = list(self.CSV_LIBRARY_TYPES.keys())
        self.REQUIRED_MODALITIES = ["GEX"]

        if not self._verify_modalities():
            raise ValueError("Unable to verify modalities for `cellranger multi` unit")

        # Check to see if a combined antibody reference is necessary (ie antibody catalog + HTO catalog)
        if "ADT" in self.modalities and "HTO" in self.modalities:
            pass
            # self._combined_adt_hto_ref()

    def _combined_adt_hto_ref(self, modality_keys=["ADT", "HTO"]):
        req_columns = pd.Series(
            ["id", "name", "read", "pattern", "sequence", "feature_type"]
        )
        reference_paths = {
            k: v.reference for k, v in self.modalities.items() if k in modality_keys
        }
        ref_dfs = []
        for k, rp in reference_paths.items():
            ref_df = pd.read_csv(rp, sep=",", header=0)
            if k == "HTO":
                self.samples = ref_df[["name", "id"]]
                self.samples.columns = ["sample_id", "hashtag_ids"]
            if not (ref_df.columns == req_columns).all():
                raise ValueError(
                    f"Antibody reference CSVs must have the columns: {','.join(req_columns)}"
                )
            else:
                ref_dfs.append(ref_df)

        with tempfile.NamedTemporaryFile(
            "w", delete=False, delete_on_close=False
        ) as th:
            self._tempfile_paths["adt_hto_combined_ref"] = th.name
            pd.concat(ref_dfs, axis=0).to_csv(th, sep=",", index=False)
            for mk in modality_keys:
                self.modalities[mk].reference = th.name

    def config_csv(self) -> io.StringIO:
        csv = io.StringIO()
        csv.write("[libraries]\n")
        self._csv_libraries_section().to_csv(csv, sep=",", index=False)
        csv.write("\n")
        for k, v in self._csv_feature_type_sections().items():
            csv.write(f"[{k}]\n")
            v.to_csv(csv, sep=",", index=False, header=False)
            csv.write("\n")

        if self.samples is not None:
            csv.write("[samples]\n")
            self.samples.to_csv(csv, sep=",", index=False)
            csv.write("\n")

        csv.seek(0)
        return csv

    def _csv_libraries_section(self) -> pd.DataFrame:
        """Generate libraries section of config CSV"""
        libraries = (
            pd.DataFrame(
                {
                    "fastq_id": self.table["sample"].str.cat(
                        others=[
                            self.table["modality"],
                        ],
                        sep="_",
                    ),
                    "fastqs": self.table["read_dir"],
                    "feature_types": self.table["modality"].apply(
                        lambda m: self.CSV_LIBRARY_TYPES[m]
                    ),
                }
            )
            .drop_duplicates()
            .reset_index(drop=True)
        )
        return libraries

    def _csv_feature_type_sections(self) -> dict[str, pd.DataFrame]:
        section_headers = set([self.CSV_SECTION_HEADERS[m] for m in self.modalities])
        sections = {}
        skip_keys = ["name"]
        for s in section_headers:
            section_values = {}
            mods = toolz.valfilter(lambda v: v == s, self.CSV_SECTION_HEADERS)

            # Only allow multi-modality feature sections in the case Of ADT+HTO
            # INFO: This will need to be changed if I ever run VDJ-T+VDJ+B
            if len(mods) > 1 and (not all([k in ["ADT", "HTO"] for k in mods])):
                raise ValueError(
                    "Feature type section for >1 modality that is not ADT+HTO. If you are running VDJ-B+VDJ-T, then you need to account for that."
                )
            # unique_keys = toolz.unique(chain([x.keys() for x in mods.
            mod_config_df = pd.DataFrame([self.modalities[m].__dict__ for m in mods])
            for k in mod_config_df.columns:
                if k in skip_keys:
                    continue
                else:
                    if mod_config_df[k].drop_duplicates().shape[0] > 1:
                        raise ValueError(
                            f"For section `{s}`, configuration values that must merge have the same key but different values."
                        )
                    else:
                        val = mod_config_df[k].astype(str).to_list()[0]
                        if val.lower() in ["true", "false"]:
                            val = val.lower()
                        section_values[k] = val

            sections[s] = pd.DataFrame(
                {"key": section_values.keys(), "value": section_values.values()}
            )

        return sections


class CellrangerAtacUnit(CellrangerUnit):
    def __init__(self, table: pd.DataFrame, config: Munch, files: list[Path]):
        super().__init__(table, config, files)
        self.CSV_LIBRARY_TYPES = {
            "ATAC": "Chromatin Accessibility",
        }
        self.SUPPORTED_MODALITIES = list(self.CSV_LIBRARY_TYPES.keys())
        self.REQUIRED_MODALITIES = self.SUPPORTED_MODALITIES
        if not self._verify_modalities():
            raise ValueError("cellranger-atac requires, exclusively, ATAC modalities")

        if not hasattr(self.modalities["ATAC"], "reference"):
            raise ValueError(
                "cellranger-atac config is missing the `modalities.ATAC.reference` configuration option"
            )
        if not hasattr(self.modalities["ATAC"], "fastqs"):
            raise ValueError(
                "cellranger-atac config is missing the required `modalities.ATAC.fastqs` configuration option"
            )

        self.reference = self.modalities["ATAC"].reference
        self.fastqs = self.modalities["ATAC"].fastqs


class CellrangerArcUnit(CellrangerUnit):
    def __init__(self, table: pd.DataFrame, config: Munch, files: list[Path]):
        super().__init__(table, config, files)
        self.CSV_LIBRARY_TYPES = {
            "ATAC": "Chromatin Accessibility",
            "GEX": "Gene Expression",
        }
        self.SUPPORTED_MODALITIES = list(self.CSV_LIBRARY_TYPES.keys())
        self.REQUIRED_MODALITIES = self.SUPPORTED_MODALITIES
        if not self._verify_modalities():
            raise ValueError(
                "cellranger-arc requires, exclusively, both GEX and ATAC modalities"
            )

        if not hasattr(self.config.run, "reference"):
            raise ValueError(
                "cellranger-arc config missing run.reference config option"
            )
        else:
            self.reference = self.config.run.reference

    def config_csv(self) -> io.StringIO:
        csv = io.StringIO()
        csv.write("fastqs,sample,library_type\n")
        for modality, csv_library_type in self.CSV_LIBRARY_TYPES.items():
            fastq_dirs = self.table[self.table["modality"] == modality]["read_dir"]
            if len(fastq_dirs.unique()) != 1:
                raise ValueError(f"Single read_dir not True for  modality: {modality}")
            fastq_dir = fastq_dirs.unique().item()
            csv.write(f"{fastq_dir},{self.sample}_{modality},{csv_library_type}\n")
        csv.seek(0)
        return csv


# %%
