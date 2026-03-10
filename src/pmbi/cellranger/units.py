from __future__ import annotations

import io
import tempfile
from pathlib import Path

import pandas as pd
import toolz
from munch import Munch


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
