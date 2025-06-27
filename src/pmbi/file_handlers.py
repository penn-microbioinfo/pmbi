from __future__ import annotations

import copy
import io
import os
import re
from pathlib import Path
from typing import Any, Callable, Generator, Iterator
import subprocess

import pandas as pd
from munch import Munch
import toolz
import operator

from pmbi.s3.lib import object_key_list
from pmbi.cellranger.util import get_modality_from_string 
from pmbi.util.misc import get_substring


# %% CHUNK: SheetHandler - bad name as not inheriting from Handler {{{
class SheetHandler(object):
    def __init__(self, path: Path):
        self.path = path
        self.sheet = None
        if self.path.suffix == ".xlsx":
            self.sheet = pd.read_excel(self.path)
        elif self.path.suffix.lower() == ".csv":
            self.sheet = pd.read_csv(self.path, sep=",")
        elif self.path.suffix.lower() == ".tsv":
            self.sheet = pd.read_csv(self.path, sep="\t")
        else:
            raise IOError(f"Unrecognized sheet filename extension {self.path.suffix}")


# }}}


# %% CHUNK: Backends {{{
class Backend(object):
    def __init__(self):
        pass

    def _list_files(self, src: Path) -> list[Path]:
        return []

    def _pull_files(self, src, dest):
        pass

    def _push_files(self, src, dest):
        pass

    def _clean(self):
        pass


class LocalBackend(Backend):
    def __init__(self):
        super().__init__()

    def _list_files(self, src: Path) -> list[Path]:
        return [src.joinpath(f) for f in src.iterdir()]

    def _pull_files(self, src, dest):
        for f in src:
            os.symlink(f, os.path.join(dest, os.path.basename(f)))


class S3Backend(Backend):
    def __init__(self):
        super().__init__()

    def _list_files(self, src: Path) -> list[Path]:
        src_str = str(src)
        return object_key_list(src_str)


# }}}


# %% CHUNK: Handler def {{{
class Handler(object):
    def __init__(
        self,
        path: Path,
        pattern="[_]R[0-9][_].+[.]fastq[.]gz$",
        backend: Backend = LocalBackend(),
    ):
        self.s3_pattern = "^s3[:][/]{2}"
        self.backend = backend

        # DEPRECATED: For typing problem of converting str to Handler If no backend is given, try to auto-detect which one to use {{{
        #     if backend is None:
        #         if re.search(self.s3_pattern, path_str):
        #             self.backend = S3Backend()
        #         else:
        #             self.backend = LocalBackend()
        #     # Otherwise, assign by name
        #     else:
        #         if backend in self._supported_backends():
        #             self.backend = self._supported_backends()[backend]()
        #         else:
        #             raise NotImplementedError(f"Unsupported backend specified: {backend}")
        # # }}}

        self.filelist: list[Path] = self.backend._list_files(src=path)

        self._filter_filelist(pattern)

    def _filter_filelist(self, pattern: str) -> None:
        filtered = []
        for f in self.filelist:
            fb = os.path.basename(f)
            if re.search(pattern, fb) is not None:
                filtered.append(f)
        self.filelist = filtered

    def _files(self) -> Iterator[Path]:
        for f in self.filelist:
            yield f

    # DEPRECATED: with above, re: typing concerns {{{
    # @staticmethod
    # def _supported_backends() -> dict[str,:
    #     return {
    #         "local": LocalBackend,
    #         "s3": S3Backend,
    #     }
    # }}}


# }}}


# %% CHUNK: CellrangerHandler def {{{
class CellrangerHandler(Handler):
    def __init__(
        self,
        path: Path,
        config: Munch,
        pattern: str = "[_]R[0-9][_].+[.]fastq[.]gz$",
        backend: Backend = LocalBackend(),
    ):

        super().__init__(
            path=path,
            pattern=pattern,
            backend=backend,
        )

        self.config = config
        self.table = self._make_table()

    def _feature_type_converter(self):
        return self._modality_metadata().set_index("name")[
            ["cellranger_multi_feature_type"]
        ]

    def _modality_metadata(self):
        return pd.DataFrame(map(lambda x: x.__dict__, self.config.modalities))

    def _make_table(self):
        """
        Constructs a DataFrame containing the information needed to generate a cellranger multi config csv
        from the list of files processed by the handler. Regex patterns from config are extracted from the
        file basenames or made from combinations of extracted strings.

        Modalities are converted to celranger
        "feature_types" based on CellrangerHandler.default_feature_type_converter().


        Returns:
            pd.DataFrame: the constructed DataFrame
        """
        ldict = []
        # TODO: These filename operations should be vectorized after the DataFrame is constructed
        for f in self._files():
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
                    "backend": self.backend,
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

    def subset(self, f: Callable) -> CellrangerHandler:
        sub = copy.deepcopy(self)
        sub.table = sub.table[f(sub.table)]
        sub.filelist = [f for f in self.filelist if f in sub.table["read_path"].values]
        return sub

    def split(self, by: str) -> Generator[CellrangerHandler, None, None]:
        for by_val in self.table[by].unique():
            yield self.subset(lambda x: x[by] == by_val)


class CellrangerMultiConfigCsv(object):
    def __init__(self, handler: CellrangerHandler):
        self.handler = handler
        self.libraries = self._libraries_section()

    def _section_header_converter(self):
        return self.handler._modality_metadata().set_index("name")[
            ["cellranger_multi_csv_section_header"]
        ]

    def _libraries_section(self):
        libraries = pd.DataFrame(
            {
                "fastq_id": self.handler.table["sample"].str.cat(
                    others=[
                        self.handler.table["modality"],
                        self.handler.table["technical_replicate"],
                    ],
                    sep="_",
                ),
                "fastqs": self.handler.table["read_dir"],
                "feature_types": self.handler.table["cellranger_multi_feature_type"],
            }
        )
        assert (
            libraries["fastq_id"].unique().shape
            == libraries["feature_types"].unique().shape
        ), "Number of unique fastq_ids should be equal to number of modalities"

        return libraries


    # %% TODO: There just has to be a cleaner way of doing this
    def _feature_type_sections(self):
        sections = (
            self.handler._modality_metadata()
            .set_index("cellranger_multi_feature_type")
            .loc[
                self.handler.table["cellranger_multi_feature_type"]
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
        with open(path, "w") as out:
            out.write("[libraries]\n")
            self.libraries.drop_duplicates().to_csv(out, index=False)
            out.write("\n")
            out.write(self._feature_type_sections())


# }}}


class CellrangerMulti(object):
    def __init__(self, handler, wd=Path("."), **kwargs):
        self.handler = handler
        self.csv = CellrangerMultiConfigCsv(handler)
        self.id = self.csv.libraries["fastq_id"].unique()[0]
        self.csv_path = wd.joinpath(f"{self.id}__cr_multi_config.csv")
        self.csv.write(self.csv_path)
        self.more_args = list(toolz.concat([[f"--{k}", v] for k,v in kwargs.items()]))

    def cmd(self):
        return list(map(str, (["cellranger", "multi", "--id", self.id, "--csv", self.csv_path] + self.more_args)))

    def run(self):
        p = subprocess.run(self.cmd(), capture_output=True)
        return p
        


# %%
