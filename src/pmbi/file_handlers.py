import copy
import os
import re
from pathlib import Path

import pandas as pd

from pmbi.s3.lib import object_key_list
from pmbi.util import get_modality_from_string, get_substring


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
        backend: str | None = None,
    ):
        self.s3_pattern = "^s3[:][/]{2}"
        path_str = str(path)

        # If no backend is given, try to auto-detect which one to use
        if backend is None:
            if re.search(self.s3_pattern, path_str):
                self.backend = S3Backend()
            else:
                self.backend = LocalBackend()
        # Otherwise, assign by name
        else:
            if backend in self._supported_backends():
                backend = self._supported_backends()[backend]()
            else:
                raise NotImplementedError(f"Unsupported backend specified: {backend}")

        self.filelist = self.backend._list_files(src=path)

        self._filter_filelist(pattern)

    def _filter_filelist(self, pattern):
        filtered = []
        for f in self.filelist:
            fb = os.path.basename(f)
            if re.search(pattern, fb) is not None:
                filtered.append(f)
        self.filelist = filtered

    def _files(self):
        for f in self.filelist:
            yield f

    @staticmethod
    def _supported_backends():
        return {
            "local": LocalBackend,
            "s3": S3Backend,
        }


# }}}


# %% CHUNK: CellrangerHandler def {{{
class CellrangerHandler(Handler):
    def __init__(
        self,
        path,
        config,
        pattern="[_]R[0-9][_].+[.]fastq[.]gz$",
        backend=None,
    ):

        super().__init__(
            path=path,
            pattern=pattern,
            backend=backend,
        )

        self.config = config
        self._make_table()

    @staticmethod
    def accepted_modalities():
        return {
            "RNA": "RNA",
            "ADT": "ADT",
            "HTO": "HTO",
            "VDJ": "VDJ",
            "VDJ-T": "VDJ-T",
            "VDJ-B": "VDJ-B",
        }

    def _make_table(self):
        ldict = []
        for f in self._files():
            fb = os.path.basename(f)
            ldict.append(
                {
                    "sample": get_substring(
                        string=fb, pattern=self.config.filename_patterns.sample
                    ),
                    "modality": get_modality_from_string(
                        string=fb,
                        pattern=self.config.filename_patterns.modality,
                        accepted_modalities=self.accepted_modalities(),
                    ),
                    "technical_rep": get_substring(
                        string=fb, pattern=self.config.filename_patterns.technical_rep
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
        df["sample_rep"] = df["sample"] + "_" + df["technical_rep"]
        df = df.sort_values(["read_filename", "modality", "read_number"])
        self.table = df

    def subset(self, f):
        sub = copy.deepcopy(self)
        sub.table = sub.table[f(sub.table)]
        sub.filelist = [f for f in self.filelist if f in sub.table["read_path"].values]
        return sub

# }}}
