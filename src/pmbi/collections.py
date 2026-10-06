import os
from pathlib import Path
from typing import Union
from itertools import chain
import natsort

# from munch import Munch
import pandas as pd

from pmbi.file_handlers import Backend, LocalBackend
from pmbi.util.misc import get_substring, substring_detected
from pmbi.logging import streamLogger

# type PatternDict = dict[str,str]

logger = streamLogger(__name__)

# %%
class FileCollection:
    def __init__(self):
        pass

class ReadFileCollection(FileCollection):
    def __init__(
        self,
        path: Path | list[Path],
        patterns: dict[str, str],
        backend: Backend = LocalBackend(),
        sort_by: Union[list[str], None] = None,
        sample_name_column = "sample_name"
    ):

        self.path = path
        self.patterns = patterns
        self.backend = backend
        self.sample_name_column = sample_name_column

        if self.patterns["include"] is None and self.patterns["exclude"] is None:
            logger.warning("FYI, `include` and `exclude` missing from patterns")

        if self.patterns["exclude"] is not None and self.patterns["include"] is None:
            raise ValueError("Cannot configure an `exclude` pattern without an `include` pattern")

        if sort_by is None:
            self.sort_by = []
        else:
            self.sort_by = sort_by

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
            add_file = self.patterns["include"] is None
            fb = os.path.basename(f)
            # if "include" in self.patterns:
            if self.patterns["include"] is not None:
                add_file = substring_detected(self.patterns["include"], fb)
                # if "exclude" in self.patterns and add_file:
                if self.patterns["exclude"] is not None and add_file:
                    add_file = not substring_detected(self.patterns["exclude"], fb)
            if add_file:
                filtered.append(f)

        return filtered

    def _make_table(self) -> pd.DataFrame:
        """Create a table of file metadata"""
        assert not any([x in self.patterns for x in ["filename", "path"]])
        ldict = []
        for f in self.filelist:
            fb = f.name
            row = {}
            for k, v in self.patterns.items():
                if k in ["include", "exclude"]:
                    continue

                row[k] = get_substring(string=fb, pattern=v, group=1)

            row["filename"] = fb
            row["path"] = f
            ldict.append(row)

        df = pd.DataFrame(ldict)

        df = df.sort_values(
            self.sort_by, key=natsort.natsort_keygen()
        ).reset_index(drop=True)

        return df

    def split_by(self, column) -> dict[str, pd.DataFrame]:
        return {val: grp for val, grp in self.table.groupby(column)}

    def samples(self):
        return self.table[self.sample_name_column].unique()

    def get_units(self, unit_type, **kwargs) -> list:
        units = []
        for samp,table in self.split_by(self.sample_name_column).items():
            units.append(
                    unit_type(
                        table=table,
                        **kwargs
                        )
                    )
        return units
# %%

