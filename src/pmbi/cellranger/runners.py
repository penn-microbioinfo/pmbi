import os
import re
import shutil
from pathlib import Path
from subprocess import SubprocessError
from typing import TYPE_CHECKING

import toolz

import pmbi.subproc as pmbiproc
from pmbi.cellranger.units import (
    CellrangerArcUnit,
    CellrangerMultiUnit,
    CellrangerUnit,
)
import pmbi.illumina.fastq as pfq
from pmbi.logging import streamLogger

# if TYPE_CHECKING:
#     from pmbi.cellranger.cellranger_command import CellrangerUnit, CellrangerArcUnit, CellrangerMultiUnit


class CellrangerRunner:

    def __init__(self, unit: CellrangerUnit, wd: Path = Path("."), **kwargs):
        self.unit = unit
        self.wd = Path(wd)
        self.id = self.unit.sample
        self.outs_dest = self.wd / f"{self.id}__outs"
        self.logger = streamLogger(self.__class__.__name__)
        self.kwargs = kwargs

    def more_args(self) -> list[str]:
        """Convert kwargs to command line arguments"""
        return list(toolz.concat([[f"--{k}", v] for k, v in self.kwargs.items()]))

    def cmd(self) -> list[str]:
        return []

    def setup(self) -> None:
        pass

    def cleanup(self) -> None:
        pass

    def run(self) -> None:
        pass


class CellrangerAtacRunner(CellrangerRunner):
    def __init__(self, unit: CellrangerArcUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.logger = streamLogger("CellrangerAtacRunner")

    def _setup(self):
        pass

    def _cmd(self) -> list[str]:
        return list(
            map(
                str,
                (
                    [
                        "cellranger-atac",
                        "count",
                        "--id",
                        self.id,
                        "--reference",
                        self.unit.reference,
                        "--fastqs",
                        self.unit.fastqs,
                        "--sample",
                        self.unit.sample
                    ]
                    + self.more_args()
                ),
            )
        )
        
    def run(self):
        if self.outs_dest.exists():
            mes = f"Output directory exists. Skipping myself: {self.id}"
            self.logger.critical(mes)
            raise OSError(mes)
        else:
            cmd_logger = streamLogger("cellranger-atac")
            self._setup()
            self.logger.info(f"Changing directory: {self.wd}")
            os.chdir(self.wd)
            pmbiproc.run_and_log(cmd=self._cmd(), logger=cmd_logger)
            # self.cleanup()

class CellrangerArcRunner(CellrangerRunner):
    def __init__(self, unit: CellrangerArcUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.logger = streamLogger("CellrangerArcRunner")
        self.csv_suffix = "__libraries.csv"
        self.csv_path = self.wd / f"{self.id}{self.csv_suffix}"

    def _setup(self):
        self.wd = Path(self.wd)
        self.logger.info(f"Writing configuration CSV to: {self.csv_path}")
        with open(self.csv_path, "w") as csv:
            csv.write(self.unit.config_csv().read())

    def _cmd(self) -> list[str]:
        more_args = list(toolz.concat([[f"--{k}", v] for k, v in self.kwargs.items()]))
        return list(
            map(
                str,
                (
                    [
                        "cellranger-arc",
                        "count",
                        "--id",
                        self.id,
                        "--reference",
                        self.unit.reference,
                        "--libraries",
                        self.csv_path,
                        "--create-bam",
                        "true",
                    ]
                    + more_args
                ),
            )
        )

    def run(self):
        if self.outs_dest.exists():
            mes = f"Output directory exists. Skipping myself: {self.id}"
            self.logger.critical(mes)
            raise OSError(mes)
        else:
            cmd_logger = streamLogger("cellranger-arc")
            self._setup()
            self.logger.info(f"Changing directory: {self.wd}")
            os.chdir(self.wd)
            pmbiproc.run_and_log(cmd=self._cmd(), logger=cmd_logger)
            # self.cleanup()

class CellrangerMultiRunner(CellrangerRunner):
    def __init__(self, unit: CellrangerMultiUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.logger = streamLogger("CellrangerMultiRunner")
        self.csv_suffix = "__multi_config.csv"
        self.csv_path = self.wd / f"{self.id}{self.csv_suffix}"

    def _setup(self):
        self.wd = Path(self.wd)
        self.logger.info(f"Writing configuration CSV to: {self.csv_path}")
        with open(self.csv_path, "w") as csv:
            csv.write(self.unit.config_csv().read())

    def _cmd(self) -> list[str]:
        cmd = [
            "cellranger",
            "multi",
            "--id",
            self.id,
            "--csv",
            self.csv_path,
        ] + self.more_args()
        return list(map(str, cmd))

    def _cleanup(self):
        # Delete temporary files
        for tf in self.unit._tempfile_paths:
            try:
                os.remove(tf)
            except OSError:
                self.logger.warning(
                    f"Temporary file marked for deletion does not exist: {tf}"
                )

    def run(self):
        if self.outs_dest.exists():
            mes = f"Output directory exists. Skipping myself: {self.id}"
            self.logger.critical(mes)
            raise OSError(mes)
        else:
            cmd_logger = streamLogger("cellranger-multi")
            self._setup()
            self.logger.info(f"Changing directory: {self.wd}")
            os.chdir(self.wd)
            try:
                pmbiproc.run_and_log(cmd=self._cmd(), logger=cmd_logger)
                self._cleanup()
            except SubprocessError:
                self._cleanup()


            # self.cleanup()


