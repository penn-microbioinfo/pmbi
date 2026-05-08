import os
import re
import shutil
from contextlib import nullcontext
from pathlib import Path
from subprocess import SubprocessError
from typing import TYPE_CHECKING

import toolz

import pmbi.illumina.fastq as pfq
import pmbi.subproc as pmbiproc
from pmbi.cellranger.units import CellrangerArcUnit, CellrangerMultiUnit, CellrangerUnit, StarSoloUnit
from pmbi.logging import streamLogger
from pmbi.sync import PathSync

# if TYPE_CHECKING:
#     from pmbi.cellranger.cellranger_command import CellrangerUnit, CellrangerArcUnit, CellrangerMultiUnit


class CellrangerRunner:

    def __init__(self, unit: CellrangerUnit, wd: Path = Path("."), **kwargs):
        self.unit = unit
        self.wd = Path(wd)
        self.id = self.unit.sample
        self.binary = "generic-cellranger-runner"
        self.outs_src = "outs"
        self.outs_dest = self.wd / f"{self.id}__outs"
        self.logger = streamLogger(self.__class__.__name__)
        self.kwargs = kwargs

    def more_args(self) -> list[str]:
        """Convert kwargs to command line arguments"""
        return list(toolz.concat([[f"--{k}", v] for k, v in self.kwargs.items()]))

    def _cmd(self) -> list[str]:
        return []

    def _setup(self) -> None:
        pass

    def _cleanup(self) -> None:
        # Delete temporary files
        for tf in self.unit._tempfile_paths:
            try:
                self.logger.debug(f"Removing temporary file: {tf}")
                os.remove(tf)
            except OSError:
                self.logger.warning(
                    f"Temporary file marked for deletion does not exist: {tf}"
                )

        return None
        # Sync outs to final resting place
        rundir = self.wd.joinpath(self.id)
        outs_src = rundir / self.outs_src
        ps = PathSync(src_dir=outs_src, dest_dir=self.outs_dest, files_to_sync=["outs"])
        try:
            self.logger.info(f"Syncing `outs` directory to: {self.outs_dest}")
            ps.sync()

            self.logger.info(f"Removing cellranger intermediate run files under: {rundir}")
            shutil.rmtree(rundir)
            
        except (NotADirectoryError, FileNotFoundError, OSError) as e:
            self.logger.critical(
                "Errors occurred during syncing of outs dir..."
                f"{e}"
                "Continuing with any remaining cleanup despite errors."
            )

    def run(self) -> None:
        if self.wd.joinpath(self.id).exists():
            mes = f"Cellranger output directory already exists. Skipping myself: {self.id}"
            self.logger.critical(mes)
            raise OSError(mes)
        if self.outs_dest.exists():
            mes = f"Future cleanup `outs` directory already exists. Skipping myself: {self.id}"
            self.logger.critical(mes)
            raise OSError(mes)
        cmd_logger = streamLogger(self.binary)
        self._setup()
        self.logger.info(f"Changing directory: {self.wd}")
        os.chdir(self.wd)
        try:
            pmbiproc.run_and_log(cmd=self._cmd(), logger=cmd_logger)
            self._cleanup()
        except SubprocessError:
            self.logger.critical("command run failed with SubprocessError")
            self._cleanup()


class CellrangerAtacRunner(CellrangerRunner):
    def __init__(self, unit: CellrangerArcUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.logger = streamLogger("CellrangerAtacRunner")
        self.binary = "cellranger-atac"

    def _setup(self):
        pass

    def _cmd(self) -> list[str]:
        return list(
            map(
                str,
                (
                    [
                        self.binary,
                        "count",
                        "--id",
                        self.id,
                        "--reference",
                        self.unit.reference,
                        "--fastqs",
                        self.unit.fastqs,
                        "--sample",
                        self.unit.sample,
                    ]
                    + self.more_args()
                ),
            )
        )

    def cleanup(self):
        super()._cleanup()

    def run(self):
        super().run()


class CellrangerArcRunner(CellrangerRunner):
    def __init__(self, unit: CellrangerArcUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.logger = streamLogger("CellrangerArcRunner")
        self.binary = "cellranger-arc"
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
                        self.binary,
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
        super().run()


class CellrangerMultiRunner(CellrangerRunner):
    def __init__(self, unit: CellrangerMultiUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.logger = streamLogger("CellrangerMultiRunner")
        self.binary = "cellranger"
        self.csv_suffix = "__multi_config.csv"
        self.csv_path = self.wd / f"{self.id}{self.csv_suffix}"

    def _setup(self):
        self.wd = Path(self.wd)
        self.logger.info(f"Writing configuration CSV to: {self.csv_path}")
        with open(self.csv_path, "w") as csv:
            csv.write(self.unit.config_csv().read())

    def _cmd(self) -> list[str]:
        cmd = [
            self.binary,
            "multi",
            "--id",
            self.id,
            "--csv",
            self.csv_path,
        ] + self.more_args()
        return list(map(str, cmd))

    def _cleanup(self):
        super()._cleanup()

    def run(self):
        super().run()

        # self.cleanup()

class StarSoloRunner(CellrangerRunner):
    def __init__(self, unit: StarSoloUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.logger = streamLogger("StarSoloRunner")
        self.binary = "STAR"

    def _setup(self):
        pass

    def _cmd(self) -> list[str]:
        return list(
            map(
                str,
                (
                    [
                        self.binary,
                        "--runMode", "alignReads",
                        "--genomeDir", self.unit.reference,
                        "--outFileNamePrefix", self.id,
                        "--readFilesCommand", "zcat",
                        "--readFilesIn", self.unit.read2, self.unit.read1,
                        "--soloType", "CB_UMI_Simple",
                        "--soloCBwhitelist", self.unit.config.run.barcode_whitelist,
                        "--soloUMIlen", 12,
                        "--soloCBlen", 16,
                        "--clipAdapterType", "CellRanger4",
                        "--outFilterScoreMin", 30,
                        "--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts",
                        "--soloUMIfiltering", "MultiGeneUMI_CR",
                        "--soloUMIdedup", "1MM_CR",
                        "--soloCellFilter", "EmptyDrops_CR",
                        "--soloFeatures", "Gene", "Velocyto", 
                        "--outSAMattributes", "CB", "UB", "sS", "sM",
                        "--outSAMtype", "BAM", "SortedByCoordinate",
                    ]
                        + self.more_args()
                ),
            )
        )

    def cleanup(self):
        pass
        # super()._cleanup()

    def run(self):
        super().run()
