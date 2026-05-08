import importlib
from pathlib import Path
from typing import Any, Callable

import pandas as pd
import toolz
from munch import Munch

import pmbi.config._config as pmbiconf
import pmbi.config.item as item
from pmbi.collections import FileCollection
from pmbi.logging import streamLogger


class ProcessUnit:
    def __init__(
        self,
        table: pd.DataFrame,
        sample_name_column: str = "sample_name",
        file_path_column: str = "path",
    ):
        self._clsname = self.__class__.__name__
        self.table = table
        self.files = self.table[file_path_column].to_list()
        self.sample_name_column = sample_name_column
        self.file_path_column = file_path_column

        self.sample = self._validate_single_sample()

    def _validate_single_sample(self) -> str:
        """Ensure unit contains exactly one sample"""
        unique_samples = self.table[self.sample_name_column].unique()
        if len(unique_samples) != 1:
            raise ValueError(
                f"ProcessUnit must contain exactly one sample. Found: {unique_samples}"
            )
        return unique_samples[0]


class KallistoUnit(ProcessUnit):
    def __init__(
        self,
        table: pd.DataFrame,
        config: Munch,
        sample_name_column: str = "sample_name",
        file_path_column: str = "path",
        read_number_column: str = "read_number",
    ):
        super().__init__(table, sample_name_column, file_path_column)
        self.config: Munch = (
            pmbiconf.Config.from_items(
                [
                    item.PathItem("kallisto.reference", must_exist=True),
                    item.PathItem("kallisto.outdir", must_exist=True),
                    item.Item("kallisto.mode", str),
                    item.Item("kallisto.nproc", int, optional=True, default=1),
                ]
            ).set_values_from_config(config, inplace=False)
        ).to_munch()

        # Make sure read numbers look right
        assert all(
            [
                rn in [f"R{n}" for n in range(1, 4)]
                for rn in self.table[read_number_column].to_list()
            ]
        )

        read_indexed_table = self.table.set_index(read_number_column, drop=False)
        self.read1 = read_indexed_table.loc["R1", file_path_column]
        self.read2 = read_indexed_table.loc["R2", file_path_column]

class ProcessRunner:
    def __init__(self, unit: ProcessUnit, wd: Path = Path("."), **kwargs):
        self.logger = streamLogger(self.__class__.__name__)
        self.unit = unit
        self.wd = wd
        self.kwargs = kwargs
        self.binary = None

    def _stringify(f: Callable) -> Callable:
        def _inner(self):
            return list(map(str, f(self)))

        return _inner

    @_stringify
    def _cmd(self) -> list[str]:
        return []

    def _setup(self) -> None:
        pass

    def _cleanup(self) -> None:
        pass

    def more_args(self) -> list[str]:
        """Convert kwargs to command line arguments"""
        return list(toolz.concat([[f"--{k}", v] for k, v in self.kwargs.items()]))

    def run(self) -> None:
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

class KallistoRunner(ProcessRunner):
    def __init__(self, unit: KallistoUnit, wd: Path = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.binary = "kallisto"

    @ProcessRunner._stringify
    def _cmd(self) -> list[str]:
        return [
                self.binary,
                "quant",
                "-i", self.unit.config.kallisto.reference,
                "-o", self.unit.config.kallisto.outdir,
                "-t", self.unit.config.kallisto.nproc,
                self.unit.read1,
                self.unit.read2,
                ]

    def _setup(self) -> None:
        super()._setup()

    def _cleanup(self) -> None:
        super()._cleanup()

    def run(self):
        super().run()
