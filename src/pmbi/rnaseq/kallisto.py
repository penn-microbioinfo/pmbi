import importlib
import logging
import os
import shutil
from pathlib import Path
from typing import Any, Callable, Union

import pandas as pd
import toolz
from munch import Munch

import pmbi.config._config as pmbiconf
import pmbi.config.item as item
from pmbi.collections import FileCollection
from pmbi.logging import stream_file_logger, streamLogger
from pmbi.subproc.process_runner import ProcessRunner
from pmbi.subproc.process_unit import ProcessUnit


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
                    item.Item("kallisto.mode", str),
                    item.Item("kallisto.nproc", int, optional=True, default=1),
                    item.Item("kallisto.force", bool, optional=True, default=False),
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


class KallistoRunner(ProcessRunner):
    def __init__(self, unit: KallistoUnit, wd: Union[Path, str] = Path("."), **kwargs):
        super().__init__(unit, wd, **kwargs)
        self.binary = "kallisto"
        self.cmd_logger = stream_file_logger(
            self.clsname, self.wd.joinpath("kallisto.log")
        )

    @ProcessRunner._stringify
    def _cmd(self) -> list[str]:
        return [
            self.binary,
            "quant",
            "-i",
            self.unit.config.kallisto.reference,
            "-o",
            self.sample_wd(),
            "-t",
            self.unit.config.kallisto.nproc,
            self.unit.read1,
            self.unit.read2,
        ] + self.more_args()

    def _setup(self) -> None:
        super()._setup()
        if self.sample_wd().exists():
            if self.unit.config.kallisto.force:
                self.logger.warning(
                    f"Removing existing output directory because force==True: {self.sample_wd()}"
                )
                shutil.rmtree(self.sample_wd())
            else:
                raise OSError(
                    f"Sample output directory exists and force==False. Skipping sample: {self.unit.sample_id}"
                )

        os.mkdir(self.sample_wd())

    def _cleanup(self) -> None:
        super()._cleanup()

    def run(
        self,
    ):
        super().run()
