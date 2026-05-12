import os
from pathlib import Path
from subprocess import SubprocessError
from typing import Callable, Union
import logging

import toolz

import pmbi.subproc as pmbiproc
from pmbi.logging import streamLogger
from pmbi.subproc.process_unit import ProcessUnit


class ProcessRunner:
    def __init__(self, unit: ProcessUnit, wd: Union[Path,str] = Path("."), **kwargs):
        self.clsname = self.__class__.__name__
        self.logger = streamLogger(self.clsname)
        self.unit = unit
        self.wd = Path(wd)
        self.kwargs = kwargs
        self.binary = None
        self.cmd_logger = None

    @staticmethod
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

    def sample_wd(self):
        return self.wd.joinpath(self.unit.sample_id)

    def run(self) -> None:
        if self.cmd_logger is None:
            self.cmd_logger = streamLogger(self.binary)
        self._setup()
        self.logger.info(f"Changing directory: {self.wd}")
        os.chdir(self.wd)
        try:
            pmbiproc.run_and_log(cmd=self._cmd(), logger=self.cmd_logger)
            self._cleanup()
        except SubprocessError:
            self.logger.critical("command run failed with SubprocessError")
            self._cleanup()
