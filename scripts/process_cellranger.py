# %%
import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed
from munch import Munch

from pmbi.cellranger.collections import CellrangerCollection
from pmbi.cellranger.runners import (
    CellrangerArcRunner,
    CellrangerAtacRunner,
    CellrangerMultiRunner,
)
from pmbi.config import import_config, table_to_dataframe
from pmbi.file_handlers import LocalBackend
from pmbi.logging import streamLogger


def log_misc_stats(collection):
    logger = logging.getLogger("process_cellranger")
    logger.info(f"Unique samples: {collection.samples().shape[0]}")
    logger.info(f"CellrangerCollection table dimensions: {collection.table.shape}")
    n_modals = collection.modalities_per_sample()
    if n_modals.value_counts().shape[0] == 1:
        logger.info(
            f"All {collection.samples().shape[0]} samples are represented by {n_modals.value_counts().index[0]} modalities: "
        )
    else:
        logger.warning(
            f"Some of the {collection.samples().shape[0]} samples are represented by different numbers of modalities: {", ".join(pd.Series(n_modals.value_counts().index).astype(str).tolist())} "
        )
    logger.info(
        f"Modalities represented: {", ".join(collection.table["modality"].unique().tolist())}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", help="Path to PMBI cellranger configuration file."
    )
    parser.add_argument("-n", "--dry", action="store_true", help="Do a dry run")
    parser.add_argument(
        "-j",
        "--njobs",
        action="store",
        required=False,
        help="Number of parallel jobs to run; will overwrite config.run.njobs",
    )
    args = parser.parse_args()

    # %% Create Collection object, etc
    logger = streamLogger("process_cellranger")
    config = import_config(args.config)
    if args.njobs is not None:
        config.run.njobs = args.njobs

    fastq_paths = (
        table_to_dataframe(config.modalities)["fastqs"].dropna().apply(Path).tolist()
    )
    if len(fastq_paths) == 0:
        raise ValueError("Unable to identify any fastq file paths in config")

    logger.info(
        f"Pulling fastq files from {len(fastq_paths)} path:\n {'\n'.join([os.fspath(f) for f in fastq_paths])}"
    )

    collection = CellrangerCollection(
        path=fastq_paths, config=config, backend=LocalBackend()
    )

    units = collection.get_units()

    runners = [
        collection.RunnerType(unit=u, wd=config.run.wd, **config.command_line.__dict__)
        for u in units
    ]

    if args.dry:
        log_misc_stats(collection)

        logger.info("-" * 25)
        logger.info("Commands that will be run:\n")
        logger.info("\n".join([" ".join(r._cmd()) for r in runners]))
        logger.info("-" * 25)

    else:
        _out = Parallel(n_jobs=config.run.njobs)(delayed(r.run)() for r in runners)
