# %%
import sys
import argparse
from pathlib import Path
import pandas as pd

from joblib import Parallel, delayed

import pmbi.cellranger.cellranger_command as crc
from pmbi.config import import_config
from pmbi.file_handlers import LocalBackend
from pmbi.logging import streamLogger

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", help = "Path to FASTQ files to run `cellranger multi` on")
    parser.add_argument("-c", "--config", help = "Path to PMBI cellranger configuration file.")
    parser.add_argument("-o", "--output", help = "Path to generate cellranger outputs.")
    parser.add_argument("-n", "--dry", action="store_true", help = "Do a dry run")
    args = parser.parse_args()

    fastq_dir = args.fastq
    config_path = args.config
    output_dir = args.output

# %% Create Collection object, etc
    logger = streamLogger("process_cellranger")
    config = import_config( config_path )

    h = crc.CellrangerCollection(
        path=Path(fastq_dir),
        config=config,
        pattern=config.filename_patterns.include,
        backend=LocalBackend(),
    )

# %% Log some validating stats 
    logger.info(f"Unique samples: {h.samples().shape[0]}")
    logger.info(f"CellrangerCollection table dimensions: {h.table.shape}")
    n_modals = h.modalities_per_sample()
    if n_modals.value_counts().shape[0] == 1:
        logger.info(f"All {h.samples().shape[0]} samples are represented by {n_modals.value_counts().index[0]} modalities: ")
    else:
        logger.warning(f"Some of the {h.samples().shape[0]} samples are represented by different numbers of modalities: {", ".join(pd.Series(n_modals.value_counts().index).astype(str).tolist())} ")
    logger.info(f"Modalities represented: {", ".join(h.table["modality"].unique().tolist())}")

# %%

# %%
    if not args.dry:
        def _run(unit, wd, **kwargs):
            unit.create_runner(wd, **kwargs).run()

        _out = Parallel(n_jobs=args.run.n_jobs)(delayed(_run)(unit,
                                                wd=Path(output_dir),
                                                localcores=config.commands.cores,
                                                localmem=config.commands.memory
                                                ) for unit in h.get_units())

    else:
        logger.info("Dry run, so doing nothing and exiting")

