import argparse
import json
import subprocess
import time
import typing
import xml.etree.ElementTree as ET
from io import StringIO, TextIOWrapper
from pathlib import Path

import pandas as pd

import pmbi.error as perr
from pmbi.logging import list_loggers, streamLogger


def check_avail(binary, help="-h"):
    try:
        subprocess.check_call(
            [binary, help], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        return True
    except FileNotFoundError:
        return False


def get_run_id(rundir: Path, tag: str = "BaseSpaceRunId"):
    params_xml = rundir.joinpath("RunParameters.xml")
    if not params_xml.exists():
        raise FileNotFoundError
    tree = ET.parse(params_xml)
    run_id = tree.findtext(tag)
    if run_id is None:
        raise ValueError(f"Unable to find {tag} in RunParameters.xml at: {params_xml}")
    return run_id


# %%
def read_samplesheet_data_section(sample_sheet: typing.TextIO) -> pd.DataFrame:
    """
    Read sample information from an Illumina sample sheet.

    Args:
        sample_sheet_path (TextIO): Path to the sample sheet

    Returns:
        pd.DataFrame: DataFrame from the comma-separated data in the [Data] section of the sample sheet. Raise ValueError if the [Data] section is missing.
    """
    data_section = None

    for line in sample_sheet:
        line = line.strip()
        if line == "[Data]":
            data_section = StringIO()
            continue

        if data_section and line and not line.startswith("["):
            data_section.write(f"{line}\n")

    if data_section:
        data_section.seek(0)
        return pd.read_csv(data_section)
    else:
        raise ValueError("Data section missing")


def fastq_missing_for(output_dir: Path, sample_sheet: typing.TextIO):
    """
    Check if all expected fastq files based on the sample sheet are present in the output directory.

    Args:
        output_dir (Path): Path to the output directory containing fastq files
        sample_sheet (TextIO): File-like object of smaple sheet

    Returns:
        None
    """

    logger = streamLogger("check_for_missing_fastq")

    # Parse sample sheet to get expected samples
    expected_samples = read_samplesheet_data_section(sample_sheet)["Sample_ID"]

    # Check if the output directory exists
    if not output_dir.exists():
        mes = f"Output directory does not exist: {output_dir}"
        logger.error(mes)
        raise OSError(mes)

    # Look for sample-matching fastq files in output directory
    missing_samples = []
    for sample_id in expected_samples:

        r1_files = list(output_dir.glob(f"{sample_id}_S*_*_R1_*.fastq.gz"))
        r2_files = list(output_dir.glob(f"{sample_id}_S*_*_R2_*.fastq.gz"))

        if len(r1_files) == 0 or len(r2_files) == 0:
            missing_samples.append(sample_id)

    if len(missing_samples) > 0:
        logger.warning(f"Missing fastq files for: {', '.join(missing_samples)}")

    else:
        logger.info(
            f"All expected fastq files found for {len(expected_samples)} samples"
        )

    return missing_samples


def bcl_convert_cmd():
    raise NotImplementedError


def bcl2fastq_cmd(
    runfolder_dir,
    output_dir,
    sample_sheet,
    loading_threads=4,
    processing_threads=4,
    writing_threads=4,
    create_fastq_for_index_reads=False,
):
    cmd = [
        "bcl2fastq",
        "-R",
        runfolder_dir,
        "--sample-sheet",
        sample_sheet,
        "-o",
        output_dir,
        "-r",
        loading_threads,
        "-p",
        processing_threads,
        "-w",
        writing_threads,
    ]
    if create_fastq_for_index_reads:
        cmd.append("--create-fastq-for-index-reads")

    return cmd


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("rundir", help="Path to Illumina run directory.")
    parser.add_argument(
        "-s",
        "--sample_sheet",
        required=False,
        help="Path to sample sheet to use for demuxing. Must conform to Illumina SampleSheet.csv format.",
    )
    parser.add_argument(
        "-o",
        "--fastq_dir",
        default="./",
        required=False,
        help="Directory within which to place demux output directory/files. Directory will be automatically named by BaseSpaceRunId, by default.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=4,
        required=False,
        help="Number of threads to apply to all thread options equally.",
    )
    parser.add_argument(
        "-b",
        "--backend",
        required=False,
        default="bcl2fastq",
        choices=["bcl2fastq", "bcl-convert"],
        help="Backend to use for demuxing: bcl2fastq or bcl-convert",
    )
    parser.add_argument(
        "--create_fastq_for_index_reads",
        action="store_true",
        required=False,
        help="Pass this flag to create index read fastq files",
    )

    args = parser.parse_args()

    fastq_dir = Path(args.fastq_dir)

    # Default Illumina SampleSheet.csv path
    if args.sample_sheet is None:
        sample_sheet = fastq_dir.joinpath("SampleSheet.csv")
    else:
        sample_sheet = args.sample_sheet

    # Check for backend availability
    if not check_avail(args.backend):
        raise OSError(f"Backend not available: {args.backend}")

    run_id = get_run_id(Path(args.rundir))

    if args.backend == "bcl2fastq":
        cmd = bcl2fastq_cmd(
            runfolder_dir=args.rundir,
            output_dir=fastq_dir.joinpath(run_id),
            sample_sheet=sample_sheet,
            loading_threads=args.threads,
            processing_threads=args.threads,
            writing_threads=args.threads,
            create_fastq_for_index_reads=args.create_fastq_for_index_reads,
        )
    elif args.backend == "bcl-convert":
        raise NotImplementedError

    else:
        # Should not happen because of `choices` in argparser
        raise ValueError

    # Run the demux command
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,  # Affirms that proc.stdout will not be None
        stderr=subprocess.STDOUT,
        text=True,
    ) as proc:
        logger = streamLogger("bcl2fastq")
        while True:
            line = proc.stdout.readline()
            if not line and proc.poll() is not None:
                break

            if line:
                logger.info(line.strip())

        if proc.returncode != 0:
            mes = f"Process failed with returncode: {proc.returncode}"
            logger.critical(mes)
            raise subprocess.SubprocessError(mes)

    # %% TODO: Swap for `fastq_missing_for`
    # This function only writes to log
   # check_for_missing_fastq(
    #     output_dir=fastq_dir.joinpath(run_id), sample_sheet=open(sample_sheet, "r")
    # )
