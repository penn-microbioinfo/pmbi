import argparse
import sys
import time
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
import logging


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


def validate_sample_sheet():
    raise NotImplementedError


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
        cmd.append("--create_fastq_for_index_reads")

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

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
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
        )
        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        ) as proc:
            
            logger = logging.getLogger("bcl2fastq")
            logger.setLevel(logging.INFO)
            # logger.handlers.clear()
            handler = logging.StreamHandler(stream=sys.stdout)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
            handler.setFormatter(logging.Formatter("%(message)s"))
            logger.addHandler(handler)
            while True:
                line = proc.stdout.readline()
                if not line and proc.poll() is not None:
                    break
                    
                if line:
                    logger.info(line.strip())

            print(proc.returncode)
                

    elif args.backend == "bcl-convert":
        raise NotImplementedError

    else:
        # Should not happen because of `choices` in argparser
        raise ValueError
