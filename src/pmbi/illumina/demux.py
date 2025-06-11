import argparse
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


def validate_sample_sheet(sample_sheet_path, logger=None):
    """
    Validate that a sample sheet has the required format and data.
    
    Args:
        sample_sheet_path (str or Path): Path to the sample sheet
        logger (logging.Logger, optional): Logger to use for logging messages
        
    Returns:
        bool: True if the sample sheet is valid, False otherwise
    """
    if logger is None:
        logger = logging.getLogger("validate_samplesheet")
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(name)s::%(levelname)s --> %(message)s"))
        logger.addHandler(handler)
    
    # Read samples from the sample sheet
    samples, success = read_samples_from_samplesheet(sample_sheet_path, logger)
    if not success:
        return False
    
    # Check if we have at least one sample
    if not samples:
        logger.error("No samples found in sample sheet")
        return False
    
    # Check for duplicate sample IDs
    sample_ids = [sample["sample_id"] for sample in samples]
    if len(sample_ids) != len(set(sample_ids)):
        logger.error("Duplicate sample IDs found in sample sheet")
        return False
    
    # Check for empty index fields if any sample has an index
    has_index = any(sample["index"] for sample in samples)
    if has_index:
        for sample in samples:
            if not sample["index"]:
                logger.warning(f"Sample {sample['sample_id']} is missing index information")
    
    logger.info(f"Sample sheet validation passed with {len(samples)} samples")
    return True

def read_samples_from_samplesheet(sample_sheet_path, logger=None):
    """
    Read sample information from an Illumina sample sheet.
    
    Args:
        sample_sheet_path (str or Path): Path to the sample sheet
        logger (logging.Logger, optional): Logger to use for logging messages
    
    Returns:
        list: List of dictionaries containing sample information (sample_id, index, index2)
        bool: Success status
    """
    if logger is None:
        logger = logging.getLogger("samplesheet_reader")
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(name)s::%(levelname)s --> %(message)s"))
        logger.addHandler(handler)
    
    samples = []
    data_section = False
    header = None
    
    try:
        with open(sample_sheet_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "[Data]":
                    data_section = True
                    continue
                
                if data_section and line and not line.startswith("["):
                    if "Sample_ID" in line:  # Header line
                        header = [col.strip() for col in line.split(',')]
                        continue
                    
                    if header:
                        parts = [part.strip() for part in line.split(',')]
                        if len(parts) >= 1 and parts[0]:
                            sample_info = {"sample_id": parts[0]}
                            
                            # Get index information if available
                            if "index" in header and len(parts) > header.index("index"):
                                sample_info["index"] = parts[header.index("index")]
                            else:
                                sample_info["index"] = ""
                                
                            if "index2" in header and len(parts) > header.index("index2"):
                                sample_info["index2"] = parts[header.index("index2")]
                            else:
                                sample_info["index2"] = ""
                                
                            samples.append(sample_info)
                    else:
                        # If we're in the data section but haven't found a header yet,
                        # and this line isn't a header, it might be malformed
                        logger.warning(f"Found data line before header: {line}")
    except Exception as e:
        logger.error(f"Error parsing sample sheet: {e}")
        return [], False
    
    if not samples:
        logger.error("No samples found in sample sheet")
        return [], False
        
    logger.info(f"Found {len(samples)} samples in sample sheet")
    return samples, True

def check_for_missing_fastq(output_dir, sample_sheet_path, logger=None):
    """
    Check if all expected fastq files based on the sample sheet are present in the output directory.
    
    Args:
        output_dir (Path): Path to the output directory containing fastq files
        sample_sheet_path (str or Path): Path to the sample sheet
        logger (logging.Logger, optional): Logger to use for logging messages
    
    Returns:
        bool: True if all expected files are present, False otherwise
        list: List of missing files if any
    """
    if logger is None:
        logger = logging.getLogger("check_fastq")
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(name)s::%(levelname)s --> %(message)s"))
        logger.addHandler(handler)
    
    # Parse sample sheet to get expected samples
    samples, success = read_samples_from_samplesheet(sample_sheet_path, logger)
    if not success:
        return False, []
    
    expected_samples = [sample["sample_id"] for sample in samples]
    
    # Look for fastq files in output directory
    # Typical pattern: SampleName_S1_L001_R1_001.fastq.gz
    missing_files = []
    
    # Check if the output directory exists
    if not output_dir.exists():
        logger.error(f"Output directory does not exist: {output_dir}")
        return False, [f"Output directory not found: {output_dir}"]
    
    # Get all fastq files in the output directory and its subdirectories
    all_fastq_files = list(output_dir.glob("**/*.fastq.gz"))
    
    # Check if each expected sample has corresponding fastq files
    for sample_id in expected_samples:
        # Look for R1 and R2 files for each sample
        r1_pattern = f"{sample_id}_S*_*_R1_*.fastq.gz"
        r2_pattern = f"{sample_id}_S*_*_R2_*.fastq.gz"
        
        r1_files = list(output_dir.glob(f"**/{r1_pattern}"))
        r2_files = list(output_dir.glob(f"**/{r2_pattern}"))
        
        if not r1_files:
            missing_files.append(f"{sample_id} (R1)")
        
        if not r2_files:
            missing_files.append(f"{sample_id} (R2)")
    
    if missing_files:
        logger.warning(f"Missing fastq files for: {', '.join(missing_files)}")
        return False, missing_files
    
    logger.info(f"All expected fastq files found for {len(expected_samples)} samples")
    return True, []

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
        stdout=subprocess.PIPE, # Affirms that proc.stdout will not be None
        stderr=subprocess.STDOUT,
        text=True
    ) as proc:
        
        logger = logging.getLogger("bcl2fastq")
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        # handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        handler.setFormatter(logging.Formatter("%(name)s::%(levelname)s --> %(message)s"))
        logger.addHandler(handler)
        while True:
            line = proc.stdout.readline()
            if not line and proc.poll() is not None:
                break
                
            if line:
                logger.info(line.strip())

        if proc.returncode != 0:
            mes = f"Demux process failed with returncode: {proc.returncode}"
            logger.critical(mes)
            raise subprocess.SubprocessError(mes)

    # Verify that the output fastq files match what is expected 
    # according to the sample sheet
    output_dir = fastq_dir.joinpath(run_id)
    success, missing_files = check_for_missing_fastq(output_dir, sample_sheet, logger)
    
    if not success:
        logger.warning(f"Some expected fastq files are missing: {missing_files}")
        logger.warning("Demultiplexing completed with warnings.")
    else:
        logger.info("Demultiplexing completed successfully. All expected files found.")

