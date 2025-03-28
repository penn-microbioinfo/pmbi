import subprocess
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import argparse


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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("rundir", help = "Path to Illumina run directory.")
    args = parser.parse_args()

    print("bcl2fastq", check_avail("bcl2fastq"))
    print("bcl-convert", check_avail("bcl-convert"))

    print(get_run_id(Path(args.rundir)))
