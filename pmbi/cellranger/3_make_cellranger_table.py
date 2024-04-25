# %% imports
import boto3
import argparse
import re
import os
from pmbi.s3.lib import split_s3_uri, object_key_list

# %%
object_key_list("s3://upenn-research.thaiss-lab-psom-01.us-east-1/snATAC/demux/fastq")

# %%
def get_readnum_from_filename(pattern: str, objkey: str) -> str:
    s = re.search(pattern, os.path.basename(objkey))
    if s is not None:
        if len(s.groups()) > 1:
            raise ValueError(f"More than one matching read number: {objkey}")
        else:
            return s.group(1)
    else:
        raise ValueError(f"No matching read number: {objkey}")

# This function defaults to assuming that the protocol type is present in the second underscore-separated
# column of the file name, which should be the case if using `make_sample_sheet.bash` 
def get_modality_from_filename(fname: os.PathLike, sep: str = "_", pos: int = 1, protocol_converter: dict[str,str]|None = None) -> str:
    protocol = fname.split(sep)[pos]
    if protocol_converter is not None:
        assert isinstance(protocol_converter, dict), "protocol_converter passed to object_protocol_type needs to be type dict"
        try:
            return protocol_converter[protocol]
        except KeyError:
            raise
    else:
        return protocol

emtab_style_protocol_converter = {
        "RNA": "scRNA-seq",
        "VDJ": "scVDJ-seq",
        "ADT": "scADT-seq",
        "snRNA-seq": "scADT-seq",
        "HTO": "HTO",
        }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--s3uri", required = True, type = str, help = "S3 URI where reads are located.")
    parser.add_argument("-i", "--include", required = False, default = "[.]fastq[.]gz$", type = str, help = "Python-style regex that matches some of all the filenames to include.")
    parser.add_argument("-r", "--readnum_pattern", required = False, default = "[_](R[0-9])[_]" , type = str, help = "Python-style regex that matches the read number.")
    args = parser.parse_args()

    print('\t'.join(["modality", "read_number", "read_filename", "read_s3_uri"]))
    for objkey in object_key_list(args.s3uri):
        if re.search(args.include, objkey) is not None:
             _dirname, basename = os.path.split(objkey)
             object_uri = os.path.join(args.s3uri, basename)
             modality = get_modality_from_filename(basename)
             readnum = get_readnum_from_filename(pattern = args.readnum_pattern, objkey = basename)
             print('\t'.join([modality, readnum, basename, object_uri]))
