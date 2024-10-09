import argparse
import io
import logging
import gzip
import tarfile
import re
import os
import glob
logging.basicConfig(level=logging.INFO)

'''
important_outs/
important_outs/raw_molecule_info.h5
important_outs/web_summary.html
important_outs/vdj_t/
important_outs/vdj_t/filtered_contig_annotations.csv
important_outs/vdj_t/vloupe.vloupe
important_outs/vdj_t/airr_rearrangement.tsv
important_outs/vdj_t/concat_ref.fasta.fai
important_outs/vdj_t/concat_ref.bam
important_outs/vdj_t/filtered_contig.fasta
important_outs/vdj_t/consensus.fasta.fai
important_outs/vdj_t/consensus.bam.bai
important_outs/vdj_t/vdj_contig_info.pb
important_outs/vdj_t/consensus.bam
important_outs/vdj_t/concat_ref.fasta
important_outs/vdj_t/clonotypes.csv
important_outs/vdj_t/consensus.fasta
important_outs/vdj_t/filtered_contig.fastq
important_outs/vdj_t/concat_ref.bam.bai
important_outs/vdj_t/consensus_annotations.csv
important_outs/vdj_t/cell_barcodes.json
important_outs/sample_filtered_feature_bc_matrix.h5
important_outs/raw_feature_bc_matrix.h5
'''

PATHS = [
        "important_outs/raw_molecule_info.h5",
        "important_outs/sample_filtered_feature_bc_matrix.h5",
        "important_outs/raw_feature_bc_matrix.h5"
        ]

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample_name_pattern", action = "store", required = False, default = "^[0-9]+[_][0-9]+[a-z]*", help = "python re pattern for extracting the sample name from the sample tarballs")
parser.add_argument("-t", "--tarballs", action = "store", required=True, help = "Path to directory containing the sample tarballs")
parser.add_argument("-p", "--paths", action = "store", required = False, help = "Comma-separated list of paths in tarball to make uniuqely named. If none, will read list from this file.")
args = parser.parse_args()

if args.paths is None:
    paths = PATHS
else:
    paths = args.paths.strip().split(",")

tarballs = glob.glob(os.path.join(args.tarballs, "*.tar.gz"))
for tarball in tarballs:
    sn = re.search(args.sample_name_pattern, os.path.basename(tarball))
    if sn is None:
        raise ValueError(f"Unable to extract sample name from: {tarball}")
    else:
        sn = sn.group(0)
        
        with tarfile.open(tarball, "r:gz") as tf:
            logging.info(f"Opening {tarball}")
            for path in paths:
                new_name = f"{sn}_{os.path.basename(path)}"
                with open(new_name, 'wb') as renamed:
                    logging.info(f"Extracting {path} and writing to {new_name}.")
                    try:
                        maybe_buff = tf.extractfile(path)
                        if maybe_buff is not None:
                            renamed.write(maybe_buff.read())
                        else:
                            logging.warning(f"Member {path} is not a regular file. Skipping")
                    except KeyError:
                        logging.warning(f"Member {path} is not a member of archive, skipping: {tarball}")


