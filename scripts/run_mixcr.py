import argparse
import subprocess
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fastq", action = "store", help = "Path to directory containing fastq files.")
parser.add_argument("-s", "--sample_id_pattern", action = "store", help = "Python regex expression to pull sameple ids from fastq filenames.")
parser.add_argument("-p", "--preset", action = "store", help = "`mixcp analyze` preset", default = "takara-human-rna-tcr-smarter")
# parser.add_argument("-o", "--output", action = "store", help = "Prefix for output")
args = parser.parse_args()

fastqs = [os.path.join(args.fastq, f) for f in os.listdir(args.fastq) if f.endswith(".fastq") or f.endswith(".fastq.gz")]
matches = [re.search(args.sample_id_pattern, os.path.basename(f)) for f in fastqs]
if any([m is None for m in matches]): 
    raise ValueError("Sample id pattern not working")
else:
    sample_ids = [m.group(1) for m in matches]
    for sid in set(sample_ids):
        print(sid)
        sample_fastqs = [f for f in fastqs if sid in os.path.basename(f)]
        sample_fastqs.sort()
        cmd = [
                "mixcr",
                "analyze",
                args.preset,
                sample_fastqs[0], sample_fastqs[1],
                sid,
                "-f",
                "--assemble-clonotypes-by", "CDR3"
                ]
        print(' '.join([str(x) for x in cmd]))
        subprocess.run(cmd)
