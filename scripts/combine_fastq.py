import re
import logging
from joblib import Parallel, delayed
import os

fastq_paths = [
        "/home/amsesk/super1/tcsl/breastvax_mouse/DSP-1001660022765/fastq/286585325/",
        "/home/amsesk/super1/tcsl/breastvax_mouse/DSP-1001660022765/fastq/286604323/"
        ]
sample_name_pat = re.compile("(DSP[-][0-9]+[-][A-Z][-][A-Z][0-9]+[_]S[0-9]+[_]L[0-9]+[_][R][0-9]_[0-9]+)[.]fastq[.]gz")
output_suffix = "fastq.gz"
outdir = "/home/amsesk/super1/tcsl/breastvax_mouse/test/combined_fastq/"

to_combine = {}
for path in fastq_paths:
    fss = os.listdir(path)
    for fs in fss:
        s = re.search(sample_name_pat, fs)
        if s is not None:
            sn = s.group(1)
            fs_full_path = os.path.join(os.path.abspath(path), fs)
            if sn in to_combine:
                to_combine[sn].append(fs_full_path)
            else:
                to_combine[sn] = [fs_full_path]
n_jobs = len(to_combine)
def _combine(prefix, file_list, outdir, suffix):
    outpath = os.path.join(outdir, f"{prefix}.{suffix}")
    if os.path.exists(outpath):
        logging.warning(f"File exists and will be overwritten: {outpath}")
    with open(outpath, 'wb') as output:
        for fs in file_list:
            with open(fs, 'rb') as input:
                output.write(input.read())
    logging.warning(f"Wrote comined file: {outpath}")
                

Parallel(n_jobs=n_jobs)(
    delayed(_combine)(prefix=k, file_list=v, outdir=outdir, suffix=output_suffix) for k,v in to_combine.items()
)

