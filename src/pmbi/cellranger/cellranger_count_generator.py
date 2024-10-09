import scriptgen
import numpy as np
import os
import argparse
import re

illumina_suffix = re.compile("[_]S[0-9]+[_]L[0-9]+[_]R[0-9][_][0-9]+[.]fastq[.]gz$")
r1_pat = re.compile("[_]R1[_]")

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--fastqs", help = "Path to fastq directory.")
group.add_argument("--commands", help = "Path to text file containing commands to run. One per line.")
parser.add_argument("--splits", default = 1, type = int, help = "Number of job scripts to split jobs across.")
args = parser.parse_args()

if args.fastqs is not None:
    r1 = [x for x in os.listdir(args.fastqs) if x.endswith(".fastq.gz") and re.search(r1_pat, x) is not None] 
    samples = np.array([re.sub(illumina_suffix, "", x) for x in r1]) 

    for i,chunk in enumerate(np.array_split(samples, args.splits)):
        sg = scriptgen.SlurmScriptGenerator(
                jobname=f"crcount_{i}",
                nodes=1,
                tasks_per_node=1,
                cpus_per_task=8,
                mem=64,
                time=168,
                partition="r6a"
                )
        for sample in chunk:
            sg.add_command(f"cellranger count --id {sample} --fastqs {args.fastqs} --sample {sample} --transcriptome /scratch/ref/refdata-gex-GRCh38-2020-A --") 
else:

    for i,cmd_chunk in enumerate(np.array_split(open(args.commands, 'r').readlines(), args.splits)):
        sg = scriptgen.SlurmScriptGenerator(
                jobname=f"crcount_{i}",
                nodes=1,
                tasks_per_node=1,
                cpus_per_task=8,
                mem=64,
                time=168,
                partition="r6a"
                )
        for cmd in cmd_chunk:
            sg.add_command(cmd)

sg.write()
