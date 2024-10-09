from pmbi.scriptgen.slurm import SlurmScriptGenerator
import numpy as np
import os
import argparse
import re
import pmbi.cellranger.cellranger_command as cc
import sys
BATCH_DIR="/shared-ebs/microbioinfo-aws/scripts/batch/betts_cellranger_count_batch"

def is_empty(cmd_chunk): 
    if len(cmd_chunk) == 0:
        return True
    else:
        return False

modality_to_ref = {
    "RNA": "/cellranger-ref/refdata-gex-GRCh38-2020-A",
    "VDJ": "/cellranger-ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0",
    "ADT": "/cellranger-ref/TotalSeq_C_Human_Universal_Cocktail_399905_Antibody_reference_UMI_counting.csv" 
    }
illumina_suffix = re.compile("[_]S[0-9]+[_]L[0-9]+[_]R[0-9][_][0-9]+[.]fastq[.]gz$")
r1_pat = re.compile("[_]R1[_]")

parser = argparse.ArgumentParser()
parser.add_argument("--table", required = True)
parser.add_argument("--splits", default = 1, type = int, help = "Number of job scripts to split jobs across.")
parser.add_argument("--s3_output_prefix", default = "cellranger_matrices", help = "Where to deposit cellranger output matrices on S3.")
parser.add_argument("--multiome", action = "store_true", help = "Pass this flag if the data should be processed as a multiome.")  
args = parser.parse_args()


setup_cmds = ["source /shared-ebs/source_me.bash",
              "unset AWS_PROFILE; unset AWS_SHARED_CREDENTIALS_FILE", # Make sure boto3 uses IAM profile creds
              "mkdir -p /scratch/cellranger/fastqs",
              "cd /scratch/cellranger/fastqs"]
if args.multiome:
    with open(args.table, 'r') as table:
        table_rows = [cc.newDelimRow(line.strip(), cc.cellranger_table_column_dict(), delim = '\t') for line in table.readlines()] 
        keys = [os.path.join(row.awsdir, row.filename) for row in table_rows] 
        experiment_objects = []
        for exp in cc.parse_keys_multiome(keys):
            config = cc.CellrangerMultiConfig.from_multi_experiment("/local-ebs/cellranger/fastqs", exp, modality_to_ref)
            multi = cc.CellrangerMulti(fastqs = "/local-ebs/cellranger/fastqs", config =config, experiment = exp, output_key_prefix = args.s3_output_prefix)
            experiment_objects.append(multi)
        chunks = np.array_split(experiment_objects, args.splits) 
        print(len(chunks))

    for i,chunk in enumerate(chunks):
        if is_empty(chunk):
            raise ValueError(f"Empty command chunk passed at index {i} (0-based). Check input table or number of splits requested. There are only {len(experiment_objects)} experiments.")
        sg = SlurmScriptGenerator(
                jobname=f"crcount_{i}",
                nodes=1,
                tasks_per_node=1,
                cpus_per_task=16,
                mem=60,
                time=168,
                partition="m6id"
                )
        for sc in setup_cmds:
            sg.add_command(sc)
        for exp in chunk:
            sg.add_command(exp.print_s3_download_cmds())
            sg.add_command("cd /local-ebs/cellranger/")
            sg.add_command(exp.print_write_config_cmd())
            sg.add_command(exp.print_cellranger_cmd())
            sg.add_command(exp.print_move_outs_delete_wd())
            sg.add_command(exp.print_pull_important_outs_cmd())
            sg.add_command(exp.print_compress_outs_cmd())
            sg.add_command(exp.print_s3_upload_cmd())
        sg.write()
else:
    chunks = np.array_split(list(cc.cellranger_commands(args.table, s3_output_prefix = args.s3_output_prefix).values()), args.splits)

    for i,cmd_chunk in enumerate(chunks):
        if is_empty(cmd_chunk):
            raise ValueError(f"Empty command chunk passed at index {i}. Check input table or number of splits requested.")
        sg = SlurmScriptGenerator(
                jobname=f"crcount_{i}",
                nodes=1,
                tasks_per_node=1,
                cpus_per_task=8,
                mem=30,
                time=168,
                partition="m6id"
                )
        for sc in setup_cmds:
            sg.add_command(sc)

        for sample in cmd_chunk:
            sg.add_command(sample[0])
            sg.add_command("cd /scratch/cellranger/ && \\")
            sg.add_command(sample[1])
            sg.add_command(sample[2])
            sg.add_command(sample[3])
            sg.add_command(sample[4])
        sg.write()
