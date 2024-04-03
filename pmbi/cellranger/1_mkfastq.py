from pmbi.scriptgen.slurm import SlurmScriptGenerator
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--run", action = "store", help = "Path to Illumina run directory.")
parser.add_argument("-c", "--csv", action = "store", help = "Path to the sample sheet csv.")
parser.add_argument("-o", "--output_dir", action = "store", help = "Path to write out fastq files.")
parser.add_argument("-t", "--threads", action = "store", help = "Threads for cellranger to use.")
parser.add_argument("-m", "--memory", action = "store", help = "Memory (Gb) for cellranger to use.")
args = parser.parse_args()

sg = SlurmScriptGenerator(
        jobname = "crmkfastq",
        cpus_per_task = args.threads,
        mem = args.memory,
        time = 48,
        partition = "m6id"
        )

sg.add_command("source /shared-ebs/source_me.bash")
sg.add_command(f"cellranger mkfastq --run={args.run} --csv={args.csv} --output-dir={args.output_dir} --jobmode=local --localcores={args.threads} --localmem={args.memory}")

sg.write()
