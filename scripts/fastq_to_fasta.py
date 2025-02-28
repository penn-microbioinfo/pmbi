from Bio import SeqIO
import numpy as np
import argparse

# %%
parser = argparse.ArgumentParser()
parser.add_argument("-fq", "--fastq", action = "store", help = "Path to input FASTQ file.")
parser.add_argument("-fa", "--fasta", action = "store", help = "Path to output FASTA file.")
args = parser.parse_args()

# %%
seqs = []
n_zero_len = 0
n_reads = 0
with open(args.fasta, 'w') as out:
    for record in SeqIO.parse(args.fastq, "fastq"):
        if len(record.seq) == 0:
            n_zero_len+=1
        else:
            out.write(record.format("fasta"))
        n_reads+=1
print(f"Number of reads: {n_reads}")
print(f"Reads with zero-len sequence: {n_zero_len}")
        
