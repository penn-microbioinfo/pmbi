from Bio import SeqIO
import numpy as np
import argparse

# %%
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", action = "store", help = "Path to FASTA file.")
args = parser.parse_args()

# %%
seqs = []
for record in SeqIO.parse(args.fasta, "fasta"):
    print(len(record.seq))

