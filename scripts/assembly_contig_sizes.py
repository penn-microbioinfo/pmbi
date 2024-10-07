from Bio import SeqIO
import argparse

# %%
parser = argparse.ArgumentParser()
parser.add_argument("fasta", action = "store", help = "Path to FASTA file.")
args = parser.parse_args()

# %%
seqs = []
for record in SeqIO.parse(args.fasta, "fasta"):
    print(f"{record.id}\t{len(record.seq)}")
