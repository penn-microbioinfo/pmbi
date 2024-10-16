from Bio import SeqIO
import argparse
import io
from numpy.random import randint

# %%
parser =argparse.ArgumentParser()
parser.add_argument("fasta", action = "store", help = "Path to FASTA file.")
args = parser.parse_args()

degen = {
        "A": ["A"],
        "G": ["G"],
        "C": ["C"],
        "T": ["T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "G", "C", "T"]
        }
# %%
for record in SeqIO.parse(args.fasta, "fasta"):
    new = ""
    for char in record.seq:
        new += degen[char.upper()]
    record.seq = new
    print(record.format("fasta"))

