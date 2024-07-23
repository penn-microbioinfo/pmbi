from Bio import SeqIO
import numpy as np
import argparse

# %%
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", action = "store", help = "Path to FASTA file.")
parser.add_argument("-o", "--outfmt", action = "store", help = "", choices = ["tab", "fancy"], default = "fancy")
# parser.add_argument("-", "--", action = "store", help = "")
args = parser.parse_args()

# %%
seqs = []
for record in SeqIO.parse(args.fasta, "fasta"):
    seqs.append(record)

# %%
sorted_seqs = sorted(seqs, key = lambda x: len(x.seq))
lengths = [len(x.seq) for x in sorted_seqs]
total_len = np.sum(lengths)
running_sum = 0
n50 = None
for l in lengths:
    running_sum += l
    if running_sum >= np.divide(total_len, 2):
        n50 = l
        break


# %%
if args.outfmt == "fancy":
    print(f"Fasta file: {args.fasta}")
    print(f"Number of contigs: {len(seqs)}")
    print(f"Mean contig length: {np.mean(lengths)}")
    print(f"Longest contig: {np.max(lengths)}")
    print(f"Shortest contig: {np.min(lengths)}")
    print(f"N50: {n50}")

else:
    print('\t'.join([str(x) for x in [args.fasta, len(seqs), np.mean(lengths), np.max(lengths), np.min(lengths), n50]]))
