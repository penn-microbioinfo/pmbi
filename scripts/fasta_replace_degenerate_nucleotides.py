from Bio import SeqIO
import argparse
import io

# %%
parser =argparse.ArgumentParser()
parser.add_argument("fasta", action = "store", help = "Path to FASTA file.")
args = parser.parse_args()

degen = {
        ""
        }
# %%
f = io.StringIO(">test\nATCKMWSYRATCG")
seqs = []
for record in SeqIO.parse(args.fasta, "fasta"):
    print(f"{record.id}\t{len(record.seq)}")

# %%
SeqIO.
seqs[0].molecule_type
    seqs.append(record)

# %%

