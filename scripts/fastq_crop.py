# %%
import gzip
import copy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import SeqIO

# %%
pos
# %%
# %%
def count_nucleotides_at_positions(stream):
    pos_counts = {
            "A": 0,
            "T": 0,
            "C": 0,
            "G": 0,
            }
    pos = {x: copy.deepcopy(pos_counts) for x in range(0,1000)}
    for record in SeqIO.parse(stream, "fastq"):
        for id,nucl in enumerate(record.seq):
            nucl = nucl.upper()
            if nucl != "N":
                pos[id][nucl]+=1
    return pos

# %%

def calculate_props(pos):
    posdf = pd.DataFrame(pos)
    posdf_prop = posdf.apply(lambda x: x/sum(x), axis = 0).transpose()
    posdf_prop = posdf_prop[~posdf_prop.apply(lambda x: all([pd.isnull(xx) for xx in x]), axis = 1)]
    posdf_prop = posdf_prop.reset_index(names = "position")
    posdf_prop = pd.melt(posdf_prop, id_vars = "position", value_vars = list(posdf_prop.columns)[1:], var_name = "base", value_name= "proportion")
    posdf_prop = posdf_prop.sort_values("position")
    return posdf_prop

# %%
with gzip.open(
    "/home/ubuntu/projmnt/betts/coculture/raw_data/hpap/fastq/raw/17p-004-Betts-2023-05-05-human-TCRb-HPAP135-rep1-100p0ng_S46_L001_R1_001.fastq.gz",
    mode="rt",
    ) as fastq:
    pos = count_nucleotides_at_positions(fastq)
posdf_prop = calculate_props(pos)

# %%
fig,ax = plt.subplots(figsize=(15,5))
lp = sns.lineplot(posdf_prop, x = "position", y = "proportion", hue = "base", ax = ax)
fig.savefig("/srv/http/betts/coculture/basecomp_r1.pdf")
plt.clf()
# %%

