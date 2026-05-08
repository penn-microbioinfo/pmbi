from Bio import SeqIO
from Bio.Seq import Seq
from scipy.spatial.distance import hamming
import gzip
import numpy as np

# %%
read_f = {
        "R1": "/home/amsesk/super2/jayme_shiv/data/fastq/314920619/adt_hto_maskDemux/HTO1_S25_R1_001.fastq.gz",
        "R2": "/home/amsesk/super2/jayme_shiv/data/fastq/314920619/adt_hto_maskDemux/HTO1_S25_R2_001.fastq.gz",
        "R3": "/home/amsesk/super2/jayme_shiv/data/fastq/314920619/adt_hto_maskDemux/HTO1_S25_R3_001.fastq.gz"
}

# %%
bus_str = {
  "umi": [0,0,10],
  "cbc": [1,0,16],
  "tag": [2,0,15]
}

# %%
target_bc = Seq("CTACAGACAACAAACA").reverse_complement()

# %%
def _subseq(seq,start,end):
    return seq[start:end]


# %%
umis=[]
with (
    SeqIO.parse(gzip.open(read_f["R1"], 'rt'), "fastq") as r1fh,
    SeqIO.parse(gzip.open(read_f["R2"], 'rt'), "fastq") as r2fh,
    SeqIO.parse(gzip.open(read_f["R3"], 'rt'), "fastq") as r3fh
):
    for r1,r2,r3 in zip(r1fh,r2fh,r3fh):
        start,end = bus_str["cbc"][1:3]
        subseq_len = len(range(start,end))
        assert len(target_bc)==subseq_len
        hamm_d = hamming(np.asarray(list(target_bc)), np.asarray(list(_subseq(r2.seq, start, end))))
        if hamm_d*subseq_len<=1:
            umis.append(_subseq(r1.seq, *bus_str["umi"][1:3]))

len(umis)
np.unique(np.array([str(x) for x in umis])).shape
np.asarray(umis)
# %%
# for record in SeqIO.parse("data.fastq", "fastq"):
#     print(f"ID: {record.id}")
#     print(f"Sequence: {record.seq}")
#     print(f"Quality Scores: {record.letter_annotations['phred_quality']}")
#     print("-" * 20)
["s", 1, 3, "x"]


hamming(np.asarray(list("ATCCA")), np.asarray(list("ATCCG")))*5
np.as_array("ATC")
np.asarray(list("ATC"))
list("ATC")
