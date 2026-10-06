import pandas as pd
import pysam
from os import PathLike

# %%
def read_bam_to_df(bam_path: PathLike) -> pd.DataFrame:
    """Read a BAM file using pysam and parse some of the important columns to a dataframe

    Args:
        bam_path: Path to the BAM file

    Returns:
        pd.DataFrame
    """
    ldict = []
    for aln in pysam.AlignmentFile(bam_path, "rb"):
        row = {
            "query": aln.query_name, 
            "reference": aln.reference_name,
            "ref_start": aln.reference_start, 
            "ref_end": aln.reference_end,
            "is_forward": aln.is_forward,
            "cigar_s": aln.cigarstring,
            "cigar_t": aln.cigartuples,
            "mq": aln.mapping_quality,
            "tags": (lambda tl: {k:v for k,v in tl})(aln.get_tags()),
        }
        ldict.append(row)
    aln_df = pd.DataFrame(ldict)
    return aln_df

# %%
def cigar_spec():
    rows = (
        ('M', 0, "alignment match (can be a sequence match or mismatch)", True, True),
        ('I', 1, "insertion to the reference", True, False),
        ('D', 2, "deletion from the reference", False, True),
        ('N', 3, "skipped region from the reference", False, True),
        ('S', 4, "soft clipping (clipped sequences present in SEQ)", True, False),
        ('H', 5, "hard clipping (clipped sequences NOT present in SEQ)", False, False),
        ('P', 6, "padding (silent deletion from padded reference)", False, False),
        ('=', 7, "sequence match", True, True),
        ('X', 8, "sequence mismatch", True, True)
    )
    columns = ("Op", "BAM", "Description", "Consumes query", "Consumes reference")
    return pd.DataFrame(rows, columns=columns)

# %%
def cigar_ref_consuming_ops(csp = cigar_spec()):
    return csp[csp["Consumes reference"]]["BAM"].to_list()

# %%
def bam_to_bed(bam_df: pd.DataFrame, file: PathLike):
    bam_df[["reference", "ref_start", "ref_end", "query"]].to_csv(file, sep="\t", index=False, header=False)



# %%
def cigar_get_ops(cigar_t: tuple[tuple[int, int]], ops: list[int]):
    return (t for t in cigar_t if t[0] in ops)


# %%
# loci_bam.iloc[0,:]
# a
# pos = a["ref_start"]
# relevent_ops = list(cigar_get_ops(a.cigar_t, cigar_ref_consuming_ops()))
# bed_d = {k: [] for k in [o[0] for o in relevent_ops]}
#
# for t in relevent_ops:
#     ls = pos
#     le = ls+t[1]+1 # add 1 to make sure ranges are only inclusive on start
#     bed_d[t[0]].append((ls, le))
#     pos = le
#
# print(bed_d)
# %%

def aligned_segment_get_gaps(s):
    blocks = s.get_blocks()
    n = len(blocks)
    i = 0
    pairs = []
    while i<(n-1):
        pairs.append( (i,i+1) )
        i+=1
    gaps = []
    for p,q in pairs:
        gaps.append( (blocks[p][1]+1, blocks[q][0]-1) )
    return gaps

# %%

