
# %%
def samflag_bin_repr(flag):
    """
    Convert a SAM flag to its binary string representation.

    The function takes an integer SAM flag and converts it to a 
    binary string representation. Each bit in the binary representation denotes 
    a specific property of the alignment. See https://broadinstitute.github.io/picard/explain-flags.html
    for help.

    Parameters:
    flag (int): An integer SAM flag, typically derived from sequence alignment data.

    Returns:
    str: A binary string representation of the input SAM flag.
    """
    return f"{flag:0>12b}"
