
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

def samflag_to_dict(flag):
    dict_keys = [
            "read_paired",
            "read_mapped_in_proper_pair",
            "read_unmapped",
            "mate_unmapped",
            "read_reverse_strand",
            "mate_reverse_strand",
            "first_in_pair",
            "second_in_pair",
            "not_primary_alignment",
            "read_fails_platform_vendor_quality_checks",
            "read_is_pcr_or_optical_duplicate",
            "supplementary_alignment"
            ]
    fbr = samflag_bin_repr(flag)
    return {k:v for k,v in zip(dict_keys, fbr)}
