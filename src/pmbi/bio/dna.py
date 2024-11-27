def reverse(s):
    return s[::-1]

def revcomp(s):
    comp = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "N": "N",
            }
    return ''.join([comp[c] for c in reverse(s)])
