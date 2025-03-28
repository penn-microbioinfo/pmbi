import pandas as pd
import numpy as np

# %% 
def variant_in(row, sample):
    if row[sample] is not None:
        return True
    else:
        return False

# %%
def variant_only_in(row, samples):
    variant_not_none_in = []
    for sample in row.index:
        if row[sample] is not None:
            variant_not_none_in.append(sample)
    if all([x in samples for x in variant_not_none_in]):
        return True
    else:
        return False
    
# %%
def split_snpeff_annots(annots: list):
    variant_annots = {}
    if annots is None:
        return None
    else:
        for ele in annots:
            spl = ele.split("|")
            if spl[0] in variant_annots:
                variant_annots[spl[0]].append(spl)
            else:
                variant_annots[spl[0]] = [spl]
        return variant_annots

# %%
def ad_to_af(ads, round_to=3):
    total = sum(ads)
    return [round(np.true_divide(x, total), round_to) for x in ads]
