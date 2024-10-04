import argparse
import gzip
import itertools
import sys
import timeit

import numpy as np
import pandas as pd

# %%
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fastq", action="store", help="Path to index fastq file.")
# parser.add_argument("-", "--", action = "store", help = "")
# parser.add_argument("-", "--", action = "store", help = "")
args = parser.parse_args()

fastq_path = args.fastq
index_counts_dict = {
    i: 0 for i in ["".join(x) for x in list(itertools.product("ATCGN", repeat=10))]
}
# all_possible_indices = np.array(["".join(x) for x in list(itertools.product("ATCGN", repeat=10))])
# index_counts = pd.Series(
#     data=[0] * len(all_possible_indices), index=all_possible_indices
# )
# samp = np.random.randint(0,len(all_possible_indices),1000)
# samp = all_possible_indices[samp]
# samp
#
# # %%
# def with_series(ser, samp):
#     for s in samp:
#         ser.loc[s]
#
# with_series(index_counts, samp)
# # %%
# def with_dict(d, samp):
#     for s in samp:
#         x = d[s]
#
# with_dict(index_counts_dict, samp)
# # %%
#
# timeit.timeit("with_series(index_counts, samp)", setup = "from __main__ import with_series, index_counts, samp")
# # %%
#
# timeit.timeit("with_dict(index_counts_dict, samp)", setup = "from __main__ import with_dict, index_counts_dict, samp")
#
# # %%

print(index_counts_dict)
with gzip.open(fastq_path, "rt") as fastq:
    _first = fastq.readline()
    for idx, line in enumerate(fastq):
        if idx % 4 == 0:
            line = line.strip()
            index_counts_dict[line] += 1
            # index_counts.loc[line] += 1

index_counts = pd.DataFrame(
    [{"index": k, "count": v} for k, v in index_counts_dict.items()]
)
index_counts = index_counts[index_counts["count"] > 0]
index_counts = index_counts.sort_values("count", ascending=False)
print(index_counts)
# sys.exit()
# %%

import copy
import itertools

# %%
chars = list("ATCGN")
char_combs = copy.deepcopy(chars)
i = 0
new = []
while i < 10:
    if len(new) > 0:
        char_combs = new
    new = []
    for c in char_combs:
        for cc in chars:
            new.append(f"{c}{cc}")
            print(f"{c}{cc}")
    i += 1

# %%
chars = list("ATCGN")
myd = {c: {} for c in chars}


def layergen(d, stoplen=10):
    keys = list(d.keys())
    # if (keys[0]) == stoplen:
    #     return 0
    for k in keys:
        d[k] = {newkey: {} for newkey in [f"{k}{kk}" for kk in keys]}
    return d


layergen(myd)

# %%
myd


def slice_dict(d, depth=1):
    for _i in range(0, depth):
        sub = []
        for _k, v in d.items():
            sub.append(v)
    return sub


slice_dict(myd)


# %%
class Node(object):
    def __init__(self, value, parent):
        self.value = value
        self.parent = parent
        self.children = None
    def spawn_children(self, chars):
        self.children = [Node(value=c, parent=self) for c in chars]
    def spawn_generations(self, chars, ngen, current_gen=0):
        self.spawn_children(chars)
        if current_gen < ngen-1:
            for child in self.children:
                child.spawn_generations(chars, ngen, current_gen=current_gen+1)
        else:
            pass
    def get_terminal_children(self):
        if self.children is not None:
           return [c.get_terminal_children() for c in self.children]
       else:
           return None
    def to_dict(self):
        if self.children is not None:
            return {self.value: {c.value: c.to_dict() for c in self.children}}
        else:
            return None
    def __eq__(self, other):
        if self.value == other:
            return True
        else:
            return False
    def __getitem__(self, index: str):
        n = self
        for char in index:
            # child_vals = [c.value for c in n.children]
            if n.children is None or char not in n.children:
                raise KeyError
            else:
                n = {c.value: c for c in n.children}[char]
        return n


# %%
chars = "ATCGN"
n = Node("", None)
n.spawn_generations("ATCGN", 10)

# %%
n["ATTC"].children
[x.value for x in n["ATTCA"].children]
n["TTGCCCAATA"]
# %%
pprint.pprint(d)

# %%
n.goto("ATCG")
n.show_descendants()
import pprint
pprint.pprint(n.to_dict())
# %%
n.show_descendants()
[c.value for c in n.children]
n.to_dict()
# %%
