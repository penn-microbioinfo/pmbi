from __future__ import annotations

import importlib
import io
from collections import namedtuple

import pmbi.bio.dna
import pmbi.plotting as pmbip

import numpy as np
import pandas as pd

importlib.reload(pmbip)

Mismatch = namedtuple(
    "Mismatch", ["blockPos", "queryPos", "targetPos", "queryBase", "targetBase"]
)

# %% CHUNK: PSL fields {{{
psl_fields = [
    "matches",
    "misMatches",
    "repMatches",
    "nCount",
    "qNumInsert",
    "qBaseInsert",
    "tNumInsert",
    "tBaseInsert",
    "strand",
    "qName",
    "qSize",
    "qStart",
    "qEnd",
    "tName",
    "tSize",
    "tStart",
    "tEnd",
    "blockCount",
    "blockSizes",
    "qStarts",
    "tStarts",
]


# }}}
# A class that represents a row of a PSL alignment file
# Combines and replaces PslRow and PslRowAln
# %% CHUNK: PslAlignment {{{
class PslAlignment:
    def __init__(self):
        # self.inner = {}
        super().__setattr__("inner", {})

    def __getattr__(self, key):
        if key in self.inner:
            return self.inner[key]
        else:
            raise KeyError(key)

    def __setattr__(self, key, value):
        super().__getattribute__("inner")[key] = value

    def __str__(self):
        return "\n".join([f"{k}: {v}" for k, v in self.inner.items()])

    def add_sequences(self, tSeq, qSeq):
        self.inner["tSeq"] = tSeq
        if self.strand == "+":
            self.inner["qSeq"] = qSeq
        elif self.strand == "-":
            self.inner["qSeq"] = pmbi.bio.dna.revcomp(qSeq)
        else:
            raise ValueError(f"Invalid strand: {self.strand}")

    @staticmethod
    # Some columns have commas in them. If so, split them into int lists. They should be lists of integers.
    # If there is no comma, try to convert the column value into int.
    # Finally just use whatever the value is.
    def from_str(string):
        inst = PslAlignment()
        spl = string.split()
        for i, field in enumerate(PslAlignment._fields()):
            if "," in spl[i]:
                val = [int(x) for x in spl[i].split(",") if len(x) > 0]
            else:
                try:
                    val = int(spl[i])
                except ValueError:
                    val = spl[i]
            inst.inner[field] = val

        if inst.strand == "+":
            pass
        elif inst.strand == "-":
            qStart_orig = inst.qStart
            qEnd_orig = inst.qEnd
            inst.inner["qStart"] = inst.inner["qSize"] - qEnd_orig
            inst.inner["qEnd"] = inst.inner["qSize"] - qStart_orig
        else:
            raise ValueError(f"Invalid strand: {inst.strand}")

        return inst

    @staticmethod
    def _fields():  # {{{
        return [
            "matches",
            "misMatches",
            "repMatches",
            "nCount",
            "qNumInsert",
            "qBaseInsert",
            "tNumInsert",
            "tBaseInsert",
            "strand",
            "qName",
            "qSize",
            "qStart",
            "qEnd",
            "tName",
            "tSize",
            "tStart",
            "tEnd",
            "blockCount",
            "blockSizes",
            "qStarts",
            "tStarts",
        ]
        # }}}

    def _block_coords(self):
        block_coords = []
        for qs, ts, bs in zip(self.qStarts, self.tStarts, self.blockSizes):
            block_coords.append((qs, ts, bs))
        return block_coords

    def _blocks(self):
        blocks = list()
        bc = self._block_coords()
        qSize_gapped = self.qSize
        tSize_gapped = self.tSize
        for c in range(0, len(bc)):
            qS, tS, bS = bc[c]
            if c == len(bc) - 1:
                block = Block(qS, qS + bS, tS, tS + bS, bS)
            else:
                block = Block(qS, bc[c + 1][0], tS, bc[c + 1][1], bS)
            qSize_gapped += block.q_insert_bases
            tSize_gapped += block.t_insert_bases
            blocks.append(block)
        # assert sum([b.t_insert_bases for b in blocks]) == self.tBaseInsert
        self.qSize_gapped = qSize_gapped
        self.tSize_gapped = tSize_gapped
        return blocks

    def _aln_repr_blocks(self):
        blocks = self._blocks()
        qseq = "".join([b.q_gapped(self.qSeq) for b in blocks])
        tseq = "".join([b.t_gapped(self.tSeq) for b in blocks])
        assert len(qseq) == len(tseq)
        return (qseq, tseq)

    def _aln_track(self):
        # blocks = self._blocks()
        _tlims, q_repr, t_repr  = self._aln_repr()
        aln_range = range(self.leading, self.leading + len(self._aln_repr_blocks()[0]))
        track = [" "] * len(q_repr)
        for i in aln_range:
            track[i] = "."
        return "".join(track)

    def n_leading_query(self):
        return self.qStart

    def n_trailing_query(self):
        return self.qSize - self.qEnd

    # def n_trailing_query(self):
    #     blocks = self._blocks()
    #     last_block = blocks[len(blocks) - 1]
    #     return self.qSize - last_block.qE

    def leading_query_pos(self):
        first_block = self._blocks()[0]
        return range(first_block.tS - first_block.qS, first_block.tS)

    def trailing_query_pos(self):
        last_block = self._blocks()[len(self._blocks()) - 1]
        return range(last_block.tE, last_block.tE + (self.qSize - last_block.qE))

    def _aln_leading(self):
        qb_repr, tb_repr = self._aln_repr_blocks()
        blocks = self._blocks()
        leading = 0
        # Add leading parts of sequences to repr
        q_leading = self.qSeq[0 : blocks[0].qS]
        t_leading = self.tSeq[0 : blocks[0].tS]
        # Add leading spaces
        if len(q_leading) > len(t_leading):
            diff = len(q_leading) - len(t_leading)
            t_leading = " " * diff + t_leading
            leading += len(t_leading)
        elif len(q_leading) < len(t_leading):
            diff = len(t_leading) - len(q_leading)
            q_leading = " " * diff + q_leading
            leading += len(q_leading)
        else:
            # Nothing to do
            pass
        assert( all([b==' ' for b in q_leading]) or [b==' 'for b in t_leading] )
        return (q_leading, t_leading)

    def _aln_trailing(self):
        trailing = 0
        # Add trailing parts of sequences to repr
        blocks = self._blocks()
        print(blocks[len(blocks) - 1].qE, self.qSize)
        q_trailing = self.qSeq[blocks[len(blocks) - 1].qE : self.qSize]

        print(blocks[len(blocks) - 1].tE, self.tSize)
        t_trailing = self.tSeq[blocks[len(blocks) - 1].tE : self.tSize]
        # Add trailing spaces
        if len(q_trailing) > len(t_trailing):
            diff = len(q_trailing) - len(t_trailing)
            t_trailing = t_trailing + " " * diff
            self.trailing = len(t_trailing)
        elif len(q_trailing) < len(t_trailing):
            diff = len(t_trailing) - len(q_trailing)
            q_trailing = q_trailing + " " * diff
            self.trailing = len(q_trailing)
        else:
            # Nothing to do
            pass
        assert( all([b==' ' for b in q_trailing]) or [b==' 'for b in t_trailing] )
        return (q_trailing, t_trailing)

    def _limits(self):
        q_leading, t_leading = self._aln_leading()
        q_trailing, t_trailing = self._aln_trailing()
    def _aln_repr(self):
        qb_repr, tb_repr = self._aln_repr_blocks()
        blocks = self._blocks()
        self.leading = 0
        tlims = {"min": 0, "max": self.tSize}
        # Add leading parts of sequences to repr
        q_leading = self.qSeq[0 : blocks[0].qS]
        t_leading = self.tSeq[0 : blocks[0].tS]
        # Add leading spaces
        if len(q_leading) > len(t_leading):
            diff = len(q_leading) - len(t_leading)
            t_leading = " " * diff + t_leading
            self.leading += len(t_leading)
            tlims["min"] =  0 - diff
        elif len(q_leading) < len(t_leading):
            diff = len(t_leading) - len(q_leading)
            q_leading = " " * diff + q_leading
            self.leading += len(q_leading)
        else:
            # Nothing to do
            pass
        self.trailing = 0
        # Add trailing parts of sequences to repr
        q_trailing = self.qSeq[blocks[len(blocks) - 1].qE : self.qSize]
        t_trailing = self.tSeq[blocks[len(blocks) - 1].tE : self.tSize]
        # Add trailing spaces
        if len(q_trailing) > len(t_trailing):
            diff = len(q_trailing) - len(t_trailing)
            t_trailing = t_trailing + " " * diff
            self.trailing = len(t_trailing)
            tlims["max"] = self.tSize + diff
        elif len(q_trailing) < len(t_trailing):
            diff = len(t_trailing) - len(q_trailing)
            q_trailing = q_trailing + " " * diff
            self.trailing = len(q_trailing)
        else:
            # Nothing to do
            pass
        q_repr = q_leading + qb_repr + q_trailing
        t_repr = t_leading + tb_repr + t_trailing
        return (tlims, q_repr, t_repr)

    def show(self, chunksize=100):
        _tlims, qseq, tseq = self._aln_repr()
        assert len(qseq) == len(tseq)
        breaks = list(range(0, len(tseq), chunksize))[0 : len(tseq)]
        lines = []
        for b in breaks:
            lines.append(
                "\n".join(
                    [
                        f"-- {b}:{b+chunksize}",
                        f"A: {self._aln_track()[b:b+chunksize]}",
                        f"Q: {qseq[b:b+chunksize]}",
                        f"T: {tseq[b:b+chunksize]}",
                    ]
                )
            )
        print("\n".join(lines))

# }}}

# %% CHUNK: Block class {{{
class Block(object):
    def __init__(self, qS, qE, tS, tE, psl_bs):
        self.qS = qS
        self.tS = tS
        # self.qE, self.tE = Block._correct_ends(qS, qE, tS, tE, psl_bs, qSize, tSize)
        self.qE, self.tE = (qE, tE)
        # self.qSize = qSize
        # self.tSize = tSize
        self.psl_bs = psl_bs
        self.qbS = self.qE - qS
        self.tbS = self.tE - tS
        assert self.qbS >= psl_bs and self.tbS >= psl_bs
        self.q_has_insert = self.qbS > psl_bs
        self.t_has_insert = self.tbS > psl_bs
        self.q_insert_bases = self.qbS - psl_bs
        self.t_insert_bases = self.tbS - psl_bs

    def __str__(self):
        return f"""
queryStart: {self.qS}
queryEnd: {self.qE}
targetStart: {self.tS}
targetEnd: {self.tE}
PSL block size: {self.psl_bs}
Query block size: {self.qbS}
Target block size: {self.tbS}
            """

    def q_gapped(self, qSeq):
        char = "-"
        qseq = qSeq[self.qS : self.qE]
        qseq += char * self.t_insert_bases
        return qseq

    def t_gapped(self, tSeq):
        char = "-"
        tseq = tSeq[self.tS : self.tE]
        tseq += char * self.q_insert_bases
        return tseq

    def mismatches(self, qSeq, tSeq):
        qseq = qSeq[self.qS : self.qE]
        tseq = tSeq[self.tS : self.tE]
        mm = []
        for i in range(0, min([len(qseq), len(tseq)])):
            if qseq[i] != tseq[i]:
                # mm.append((i, self.qS + i, self.tS + i, qseq[i], tseq[i]))
                mm.append(
                    Mismatch(
                        blockPos=i,
                        queryPos=self.qS + i,
                        targetPos=self.tS + i,
                        queryBase=qseq[i],
                        targetBase=tseq[i],
                    )
                )
        return mm

# %% PslAln {{{
class PslAlignmentCollection(object):
    def __init__(self, rows: list[PslAlignment]|None = None):
        if rows is None:
            self.rows = []
        else:
            self.rows = rows

    @staticmethod
    def from_file(
        stream: io.TextIOWrapper,
        target_seqs: dict,
        query_seqs: dict,
        n: int | None = None,
        targetSeq_include: list[str] | None = None,
    ) -> PslAlignmentCollection:
        pa = PslAlignmentCollection()
        while not stream.readline().startswith("-"):
            continue
        i = 0
        for line in stream:
            row = PslAlignment.from_str(string=line)
            row.add_sequences(tSeq=target_seqs[row.tName], qSeq=query_seqs[row.qName])
            if targetSeq_include is not None:
                if row.tName in targetSeq_include:
                    pa.append(row)
                    i = i+1
            else:
                pa.append(row)
                i = i+1
            if n is not None and i == n:
                break
        return pa

    def append(self, row: PslAlignment):
        self.rows.append(row)

    def filter(self, f, *args, **kwargs):
        filt_rows = []
        for r in self.rows:
            if f(r, *args, **kwargs):
                filt_rows.append(r)

        return PslAlignmentCollection(rows=filt_rows)

    def select(self, attrs):
        ldict = []
        for r in self.rows:
            ldict.append({attr: getattr(r, attr) for attr in attrs})

        return pd.DataFrame(ldict)

    def as_dataframe(self):
        return pd.DataFrame([r.inner for r in self.rows])

class TargetAlignments(object):
    def __init__(self, pac: PslAlignmentCollection, target: str):
        self.rows = [p for p in pac.rows if p.tName == target]
        self.target = target
        assert all([p.tName==self.target for p in self.rows])
        tSeq_lens = [len(x.tSeq) for x in self.rows]
        assert all([x==tSeq_lens[0] for x in tSeq_lens])
        assert(len(set([x.tSeq for x in self.rows]))==1)
        self.tSeq = self.rows[0].tSeq

    def target_coverage(self):
        arrs = list()
        for p in self.rows:
            arr = np.array([0]*len(self.tSeq))
            arr[p.tStart:p.tEnd] = 1
            arrs.append(arr)
        return np.apply_along_axis(sum, axis=0, arr=np.stack(arrs))

        






# }}}

# %%
# # %% CHUNK: PslRow{{{ # To be deprecated - replaced by PslAlignment
# class PslRow(object):
#     def __init__(
#         self,
#         row,
#     ):
#         spl = row.split()
#         for i, field in enumerate(psl_fields):
#             if "," in spl[i]:
#                 val = [int(x) for x in spl[i].split(",") if len(x) > 0]
#             else:
#                 try:
#                     val = int(spl[i])
#                 except ValueError:
#                     val = spl[i]
#             setattr(self, field, val)
#         self.tSeq = None
#         self.qSeq = None
#         if self.strand == "+":
#             pass
#         elif self.strand == "-":
#             self.qStart = self.qSize - self.qEnd
#             self.qEnd = self.qSize - self.qStart
#         else:
#             raise ValueError(f"Invalid strand: {self.strand}")
#
#     def __str__(self):
#         return "\n".join([f"{f}: {getattr(self, f)}" for f in psl_fields])
#
#     def add_sequences(self, tSeq, qSeq):
#         self.tSeq = tSeq
#         if self.strand == "+":
#             self.qSeq = qSeq
#         elif self.strand == "-":
#             self.qSeq = pmbi.bio.dna.revcomp(qSeq)
#             # self.qSeq = qSeq
#         else:
#             raise ValueError(f"Invalid strand: {self.strand}")
#
#     def _block_coords(self):
#         block_coords = []
#         for qs, ts, bs in zip(self.r.qStarts, self.r.tStarts, self.r.blockSizes):
#             block_coords.append((qs, ts, bs))
#         return block_coords
#
#     def _blocks(self):
#         blocks = list()
#         bc = self._block_coords()
#         qSize_gapped = self.r.qSize
#         tSize_gapped = self.r.tSize
#         for c in range(0, len(bc)):
#             qS, tS, bS = bc[c]
#             if c == len(bc) - 1:
#                 block = Block(qS, qS + bS, tS, tS + bS, bS)
#             else:
#                 block = Block(qS, bc[c + 1][0], tS, bc[c + 1][1], bS)
#             qSize_gapped += block.q_insert_bases
#             tSize_gapped += block.t_insert_bases
#             blocks.append(block)
#         # assert sum([b.t_insert_bases for b in blocks]) == self.r.tBaseInsert
#         self.qSize_gapped = qSize_gapped
#         self.tSize_gapped = tSize_gapped
#         return blocks
#
#     def _aln_repr_blocks(self):
#         blocks = self._blocks()
#         qseq = "".join([b.q_gapped(self.r.qSeq) for b in blocks])
#         tseq = "".join([b.t_gapped(self.r.tSeq) for b in blocks])
#         assert len(qseq) == len(tseq)
#         return (qseq, tseq)
#
#     def _aln_track(self):
#         blocks = self._blocks()
#         repr = self._aln_repr()
#         aln_range = range(self.leading, self.leading + len(self._aln_repr_blocks()[0]))
#         track = [" "] * len(repr[0])
#         for i in aln_range:
#             track[i] = "."
#         return "".join(track)
#
#     def n_leading_query(self):
#         blocks = self._blocks()
#         first_block = blocks[0]
#         return first_block.qS
#
#     def n_trailing_query(self):
#         blocks = self._blocks()
#         last_block = blocks[len(blocks) - 1]
#         return self.r.qSize - last_block.qE
#
#     def leading_query_pos(self):
#         first_block = self._blocks()[0]
#         return range(first_block.tS - first_block.qS, first_block.tS)
#
#     def trailing_query_pos(self):
#         last_block = self._blocks()[len(self._blocks()) - 1]
#         return range(last_block.tE, last_block.tE + (self.r.qSize - last_block.qE))
#
#     def _aln_repr(self):
#         qb_repr, tb_repr = self._aln_repr_blocks()
#         blocks = self._blocks()
#         self.leading = 0
#         # Add leading parts of sequences to repr
#         q_leading = self.r.qSeq[0 : blocks[0].qS]
#         t_leading = self.r.tSeq[0 : blocks[0].tS]
#         # Add leading spaces
#         if len(q_leading) > len(t_leading):
#             diff = len(q_leading) - len(t_leading)
#             t_leading = " " * diff + t_leading
#             self.leading += len(t_leading)
#         elif len(q_leading) < len(t_leading):
#             diff = len(t_leading) - len(q_leading)
#             q_leading = " " * diff + q_leading
#             self.leading += len(q_leading)
#         else:
#             # Nothing to do
#             pass
#         self.trailing = 0
#         # Add trailing parts of sequences to repr
#         q_trailing = self.r.qSeq[blocks[len(blocks) - 1].qE : self.r.qSize]
#         t_trailing = self.r.tSeq[blocks[len(blocks) - 1].tE : self.r.tSize]
#         # Add trailing spaces
#         if len(q_trailing) > len(t_trailing):
#             diff = len(q_trailing) - len(t_trailing)
#             t_trailing = t_trailing + " " * diff
#             self.trailing = len(t_trailing)
#         elif len(q_trailing) < len(t_trailing):
#             diff = len(t_trailing) - len(q_trailing)
#             q_trailing = q_trailing + " " * diff
#             self.trailing = len(q_trailing)
#         else:
#             # Nothing to do
#             pass
#         q_repr = q_leading + qb_repr + q_trailing
#         t_repr = t_leading + tb_repr + t_trailing
#         return (q_repr, t_repr)
#
#     def show(self, chunksize=100):
#         qseq, tseq = self._aln_repr()
#         assert len(qseq) == len(tseq)
#         breaks = list(range(0, len(tseq), chunksize))[0 : len(tseq)]
#         lines = []
#         for b in breaks:
#             lines.append(
#                 "\n".join(
#                     [
#                         f"-- {b}:{b+chunksize}",
#                         f"A: {self._aln_track()[b:b+chunksize]}",
#                         f"Q: {qseq[b:b+chunksize]}",
#                         f"T: {tseq[b:b+chunksize]}",
#                     ]
#                 )
#             )
#         print("\n".join(lines))
#
#
#
# # }}}
#
# # %% CHUNK: PslRowAln {{{ # To be deprecated - replaced by PslAlignment
# class PslRowAln(object):
#     def __init__(self, pslrow: PslRow, query_seqs, target_seqs):
#         self.r = pslrow
#
#     #     self.df = pd.DataFrame({
#     #         "aln": [" "]*1200,
#     #         "tseq": [" "]*1200,
#     #         "qseq": [" "]*1200
#     #         },
#     #         index = range(-400,800))
#     #     insert(self.df, range(0,pslrow.tSize), "tseq", self.r.tSeq)
#     #     insert(self.df, self._query_insert_point(), "qseq", self.r.qSeq)
#     #     for block in self._block_coords():
#     #         _qs,ts,bs = block
#     #         insert(self.df, range(ts,ts+bs), "aln", ["."]*bs)
#     # def _query_insert_point(self):
#     #     q_start_adj = self.r.tStart-self.r.qStart
#     #     t_range_q = range(q_start_adj, q_start_adj+self.r.qSize)
#     #     return t_range_q
#     def _block_coords(self):
#         block_coords = []
#         for qs, ts, bs in zip(self.r.qStarts, self.r.tStarts, self.r.blockSizes):
#             block_coords.append((qs, ts, bs))
#         return block_coords
#
#     def _blocks(self):
#         blocks = list()
#         bc = self._block_coords()
#         qSize_gapped = self.r.qSize
#         tSize_gapped = self.r.tSize
#         for c in range(0, len(bc)):
#             qS, tS, bS = bc[c]
#             if c == len(bc) - 1:
#                 block = Block(qS, qS + bS, tS, tS + bS, bS)
#             else:
#                 block = Block(qS, bc[c + 1][0], tS, bc[c + 1][1], bS)
#             qSize_gapped += block.q_insert_bases
#             tSize_gapped += block.t_insert_bases
#             blocks.append(block)
#         # assert sum([b.t_insert_bases for b in blocks]) == self.r.tBaseInsert
#         self.qSize_gapped = qSize_gapped
#         self.tSize_gapped = tSize_gapped
#         return blocks
#
#     def _aln_repr_blocks(self):
#         blocks = self._blocks()
#         qseq = "".join([b.q_gapped(self.r.qSeq) for b in blocks])
#         tseq = "".join([b.t_gapped(self.r.tSeq) for b in blocks])
#         assert len(qseq) == len(tseq)
#         return (qseq, tseq)
#
#     def _aln_track(self):
#         blocks = self._blocks()
#         repr = self._aln_repr()
#         aln_range = range(self.leading, self.leading + len(self._aln_repr_blocks()[0]))
#         track = [" "] * len(repr[0])
#         for i in aln_range:
#             track[i] = "."
#         return "".join(track)
#
#     def n_leading_query(self):
#         blocks = self._blocks()
#         first_block = blocks[0]
#         return first_block.qS
#
#     def n_trailing_query(self):
#         blocks = self._blocks()
#         last_block = blocks[len(blocks) - 1]
#         return self.r.qSize - last_block.qE
#
#     def leading_query_pos(self):
#         first_block = self._blocks()[0]
#         return range(first_block.tS - first_block.qS, first_block.tS)
#
#     def trailing_query_pos(self):
#         last_block = self._blocks()[len(self._blocks()) - 1]
#         return range(last_block.tE, last_block.tE + (self.r.qSize - last_block.qE))
#
#     def _aln_repr(self):
#         qb_repr, tb_repr = self._aln_repr_blocks()
#         blocks = self._blocks()
#         self.leading = 0
#         # Add leading parts of sequences to repr
#         q_leading = self.r.qSeq[0 : blocks[0].qS]
#         t_leading = self.r.tSeq[0 : blocks[0].tS]
#         # Add leading spaces
#         if len(q_leading) > len(t_leading):
#             diff = len(q_leading) - len(t_leading)
#             t_leading = " " * diff + t_leading
#             self.leading += len(t_leading)
#         elif len(q_leading) < len(t_leading):
#             diff = len(t_leading) - len(q_leading)
#             q_leading = " " * diff + q_leading
#             self.leading += len(q_leading)
#         else:
#             # Nothing to do
#             pass
#         self.trailing = 0
#         # Add trailing parts of sequences to repr
#         q_trailing = self.r.qSeq[blocks[len(blocks) - 1].qE : self.r.qSize]
#         t_trailing = self.r.tSeq[blocks[len(blocks) - 1].tE : self.r.tSize]
#         # Add trailing spaces
#         if len(q_trailing) > len(t_trailing):
#             diff = len(q_trailing) - len(t_trailing)
#             t_trailing = t_trailing + " " * diff
#             self.trailing = len(t_trailing)
#         elif len(q_trailing) < len(t_trailing):
#             diff = len(t_trailing) - len(q_trailing)
#             q_trailing = q_trailing + " " * diff
#             self.trailing = len(q_trailing)
#         else:
#             # Nothing to do
#             pass
#         q_repr = q_leading + qb_repr + q_trailing
#         t_repr = t_leading + tb_repr + t_trailing
#         return (q_repr, t_repr)
#
#     def show(self, chunksize=100):
#         qseq, tseq = self._aln_repr()
#         assert len(qseq) == len(tseq)
#         breaks = list(range(0, len(tseq), chunksize))[0 : len(tseq)]
#         lines = []
#         for b in breaks:
#             lines.append(
#                 "\n".join(
#                     [
#                         f"-- {b}:{b+chunksize}",
#                         f"A: {self._aln_track()[b:b+chunksize]}",
#                         f"Q: {qseq[b:b+chunksize]}",
#                         f"T: {tseq[b:b+chunksize]}",
#                     ]
#                 )
#             )
#         print("\n".join(lines))
#
#
# # }}}
#
#
#
#
