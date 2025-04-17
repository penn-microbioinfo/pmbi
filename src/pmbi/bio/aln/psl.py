import importlib

import pmbi.bio.dna
import pmbi.plotting as pmbip

importlib.reload(pmbip)

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
class PslAlignment(object):
    def __init__(self):
        self.inner = {}
        self.tSeq = None
        self.qSeq = None

    def __getattr__(self, key):
        if key in self.inner:
            return self.inner[key]
        else:
            raise KeyError(key)

    def __str__(self):
        return "\n".join([f"{k}: {v}" for k,v in self.inner.items()])

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
            inst.qStart = inst.qSize - inst.qEnd
            inst.qEnd = inst.qSize - inst.qStart
        else:
            raise ValueError(f"Invalid strand: {inst.strand}")

        return inst

    @staticmethod
    def _fields(): # {{{
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



# }}}

# %% CHUNK: PslRow{{{ # To be deprecated - replaced by PslAlignment
class PslRow(object):
    def __init__(
        self,
        row,
    ):
        spl = row.split()
        for i, field in enumerate(psl_fields):
            if "," in spl[i]:
                val = [int(x) for x in spl[i].split(",") if len(x) > 0]
            else:
                try:
                    val = int(spl[i])
                except ValueError:
                    val = spl[i]
            setattr(self, field, val)
        self.tSeq = None
        self.qSeq = None
        if self.strand == "+":
            pass
        elif self.strand == "-":
            self.qStart = self.qSize - self.qEnd
            self.qEnd = self.qSize - self.qStart
        else:
            raise ValueError(f"Invalid strand: {self.strand}")

    def __str__(self):
        return "\n".join([f"{f}: {getattr(self, f)}" for f in psl_fields])

    def add_sequences(self, tSeq, qSeq):
        self.tSeq = tSeq
        if self.strand == "+":
            self.qSeq = qSeq
        elif self.strand == "-":
            self.qSeq = pmbi.bio.dna.revcomp(qSeq)
            # self.qSeq = qSeq
        else:
            raise ValueError(f"Invalid strand: {self.strand}")


# }}}

# %% CHUNK: PslRowAln {{{ # To be deprecated - replaced by PslAlignment
class PslRowAln(object):
    def __init__(self, pslrow: PslRow, query_seqs, target_seqs):
        self.r = pslrow

    #     self.df = pd.DataFrame({
    #         "aln": [" "]*1200,
    #         "tseq": [" "]*1200,
    #         "qseq": [" "]*1200
    #         },
    #         index = range(-400,800))
    #     insert(self.df, range(0,pslrow.tSize), "tseq", self.r.tSeq)
    #     insert(self.df, self._query_insert_point(), "qseq", self.r.qSeq)
    #     for block in self._block_coords():
    #         _qs,ts,bs = block
    #         insert(self.df, range(ts,ts+bs), "aln", ["."]*bs)
    # def _query_insert_point(self):
    #     q_start_adj = self.r.tStart-self.r.qStart
    #     t_range_q = range(q_start_adj, q_start_adj+self.r.qSize)
    #     return t_range_q
    def _block_coords(self):
        block_coords = []
        for qs, ts, bs in zip(self.r.qStarts, self.r.tStarts, self.r.blockSizes):
            block_coords.append((qs, ts, bs))
        return block_coords

    def _blocks(self):
        blocks = list()
        bc = self._block_coords()
        qSize_gapped = self.r.qSize
        tSize_gapped = self.r.tSize
        for c in range(0, len(bc)):
            qS, tS, bS = bc[c]
            if c == len(bc) - 1:
                block = Block(qS, qS + bS, tS, tS + bS, bS)
            else:
                block = Block(qS, bc[c + 1][0], tS, bc[c + 1][1], bS)
            qSize_gapped += block.q_insert_bases
            tSize_gapped += block.t_insert_bases
            blocks.append(block)
        # assert sum([b.t_insert_bases for b in blocks]) == self.r.tBaseInsert
        self.qSize_gapped = qSize_gapped
        self.tSize_gapped = tSize_gapped
        return blocks

    def _aln_repr_blocks(self):
        blocks = self._blocks()
        qseq = "".join([b.q_gapped(self.r.qSeq) for b in blocks])
        tseq = "".join([b.t_gapped(self.r.tSeq) for b in blocks])
        assert len(qseq) == len(tseq)
        return (qseq, tseq)

    def _aln_track(self):
        blocks = self._blocks()
        repr = self._aln_repr()
        aln_range = range(self.leading, self.leading + len(self._aln_repr_blocks()[0]))
        track = [" "] * len(repr[0])
        for i in aln_range:
            track[i] = "."
        return "".join(track)

    def n_leading_query(self):
        blocks = self._blocks()
        first_block = blocks[0]
        return first_block.qS

    def n_trailing_query(self):
        blocks = self._blocks()
        last_block = blocks[len(blocks) - 1]
        return self.r.qSize - last_block.qE

    def leading_query_pos(self):
        first_block = self._blocks()[0]
        return range(first_block.tS - first_block.qS, first_block.tS)

    def trailing_query_pos(self):
        last_block = self._blocks()[len(self._blocks()) - 1]
        return range(last_block.tE, last_block.tE + (self.r.qSize - last_block.qE))

    def _aln_repr(self):
        qb_repr, tb_repr = self._aln_repr_blocks()
        blocks = self._blocks()
        self.leading = 0
        # Add leading parts of sequences to repr
        q_leading = self.r.qSeq[0 : blocks[0].qS]
        t_leading = self.r.tSeq[0 : blocks[0].tS]
        # Add leading spaces
        if len(q_leading) > len(t_leading):
            diff = len(q_leading) - len(t_leading)
            t_leading = " " * diff + t_leading
            self.leading += len(t_leading)
        elif len(q_leading) < len(t_leading):
            diff = len(t_leading) - len(q_leading)
            q_leading = " " * diff + q_leading
            self.leading += len(q_leading)
        else:
            # Nothing to do
            pass
        self.trailing = 0
        # Add trailing parts of sequences to repr
        q_trailing = self.r.qSeq[blocks[len(blocks) - 1].qE : self.r.qSize]
        t_trailing = self.r.tSeq[blocks[len(blocks) - 1].tE : self.r.tSize]
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
        q_repr = q_leading + qb_repr + q_trailing
        t_repr = t_leading + tb_repr + t_trailing
        return (q_repr, t_repr)

    def show(self, chunksize=100):
        qseq, tseq = self._aln_repr()
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
        return " ".join(
            map(
                str,
                [self.qS, self.qE, self.tS, self.tE, self.psl_bs, self.qbS, self.tbS],
            )
        )

    # @staticmethod
    # def _correct_ends(qS, qE, tS, tE, psl_bs, qSize, tSize):
    #     if ( (tS+(qE-qS)) > tSize ):
    #         corrected_qE = qS+psl_bs
    #     else:
    #         corrected_qE = qE
    #     if ( (qS+(tE-tS)) > qSize ):
    #         corrected_tE = tS+psl_bs
    #     else:
    #         corrected_tE = tE
    #     return(corrected_qE, corrected_tE)
    # def reaches_end_of_query(self):
    #     if self.qS+self.tbS > self.qSize:
    #         return True
    #     else:
    #         return False
    # def reaches_end_of_target(self):
    #     if self.tS+self.qbS > self.tSize:
    #         self.qE = self.qS+self.psl_bs
    #         return True
    #     else:
    #         return False
    def q_gapped(self, qSeq):
        # if self.reaches_end_of_query():
        #     self.tE = self.tS+self.psl_bs
        #     char = " "
        # else:
        #     char = "-"
        char = "-"
        qseq = qSeq[self.qS : self.qE]
        qseq += "-" * self.t_insert_bases
        return qseq

    def t_gapped(self, tSeq):
        # if self.reaches_end_of_target():
        #     char = " "
        # else:
        #     char = "-"
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
                mm.append((i, self.qS + i, self.tS + i, qseq[i], tseq[i]))
        return mm


# }}}


# %% PslAln {{{
class PslAln(object):
    def __init__(self, target_seqs, query_seqs, rows: list[PslRow] = []):
        self.rows = rows

    def from_file(stream, target_seqs, query_seqs, n=None):
        pa = PslAln(target_seqs=target_seqs, query_seqs=query_seqs)
        while not stream.readline().startswith("-"):
            continue
        for i, line in enumerate(stream):
            row = PslRow(line)
            row.add_sequences(tSeq=target_seqs[row.tName], qSeq=query_seqs[row.qName])
            pa.append(row)
            if n is not None and i == n:
                break
        return pa

    def append(self, row: PslRow):
        self.rows.append(row)


# }}}




