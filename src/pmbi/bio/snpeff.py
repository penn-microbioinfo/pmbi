from pathlib import Path

# %% FUNC: Class for checking and handling snpEff reference paths {{{
class snpEffRef(object):
    def __init__(self, refdir: Path):
        self.refdir = refdir
        self._check_types()
        self.sequences, self.genes, self.cds, self.proteins = self._parse_refdir(refdir)
    @staticmethod
    def _parse_refdir(refdir):
        expected_files = ["sequences.fa.gz",
                          "genes.gtf.gz",
                          "cds.fa.gz",
                          "protein.fa.gz"
                          ]
        fs = [f.name for f in refdir.iterdir()]
        if any([e not in fs for e in expected_files]):
            raise OSError("Missing files in snpEff reference directory")
        else:
            return tuple([refdir.joinpath(e) for e in expected_files])
    def _check_types(self):
        if not isinstance(self.refdir, Path):
            raise TypeError("refdir is not of type Path")

# }}}
