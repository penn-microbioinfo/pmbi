#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:23:26 2020

@author: aimzez
"""
from __future__ import annotations

import ast
import os
import re
from pathlib import Path
from typing import Callable, Union

import bs4
import numpy as np
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup as BS

import pmbi.logging as plog
import pmbi.subproc as pmbiproc

logger = plog.streamLogger(name=__name__)


class NcbiNucleotideRecords(object):

    def __init__(self):
        pass

    def from_keyword(soup):
        pass


class NcbiSraExperimentPackageSet(object):

    def __init__(self, index):
        self.index = index

    def __iter__(self):
        yield from self.index.items()

    def __len__(self):
        return len(self.index)

    def subset(self, keys: list[str]) -> NcbiSraExperimentPackageSet:
        for k in keys:
            if k not in self.index.keys():
                raise KeyError(f"invalid key: {k}")
        return NcbiSraExperimentPackageSet(
            index={k: v for k, v in self.index.items() if k in keys}
        )

    def isubset(self, indices: list[int]) -> NcbiSraExperimentPackageSet:
        if max(indices) > (len(self.index) - 1):
            raise KeyError("indices out of range")

        keys = [k for i, k in enumerate(self.index.keys()) if i in indices]
        return self.subset(keys=keys)

    def select(self, key: str) -> NcbiSraExperimentPackage:
        return self.index[key]

    def iselect(self, index: int) -> NcbiSraExperimentPackgeSet:
        return self.index[list(self.index.keys())[index]]

    def from_soup(soup, delim="EXPERIMENT_PACKAGE"):
        index = dict()
        for run in soup.find_all(delim):

            acc = run.find("EXPERIMENT").attrs["accession"]

            if acc not in index:
                index[acc] = NcbiSraExperimentPackage(acc, run)
            else:
                raise KeyError(f"Experiment {exp} duplicated.")

        return NcbiSraExperimentPackageSet(index)

    def from_bioproject(bioproject, delim="EXPERIMENT_PACKAGE"):
        handle = Entrez.esearch(db="sra", term=bioproject, retmax=1000)

        records = Entrez.read(handle)
        sra_ids = records["IdList"]
        handle.close()

        handle = Entrez.efetch(db="sra", id=sra_ids)
        sra_records = BS(handle.read(), "xml")
        handle.close()

        sra_accessions = [
            record.attrs["accession"] for record in sra_records("EXPERIMENT")
        ]

        handle = Entrez.efetch(db="sra", id=sra_accessions)
        sra_records = BS(handle.read(), "xml")
        handle.close()

        return NcbiSraExperimentPackageSet.from_soup(sra_records)

    # def filter(self, expr):
    #     s = re.search("(.*)(==|!=)(.*)", expr)
    #     if s is None:
    #         raise ValueError("Invalid filtering expression.")
    #
    #     else:
    #         attr, operand, value = (x.strip() for x in s.groups())
    #
    #     new_index = dict()
    #     for acc, run in self:
    #         if not hasattr(run, attr):
    #             raise ValueError("Invalid attribute.")
    #
    #         if operand == "==":
    #             if getattr(run, attr) == value:
    #                 new_index[acc] = run
    #         elif operand == "!=":
    #             if getattr(run, attr) != value:
    #                 new_index[acc] = run
    #
    #     return NcbiSraExperimentPackageSet(new_index)

    # def fetch_metadata(self):
    #     biosamples = [x.BioSample for _, x in self]
    #     handle = Entrez.efetch(db="biosample", id=biosamples, maxret=100)
    #     sample_soup = BS(handle.read(), "xml")
    #     handle.close()
    #
    #     for _, run in self:
    #         s = sample_soup.select(f'BioSample[accession={run.BioSample}]')
    #         s = s[0]
    #         meta_dict = dict()
    #         for m in self.META:
    #             try:
    #                 meta_dict[m] = s.select(f'Attribute[harmonized_name="{m}"]')[0].string
    #             except IndexError:
    #                 meta_dict[m] = np.nan
    #
    #         try:
    #             meta_dict["taxonomy_name"] = s.Description.Organism.attrs[
    #                 "taxonomy_name"]
    #         except:
    #             meta_dict["taxonomy_name"] = np.nan
    #
    #         try:
    #             meta_dict["title"] = s.Description.Title.string
    #         except:
    #             meta_dict["title"] = np.nan
    #
    #         run.metadata = meta_dict
    #
    #     return None
    #
    # def to_ldict(self):
    #     ldict = []
    #     for _, run in self:
    #         run_dict = {
    #             "SRX": run.srx,
    #             "SRR": run.srr,
    #             "BioSample": run.BioSample
    #         }
    #         run_dict.update(run.metadata)
    #
    #         ldict.append(run_dict)
    #
    #     return ldict
    #
    # def to_pandas(self):
    #     return pd.DataFrame(self.to_ldict())[[
    #         "SRX",
    #         "SRR",
    #         "BioSample",
    #         "taxonomy_name",
    #         "strain",
    #         "title",
    #         "host",
    #         "isolation_source",
    #         "geo_loc_name",
    #         "sample_name"
    #     ]]

    # def show_prefetch(self):
    #     for srx, _ in self:
    #         print(f"prefetch -O . {srx}")

    # def show_rename_dump_strain(self):
    #     for _, run in self:
    #         print(f"mv {run.srr}_1.fastq {run.metadata['strain']}_1.fastq")
    #         print(f"mv {run.srr}_2.fastq {run.metadata['strain']}_2.fastq")
    #
    # def show_rename_dump_full(self):
    #     for _, run in self:
    #         print(f"mv {run.srr}_1.fastq {run.metadata['taxonomy_name'].replace(' ','_')}_{run.metadata['strain']}_1.fastq")
    #         print(f"mv {run.srr}_2.fastq {run.metadata['taxonomy_name'].replace(' ','_')}_{run.metadata['strain']}_2.fastq")


class NcbiBioSample:
    def __init__(self, accession):
        self.accession = accession
        handle = Entrez.efetch(db="biosample", id=self.accession, maxret=1)
        ret = BS(handle.read(), "xml")
        handle.close()

        s = ret.select(f"BioSample[accession={self.accession}]")
        assert len(s) == 1
        self._soup = s[0]

    def find(self, term):
        s = self._soup.find(term)
        if s is None:
            raise ValueError(f"Term not found: {term}")
        else:
            return s.string

    def get_attribute(self, harmonized_name: str):
        try:
            return self._soup.select(f'Attribute[harmonized_name="{harmonized_name}"]')[
                0
            ].string
        except IndexError:
            return np.nan


class NcbiSraExperimentPackage(object):

    def __init__(self, srx, soup):
        self.srx = srx
        self._soup = soup
        self.LIBRARY_SOURCE = self._soup.LIBRARY_SOURCE.string
        self.srr = self._get_srr()
        self.biosample = self._get_biosample()

    def _get_biosample(
        self, term="EXTERNAL_ID", attr="namespace"
    ) -> Union[NcbiBioSample, None]:

        for tag in self._soup.select(term):
            namespace_lower = tag.attrs[attr].lower()
            if namespace_lower == "biosample":
                return NcbiBioSample(accession=tag.string)

        logger.info(f"[WARNING]: No BioSample associated with SRA record: {self.srx}")
        return None

    def _get_srr(self):
        runs = self._soup.find_all("RUN")
        if len(runs) > 1:
            return [run.attrs["accession"] for run in runs]
        else:
            return [runs[0].attrs["accession"]]


# %%
class FasterqDumper:
    def __init__(
        self,
        pkgset: NcbiSraExperimentPackage,
        outdir: Path = Path("./"),
        compress: bool = True,
        compress_ncores: int = 8,
    ):
        self.pkgset: NcbiSraExperimentPackgeSet = pkgset
        self.outdir: Path = Path(outdir)
        self.sra_paths: Union[dict[str, Path], None] = None
        self.runs: Union[list[str], None] = None
        # self.sra_stats = None
        self.fastq_paths: Union[dict[str, Path], None] = None
        self.renamed_fastq_paths: Union[dict[str, Path], None] = None
        self.compress: bool = compress
        self.compress_ncores: int = compress_ncores

    def _prefetch(self):
        self.sra_paths = {}
        self.runs = []
        for srx, pkg in self.pkgset:
            cmd = ["prefetch", "-O", str(self.outdir), srx]

            logger.info(" ".join(cmd))
            pmbiproc.run_and_log(cmd, logger)

            for srr in pkg.srr:
                sra_path = self.outdir.joinpath(srr, f"{srr}.sra")
                self.runs.append(srr)
                self.sra_paths[srr] = sra_path

    def _dump(self):
        if self.runs is None:
            raise ValueError("No runs in dumper. Run _prefetch() first")
        self.fastq_paths = {}
        for run, sra_path in self.sra_paths.items():
            dump_cmd = [
                "fasterq-dump",
                "--progress",
                "-O",
                str(self.outdir),
                str(sra_path),
            ]
            logger.info(" ".join(dump_cmd))
            pmbiproc.run_and_log(dump_cmd, logger)

            self.fastq_paths[run] = []
            fq_paths = [
                f
                for f in self.outdir.iterdir()
                if re.match(f"^{run}(_[0-9])*[.]fastq$", f.name) is not None
            ]
            assert (
                len(fq_paths) != 0
                ), f"{run}: fastq files aren't where they belong - something went wrong"
            for fq_path in fq_paths:
                if self.compress:
                    compr_cmd = [
                        "pigz",
                        "--processes",
                        str(self.compress_ncores),
                        "--force",
                        str(fq_path),
                    ]
                    logger.info(
                        f"Compressing {fq_path.name} with pigz using {self.compress_ncores} processes."
                    )
                    pmbiproc.run_and_log(compr_cmd, logger)
                    self.fastq_paths[run].append(Path(f"{fq_path}.gz"))
                else:
                    self.fastq_paths[run].append(fq_path)

    def dump(self):
        self._prefetch()
        self._dump()

    def rename_fastq(
        self,
        select_fn: Callable,
        srr_rename_fn: Callable,
        base_rename_fn: Callable,
        dryrun=False,
    ):
        self.renamed_fastq_paths = {}
        for srx, pkg in self.pkgset:
            select = select_fn(pkg)
            rename = srr_rename_fn(select)
            for run in pkg.srr:
                if not dryrun:
                    self.renamed_fastq_paths[run] = []
                for fq in self.fastq_paths[run]:
                    base = base_rename_fn(re.sub(run, rename, fq.name))
                    parent = fq.parent
                    dst = parent.joinpath(base)
                    logger.info(f"Renaming fastq: {fq.name} --> {base}")
                    if not dryrun:
                        os.rename(src=fq, dst=dst)
                        self.renamed_fastq_paths[run].append(dst)
                    else:
                        logger.info("But did nothing because `dryrun` == True")


if __name__ == "__main__":
    pass
