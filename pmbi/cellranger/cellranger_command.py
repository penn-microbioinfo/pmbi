from collections import namedtuple
import tomllib
import itertools
import os
import argparse
import sys
import re
import csv
import logging
import pandas as pd
import pmbi
from pmbi.s3.lib import object_key_list, split_s3_uri
''' for EMTAB sheets
coldict = {
        9: "protocol",
        10: "filename",
        11: "read_index",
        14: "awsdir"
        }
'''
#DelimRow = namedtuple("DelimRow", ["protocol", "filename", "read_index", "awsdir"])

CRCONFIG = tomllib.load(open(os.path.join(os.path.dirname(pmbi.__file__), "cellranger", "config.toml"), "rb"))

multiome_sample_p = re.compile("([0-9]+)[_]([a-zA-Z-]+)[_]([0-9]+[a-z]*)[_]")
super_sample_p = re.compile("^([0-9]+)[_]([0-9]+[a-z]*)")
def get_patient_id(super_sample):
    return re.search(super_sample_p, super_sample).group(1)

def groupby_patient_id(super_samples: list):
    grp = {}
    for ss in super_samples:
        pid = get_patient_id(ss)
        if pid in grp:
            grp[pid].append(ss)
        else:
            grp[pid] = [ss]
    return grp
class CellrangerAggrCommand(object):
    def __init__(self, pid, mode = "count", origin_mod = "coculture"):
        self.pid = pid
        self.csv = None
        self.to_fetch = None
        self.to_fetch_bases = None
        self.origin_mod = origin_mod
        if mode in ["count", "vdj", "multi"]:
            self.mode = mode
        else:
            raise ValueError(f"Invalid CellrangerAggrCommand mode: {mode}")
    def generate_csv(self):
        csv = None
        if self.mode == "count":
            csv = "sample_id,molecule_h5\n"
            for k in self.to_fetch_bases:
                ss=re.search(super_sample_p, k).group(0)
                csv += f"{ss},../multi/{ss}/outs/per_sample_outs/{ss}/count/sample_molecule_info.h5\n"
        elif self.mode == "multi":
            csv = "sample_id,sample_outs,donor,origin\n"
            for k in self.to_fetch_bases:
                ss=re.search(super_sample_p, k).group(0)
                csv += f"{ss},../multi/{ss}/outs/per_sample_outs/{ss},{self.pid},{self.pid}_{self.origin_mod}\n"
        else:
            csv = "sample_id,vdj_contig_info,donor,origin\n"
            for k in self.to_fetch_bases:
                ss=re.search(super_sample_p, k).group(0)
                csv += f"{ss},../multi/{ss}/outs/per_sample_outs/{ss}/vdj_t/,{self.pid}_{self.origin_mod}n"
        self.csv = csv 
        return self.csv
        
    def from_s3_archive_keys(pid, archive_keys, mode = "count"):
        archive_bases = [os.path.basename(x) for x in archive_keys]
        cmd = CellrangerAggrCommand(pid, mode=mode)
        cmd.to_fetch = archive_keys
        cmd.to_fetch_bases = archive_bases
        cmd.generate_csv()
        return cmd
    def awsdirs(self):
        return set([os.path.dirname(p) for p in self.to_fetch])
class CellrangerAggrGroups(object):
    def __init__(self):
        self.groups = {}
    def insert(self, pid, value):
        if pid in self.groups:
            self.groups[pid].append(value)
        else:
            self.groups[pid] = [value]

def check_multiome_experiments(experiments):
    allmodals = set(itertools.chain(*[[xx.type for xx in x.modalities()] for x in experiments]))
    nmodals = [len(x.modalities()) for x in experiments]
    max_modals = len(allmodals)
    if not all([x == max_modals for x in nmodals]):
        for exp in experiments:
            if len(exp.modalities()) != max_modals:
                maybe_missing_modals = [x for x in experiments if x.identifier == re.sub("[a-z]{1}$", "", exp.identifier)]
                if len(maybe_missing_modals) > 0:
                    for mm in maybe_missing_modals:
                        for m in mm.modalities():
                            if exp.try_add_modality(m):
                                logging.warning(f"Added missing modality to experiment {exp.identifier}: ({m.super_sample}, {m.type})")
                                mm.redundant=True
                    if len(exp.modalities()) != max_modals:
                        logging.warning(f"Unable to find missing modalities for {exp.identifier}. Continuing with only: {[m.type for m in exp.modalities()]}")
    finalized_experiments = []
    for exp in experiments:
        if not exp.redundant:
            finalized_experiments.append(exp)
        else:
            logging.warning(f"Skipping redundant experiment for {exp.identifier} with modalities: {[x.type for x in exp.modalities()]}")
    return finalized_experiments

def parse_keys_multiome(keys):
    super_samples = {}
    if not all([x == [os.path.dirname(y) for y in keys][0] for x in [os.path.dirname(x) for x in keys]]):
        logging.warning(f"awsdirs vary: {set([os.path.dirname(w) for w in keys])}")
    for key in keys:
        keyspl = os.path.split(key)
        awsdir = keyspl[0]
        key = keyspl[1]
        s = re.search(multiome_sample_p, key)
        if s is not None:
            sample = re.sub("[_]$", "", s.group(0))
            super_sample = f"{s.group(1)}_{s.group(3)}"
            modal = s.group(2)
            if not super_sample in super_samples:
                super_samples[super_sample] = [CellrangerModality(super_sample, sample, modal, awsdir)]
            else:
                super_samples[super_sample].append(CellrangerModality(super_sample, sample, modal, awsdir))
        else:
            raise ValueError(f"No matching sample found in: {key}")
    experiments = []
    for k,v in super_samples.items():
        experiments.append(CellrangerMultiomeExperiment.from_modalities(k,v))
    experiments = check_multiome_experiments(experiments)
    return experiments


class CellrangerModality(object):
    def __init__(self, super_sample, sample, modality, awsdir):
        self.super_sample = super_sample
        self.sample = sample
        self.type = modality
        self.awsdir = awsdir

class CellrangerMultiomeExperiment(object):
    def __init__(self, ident):
        self.identifier = ident
        self.redundant = False
        pass

    def modalities(self):
        return [v for v in self.__dict__.values() if isinstance(v, CellrangerModality)]

    def try_add_modality(self, modality: CellrangerModality):
        if modality.type not in self.__dict__:
            setattr(self, modality.type, modality)
            return True
        else:
            return False

    def from_modalities(ident, modalities):
        exp = CellrangerMultiomeExperiment(ident)
        for m in modalities:
            exp.try_add_modality(m)
        return exp

class CellrangerCommand(object):
    def __init__(self, protocol, reference = None, **kwargs):

        self.command = "cellranger"
        self.subcommand = None
        self.protocol = protocol
        self.command_line_args = dict(kwargs)
        self.sample = self.command_line_args["sample"]
        self.archive_name = None
        self.important_outs = []
        self.path_to_outs = None

        if reference is None:
            self.command_line_args["reference"] = CRCONFIG["references"][protocol]
        else: 
            self.command_line_args["reference"] = reference

    def compress_outs(self):
        if self.archive_name is None or self.path_to_outs is None:
            raise ValueError("Cannot compress outs whose locations are not known.")
        cmds = [f"tar cvfz {self.archive_name} {self.path_to_outs} && \\"]
        if len(self.important_outs) > 0:
            cmds.append("mkdir important_outs")
            for imp in self.important_outs:
                cmds.append(f"mv {imp} important_outs/.")
            cmds.append(f"tar cvzf {self.archive_name.replace('tar.gz', 'importantOuts.tar.gz')} important_outs/")
        return '\n'.join(cmds)

    def get_reference_path(self):
        return self.protocol_to_ref[self.protocol]

    def command_line_args_repr(self, exclude: list):
        cli_args = [f"--{k.replace('_', '-')} {v}" for k,v in self.command_line_args.items() if k not in exclude]
        if self.subcommand is not None:
            return " ".join([self.command, self.subcommand] + cli_args)
        else:
            return " ".join([self.command, self.subcommand] + cli_args)


    def __str__(self):
        return self.command_line_args_repr(exclude = []) 

class CellrangerCountRNA(CellrangerCommand):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.subcommand = "count"
        self.command_line_args["include_introns"] = str(self.should_include_introns()).lower()
        self.command_line_args["transcriptome"] = self.get_reference_path() 
        self.archive_name = f"{self.command_line_args['sample']}_filtered_feature_bc_matrix.tar"
        self.path_to_outs = f"{self.command_line_args['sample']}/outs/filtered_feature_bc_matrix"
    
    # Out of date because cellranger now ALWAYS includes introns by default
    #def should_include_introns(self):
    #    if self.protocol == "snRNA-seq":
    #        return True
    #    else:
    #        return False

class CellrangerCountADT(CellrangerCommand):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.subcommand = "count"
        self.command_line_args["feature-ref"] = self.get_reference_path()
        self.command_line_args["transcriptome"] = self.protocol_to_ref["scRNA-seq"]
        self.archive_name = f"{self.command_line_args['sample']}_filtered_feature_bc_matrix.tar"
        self.path_to_outs = f"{self.command_line_args['sample']}/outs/filtered_feature_bc_matrix"
        
        self.command_line_args["libraries"] = "libraries.csv"

    def bash_make_libraries_csv(self):
        return f"""echo fastqs,sample,library_type > libraries.csv 
        echo $(readlink -f {self.command_line_args['fastqs']}),{self.command_line_args['sample']},Antibody Capture >> libraries.csv
        """
    def __str__(self):
         return "\n".join([self.bash_make_libraries_csv(), self.command_line_args_repr(exclude = ["fastqs", "sample"])])

class CellrangerVDJ(CellrangerCommand):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.subcommand = "vdj"
        self.command_line_args["reference"] = self.get_reference_path()
        self.archive_name = f"{self.command_line_args['sample']}_vdj_outs.tar"
        self.path_to_outs = f"{self.command_line_args['sample']}/outs"

class CellrangerAtac(CellrangerCommand):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.command = "cellranger-atac"
        self.subcommand = "count"
        self.archive_name = f"{self.sample}.{self.protocol}.tar.gz"
        self.path_to_outs = f"{self.sample}/outs"
        self.important_outs = [
                f"{self.path_to_outs}/filtered_peak_bc_matrix/",
                f"{self.path_to_outs}/filtered_peak_bc_matrix.h5",
                f"{self.path_to_outs}/filtered_tf_bc_matrix/",
                f"{self.path_to_outs}/filtered_tf_bc_matrix.h5",
                f"{self.path_to_outs}/raw_peak_bc_matrix/",
                f"{self.path_to_outs}/raw_peak_bc_matrix.h5",
                f"{self.path_to_outs}/web_summary.html",
                f"{self.path_to_outs}/peak_annotation.tsv",
                f"{self.path_to_outs}/peak_motif_mapping.bed",
                f"{self.path_to_outs}/peaks.bed",
                f"{self.path_to_outs}/analysis/umap",
                f"{self.path_to_outs}/cloupe.cloupe"
                ]

class CellrangerMulti(object):
    def __init__(self, fastqs, config, experiment, output_key_prefix):
        self.subcommand = "multi"
        self.config = config
        self.experiment = experiment
        self.commands = []
        self.output_key_prefix = output_key_prefix
        self.path_to_outs = f"{self.experiment.identifier}/outs"
        self.archive_name = f"{self.experiment.identifier}.multi.outs.tar.gz"
        self.important_archive_name = f"{self.experiment.identifier}.multi.importantOuts.tar.gz"
     
    def print_s3_download_cmds(self): 
        cmds=[]
        for m in self.experiment.modalities():
            cmds.append(f"python /shared-ebs/microbioinfo-aws/scripts/s3-download-multi.py --bucket microbioinfo-storage --prefix {os.path.join(m.awsdir, m.sample)} --chunksize 1.024e8 --pattern '[_]R[0-9]+[_][0-9]+[.]fastq[.]gz$'")
        return '\n'.join(cmds)
    
    def print_write_config_cmd(self):
        return f"cat << EOF > config.csv\n{self.config}\nEOF"

    def print_cellranger_cmd(self):
        return f"cellranger {self.subcommand} --id {self.experiment.identifier} --csv config.csv --localcores 16 --localmem 58"

    def print_compress_outs_cmd(self):
        if self.archive_name is None or self.output_key_prefix is None:
            raise ValueError("Cannot compress outs whose locations are not known.")
        return f"tar cvfz {self.archive_name} {self.path_to_outs} && rm -r {self.path_to_outs} && tar cvfz {self.important_archive_name} important_outs && rm -r important_outs"
    
    def print_pull_important_outs_cmd(self):
        return f"""
                mkdir important_outs && \\
                cp -r {self.path_to_outs}/per_sample_outs/{self.experiment.identifier}/vdj_t important_outs/. && \\
                cp {self.path_to_outs}/per_sample_outs/{self.experiment.identifier}/count/sample_filtered_feature_bc_matrix.h5 important_outs/. && \\
                cp -r {self.path_to_outs}/per_sample_outs/{self.experiment.identifier}/web_summary.html important_outs/. && \\
                cp {self.path_to_outs}/multi/count/raw_feature_bc_matrix.h5 important_outs/. && \\
                cp {self.path_to_outs}/multi/count/raw_molecule_info.h5 important_outs/.
                """

    # After cellranger runs, move the outs directory out of rundir and delete rundir (for disk space)
    def print_move_outs_delete_wd(self):
        ret = f"if [ -d outs ]; then rm -r outs; fi; if [ -d important_outs ]; then rm -r important_outs; fi; mv {self.path_to_outs} .\nrm -r {self.experiment.identifier}"
        self.path_to_outs = "outs"
        return ret
    def print_s3_upload_cmd(self):
        out_up = f"python /shared-ebs/pmbi/pmbi/s3/s3-multipart-upload.py --bucket microbioinfo-storage --key {os.path.join(self.output_key_prefix, self.archive_name)} --partsize 1000000000 --nproc 8 --largefile {self.archive_name}"
        important_out_up = f"python /shared-ebs/microbioinfo-aws/scripts/s3-multipart-upload.py --bucket microbioinfo-storage --key {os.path.join(self.output_key_prefix, self.important_archive_name)} --partsize 1000000000 --nproc 8 --largefile {self.important_archive_name}"
        return '\n'.join([out_up, important_out_up])

class CellrangerMultiConfigSection(object):
    def __init__(self, name):
        self.name = name
        self.has_column_names = False
        self.colnames = []
        self.rows = []
        self.kvpairs = []
        self.type = None
    def set_colnames(self, colnames: list):
        self.has_column_names = True
        self.colnames = colnames

    def _framify(self):
        if self.type == "rows":
            colnames_repr = ','.join(self.colnames)
            rows_repr = '\n'.join([','.join(x) for x in self.rows])
            return f"{colnames_repr}\n{rows_repr}"
        elif self.type == "kvpairs":
            kvpairs_repr = '\n'.join([','.join(x) for x in self.kvpairs])
            return kvpairs_repr

    def _assess(self):
        if len(self.rows) > 0 and len(self.kvpairs) > 0:
            raise ValueError(f"Section {self.name} cannot be composed of rows and key-value pairs.")
        else:
            if len(self.rows) == 0 and len(self.kvpairs) == 0:
                raise ValueError(f"Section {self.name} has to have one of: rows or key-value pairs.")
            else:
                if len(self.rows) > 0:
                    self.type = "rows"
                elif len(self.kvpairs) > 0:
                    self.type = "kvpairs"
                else:
                    raise NotImplementedError
    def add_kvpair(self, kvpair: tuple):
        self.kvpairs.append(kvpair)

    def add_row(self, rowdict: dict):
        if len(self.colnames) == 0:
            raise ValueError("Define column names before adding rows")
        else:
            row = []
            for c in self.colnames:
                if c in rowdict:
                    row.append(rowdict.pop(c))
                else:
                    raise ValueError(f"'{c}' does not occur in the row: {rowdict}")
            if len(rowdict) != 0:
                raise ValueError(f"Some columns did not appear in colnames for row: {rowdict}")
            self.rows.append(row)

    def write(self):
        self._assess()
        return self._framify()

class CellrangerMultiConfig(object):
    def __init__(self):
        pass

    def modality_to_section_header():
        return {
                "RNA": "gene-expression",
                "VDJ": "vdj",
                "ADT": "feature"
                }
    
    def modality_to_feature_type():
        return {
            "RNA": "Gene Expression",
            "VDJ": "VDJ",
            "ADT": "Antibody Capture" 
            }

    def add_section(self, name: str, section: CellrangerMultiConfigSection):
        setattr(self, name, section) 

    def from_multi_experiment(path_to_fastqs, experiment: CellrangerMultiomeExperiment, refs: dict):
        config = CellrangerMultiConfig()
        m2sh = CellrangerMultiConfig.modality_to_section_header()
        m2ft = CellrangerMultiConfig.modality_to_feature_type()
        libraries = CellrangerMultiConfigSection("libraries")
        libraries.set_colnames(["fastq_id", "fastqs", "feature_types"])
        for modality in experiment.modalities():
            section = CellrangerMultiConfigSection(modality.type)
            section.add_kvpair(("reference", refs[modality.type]))
            config.add_section(m2sh[modality.type], section)
            libraries.add_row({
                "fastq_id": modality.sample,
                "fastqs": path_to_fastqs,
                "feature_types": m2ft[modality.type]
                })
        config.add_section("libraries", libraries)
        return config

    def __str__(self):
        sections = []
        for k,v in self.__dict__.items():
            section_header = f"[{k}]"
            section = v.write()
            sections.append(f"{section_header}\n{section}\n") 
        return '\n'.join(sections)

    def write(self, fname):
        with open("fname", 'w') as config_out:
            config_out.write(self.__str__())
 
# This is stupid - just use pandas
#def cellranger_table_column_dict(protocol: int = 0, filename: int = 1, read_index: int = 2, awsdir: int = 3):
#    coldict = {
#            protocol: "protocol",
#            filename: "filename",
#            read_index: "read_index",
#            awsdir: "awsdir"
#            }
#    return coldict

def cellranger_cmd(protocol, protocol_to_cr_command_generator):
    return protocol_to_cr_command_generator[protocol]

def cellranger_count_cmd(fastq_dir, ref, ident, sample, include_introns = False, localcores = 1, localmem = 8, chemistry = "auto"):
    return f"cellranger \
            count \
            --fastqs {fastq_dir} \
            --transcriptome {ref} \
            --id {ident} \
            --sample {sample} \
            --include-introns {str(include_introns).lower()} \
            --localcores {localcores} \
            --localmem {localmem} \
            --chemistry {chemistry}"

def cellranger_vdj_cmd(fastq_dir, ref, ident, sample, localcores = 1, localmem = 8):
    return f"cellranger \
            count \
            --fastqs {fastq_dir} \
            --reference {ref} \
            --id {ident} \
            --sample {sample} \
            --localcores {localcores} \
            --localmem {localmem}"


def cellranger_commands(table, s3_output_prefix = "cellranger_matrices"):
    cmds = {}
    crtable = pd.read_csv(table, sep = "\t")
    for row in crtable.itertuples():
        this_sample = []
        if row.read_number in ["R1", "read1"]:
            sampleid = re.sub("_S[0-9]+[_]L[0-9]+[_][IR][0-9]+[_][0-9]+[.]fastq[.]gz", "", row.read_filename)
            bucket,r1prefix = split_s3_uri(uri = row.read_s3_uri)
            read_prefix = os.path.join(os.path.dirname(r1prefix), sampleid)
            this_sample.append(f"python /shared-ebs/pmbi/pmbi/s3/s3-download-multi.py --bucket {bucket} --prefix {read_prefix} --chunksize 1.024e8 --pattern '[_]R[0-9]+[_][0-9]+[.]fastq[.]gz$' && \\")
            this_cellranger_cmd = cellranger_cmd(row.modality, modality_to_class)(
                    protocol = row.modality,
                    fastqs = "fastqs", 
                    id = sampleid, 
                    sample= sampleid, 
                    localcores = 8, 
                    localmem=58)
            this_sample.append(this_cellranger_cmd.__str__() + " && \\")

            this_sample.append(this_cellranger_cmd.compress_outs())
            output_uri = f"s3://{bucket}/{s3_output_prefix}/{this_cellranger_cmd.archive_name}"
            this_sample.append(f"python /shared-ebs/pmbi/pmbi/s3/s3-multipart-upload.py --uri {output_uri} --partsize 100000000 --nproc 8 --skip_verify_etags --largefile {this_cellranger_cmd.archive_name}")
            this_sample.append(f'python /shared-ebs/pmbi/pmbi/s3/s3-multipart-upload.py --uri {output_uri.replace(".tar.gz", ".importantOuts.tar.gz")} --partsize 100000000 --nproc 8 --skip_verify_etags --largefile {this_cellranger_cmd.archive_name.replace('.tar.gz', '.importantOuts.tar.gz')}')
            if sampleid in cmds:
                assert this_sample == cmds[sampleid]
            else:
                cmds[sampleid] = this_sample

    return cmds

modality_to_class = {
        "RNA": CellrangerCountRNA,
        "VDJ": CellrangerVDJ,
        "ADT": CellrangerCountADT,
        "ATAC": CellrangerAtac,
        }

emtab_protocol_to_cr_command_generator = {
        "scRNA-seq": CellrangerCountRNA,
        "scVDJ-seq": CellrangerVDJ,
        "scADT-seq": CellrangerCountADT
        }
def newDelimRow(row, coldict, delim = '\t'):
    spl = row.split(delim)
    cols = []
    for i in coldict:
        cols.append(spl[i])
    return DelimRow._make(tuple(cols))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required = True)
    parser.add_argument("--r1_filename", required = False)
    args = parser.parse_args()


