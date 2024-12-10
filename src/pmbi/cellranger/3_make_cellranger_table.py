# %% Imports
import argparse
import copy
import io
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import boto3
import munch
import pandas as pd
from joblib import Parallel, delayed

import pmbi.config
from pmbi.s3.lib import object_key_list


# %%






# class LocalHandler(Handler):
#     def __init__(
#         self, path, working_dir, pattern="[_]R[0-9][_].+[.]fastq[.]gz$", backend=None
#     ):
#         super().__init__(path, working_dir, pattern, backend)
#
#     def _pull_files(self):
#         for f in self._files():
#             os.symlink(f, os.path.join(self.working_dir, "fastq", os.path.basename(f)))
#
#     def _run(self, cmd):
#         pass
#
# def path_to_dict(path: Path):
#     if len(path.parts) > 2:
#         return {path.parts[0]: path_to_dict(Path(*path.parts[1::]))}
#     elif len(path.parts) == 2:
#         return {path.parts[0]: path.parts[1]}
#     else:
#         raise ValueError("Should not get a path.parts of len 1")


# def validate_path_dict(pathdict: dict):
#     if len(pathdict) > 1:
#         raise ValueError("Len of pathdict level > 1")
#     for k,v in pathdict.items():
#         if isinstance(v, dict):
#             validate_path_dict(v)
#         elif isinstance(v, str):
#             print(v)
#             pass
#         else:
#             raise ValueError("Invalid type in pathdict")
#
# def pathdict_append(pathdict: dict, p2: dict):
#     if len(pathdict) > 1:
#         raise ValueError("Len of pathdict level > 1")
#     for k,v in pathdict.items():
#         if isinstance(v, dict):
#             pathdict_append(v, p2)
#         elif isinstance(v, str):
#             return {v: p2}
#         else:
#             raise ValueError("Invalid type in pathdict")


# %% CHUNK: Command def {{{
class Command(object):
    def __init__(self, handler, backend, sample_config, config):
        self.handler = handler
        self.backend = backend
        self.sample_config = sample_config
        self.config = config
        self.inner = None
        self.subdirs = [
            Path("fastq"),
            Path("outs_archive"),
            Path("cellranger_wd"),
        ]
        self.wd = Path(self.config.run.wd).joinpath(self.sample_config.sample)
        self.wd_tree = []
        # if self.wd.exists() and self.wd.is_dir():
        #     if self.config.command_line.force:
        #         print("Deleting existing run dir because of force.")
        #     else:
        #         raise IOError("Run dir already exists.")
        self.wd_tree = [self.wd.joinpath(p) for p in self.subdirs]
        self._make_wd_tree()
        self.sample_config.write(outdir=self.wd)
        self._gen()
        print(self.inner)

    def _gen(self):
        if self.backend == "local":
            self.inner = [
                [
                    "cellranger",
                    "multi",
                    "--id",
                    self.sample_config.sample,
                    "--csv",
                    self.wd.joinpath(self.sample_config.relpath),
                    "--localcores",
                    self.config.commands.cores,
                    "--localmem",
                    self.config.commands.memory,
                ],
            ]

    def _make_wd_tree(self):
        for path in self.wd_tree:
            print(path)
            os.makedirs(path, exist_ok=False)

    def _run(self):
        os.chdir(self.wd)
        for cmd in self.inner:
            cmd = [str(e) for e in cmd]
            print(cmd)
            print(" ".join(cmd))
            p = subprocess.Popen(
                cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            out, err = p.communicate()
        return p.returncode, out, err


# }}}
# %% CHUNK: SampleConfig def {{{
class SampleConfig(object):
    def __init__(self, sample, modalities, config):
        self.sample = sample
        self.modalities = modalities
        self.relpath = f"./{self.sample}_config.csv"

    def _modality_configs(self):
        seen = []
        sections = ""
        for m in self.modalities.itertuples():
            if m.toml_header not in seen:
                sections += f"""[{m.toml_header}]
                reference,{getattr(config.references, m.modality)}

                """
                seen.append(m.toml_header)
        return sections

    def _libraries(self):
        columns = ["fastq_id", "fastqs", "feature_types"]
        section = self.modalities[columns].to_csv(sep=",", index=False)
        return f"[libraries]\n{section}"

    def write(self, outdir):
        config = ""
        config += self._modality_configs()
        config += self._libraries()
        with open(os.path.join(outdir, self.relpath), "w") as sc:
            sc.write("\n".join([l.strip() for l in config.splitlines()]))


# }}}


# %% CHUNK: Run def {{{
class Run(object):
    def __init__(self, handler: Handler, config: munch.Munch):
        self.handler = handler
        self.config = config
        self.samples = {}
        self.wd = Path(self.config.run.wd)
        self.modality_to_cr_feature_type = {
            "RNA": "Gene Expression",
            "ADT": "Antibody Capture",
            "VDJ-T": "VDJ-T",
            "VDJ-B": "VDJ-B",
        }
        self.modality_to_cr_toml_header = {
            "RNA": "gene-expression",
            "ADT": "feature",
            "VDJ-T": "vdj",
            "VDJ-B": "vdj",
        }
        if self.wd.exists() and self.wd.is_dir():
            if self.config.command_line.force:
                print("Deleting existing run dir because of force.")
                shutil.rmtree(self.wd)
            else:
                raise IOError("Run dir already exists.")

        os.makedirs(self.wd, exist_ok=False)
        self.handler.table.to_csv(
            self.wd.joinpath("Handler_table.tsv"), sep="\t", index=False
        )
        self._populate()

    def _populate(self):
        cr_config = (
            self.handler.table[
                ["sample_rep", "sample", "modality", "technical_rep", "read_dir"]
            ]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        cr_config["feature_types"] = pd.Series(
            [self.modality_to_cr_feature_type[m] for m in cr_config["modality"].values]
        )
        cr_config["toml_header"] = pd.Series(
            [self.modality_to_cr_toml_header[m] for m in cr_config["modality"].values]
        )
        cr_config["fastq_id"] = pd.Series(
            [
                f"{row.sample}_{row.modality}_{row.technical_rep}"
                for row in cr_config[
                    ["sample", "modality", "technical_rep"]
                ].itertuples()
            ]
        )
        # cr_config["sample_rep"] = cr_config["sample"] + "_" + cr_config["technical_rep"]
        cr_config = cr_config[
            [
                "sample_rep",
                "modality",
                "fastq_id",
                "read_dir",
                "feature_types",
                "toml_header",
            ]
        ]
        cr_config = cr_config.rename(columns={"read_dir": "fastqs"})
        for sample_rep, modalities in cr_config.groupby(["sample_rep"]):
            sample_rep = sample_rep[0]
            sc = SampleConfig(sample=sample_rep, modalities=modalities, config=config)
            if sample_rep in self.samples:
                raise KeyError(f"Duplicate sample_rep names: {sample_rep}")
            sub_handler = self.handler.subset(lambda r: r["sample_rep"] == sample_rep)
            # sub_handler.working_dir = Path(sub_handler.working_dir).joinpath(sample_rep)
            self.samples[sample_rep] = Command(
                handler=sub_handler,
                backend=self.config.command_line.backend,
                sample_config=sc,
                config=config,
            )

    def _run(self):
        def _compute(cmd):
            return cmd._run()

        print(self.samples.values())
        ret = Parallel(n_jobs=self.config.run.n_jobs)(
            delayed(_compute)(cmd) for cmd in self.samples.values()
        )
        print(ret)


# }}}


if __name__ == "__main__":

    # %% CHUNK: Argparser {{{
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--path",
        required=True,
        type=str,
        help="Local path or s3 URI to pointing to where the reads are located.",
    )
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        type=str,
        help="Path to pmbi cellranger config file.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        help="Force creation of output directories.",
    )
    parser.add_argument(
        "-b",
        "--backend",
        required=False,
        type=str,
        choices=["s3", "local"],
        help="Explicit backend to use. Will attempt to infer from path.",
    )
    parser.add_argument(
        "-i",
        "--include",
        required=False,
        default="[.]fastq[.]gz$",
        type=str,
        help="Python-style regex that matches some of all the filenames to include.",
    )
    args = parser.parse_args()
    # }}}

    config = pmbi.config.import_config(args.config)
    for k, v in args.__dict__.items():
        setattr(config.command_line, k, v)

    handler = LocalHandler(
        path=config.command_line.path,
        # working_dir=config.run.wd,
        pattern=config.filename_patterns.include,
    )

    run = Run(handler=handler, config=config)
    print(run.handler.table)
    for samp, cmd in run.samples.items():
        cmd.handler._pull_files()
    run._run()
