# %% Imports
import argparse
import io
import os
import re
import sys

import boto3
import munch
import pandas as pd

import pmbi.config
from pmbi.s3.lib import object_key_list


# %%
def get_substring(string: str, pattern: str) -> str:
    s = re.search(pattern, os.path.basename(string))
    if s is not None:
        if len(s.groups()) > 1:
            raise ValueError(
                f"More than one matching substring for pattern `{pattern}` in string: {string}"
            )
        else:
            return s.group(1)
    else:
        raise ValueError(
            f"No match found for pattern `{pattern}` in string: `{string}`"
        )


accepted_modalities = {
    "RNA": "RNA",
    "ADT": "ADT",
    "HTO": "HTO",
    "VDJ": "VDJ",
    "VDJ-T": "VDJ-T",
    "VDJ-B": "VDJ-B",
}


# %% CHUNK: Handler def {{{
class Handler(object):
    def __init__(self, path, working_dir, pattern="[_]R[0-9][_].+[.]fastq[.]gz$", backend=None):
        self.filelist = []
        self.s3_pattern = "^s3[:][/]{2}"
        self.working_dir = working_dir
        self.table: pd.DataFrame = pd.DataFrame()
        if backend is None:
            if re.search(self.s3_pattern, path):
                self.backend = "s3"
                self.filelist = object_key_list(path)
            else:
                self.backend = "local"
                self.filelist = [os.path.join(path, f) for f in os.listdir(path)]
        else:
            raise NotImplementedError
        self._filter_filelist(pattern)
        self._make_table()

    def _filter_filelist(self, pattern):
        filtered = []
        for f in self.filelist:
            fb = os.path.basename(f)
            if re.search(pattern, fb) is not None:
                filtered.append(f)
        self.filelist = filtered

    def _files(self):
        for f in self.filelist:
            yield f

    def _make_table(self):
        ldict = []
        for f in self._files():
            fb = os.path.basename(f)
            ldict.append(
                {
                    "sample": get_substring(
                        string=fb, pattern=config.filename_patterns.sample
                    ),
                    "modality": get_modality_from_string(
                        string=fb,
                        pattern=config.filename_patterns.modality,
                        accepted_modalities=accepted_modalities,
                    ),
                    "technical_rep": get_substring(
                        string=fb, pattern=config.filename_patterns.technical_rep
                    ),
                    "read_number": get_substring(
                        string=fb, pattern=config.filename_patterns.read_number
                    ),
                    "read_filename": fb,
                    "read_path": f,
                    "read_dir": os.path.split(f)[0],
                    "backend": self.backend,
                }
            )
        df = pd.DataFrame(ldict)
        df = df.sort_values(["read_filename", "modality", "read_number"])
        self.table = df
        df.to_csv(os.path.join(self.working_dir, "Handler_table.tsv"), sep="\t", index=False)

    def _sync_files(self):
        pass

# }}}


# %% CHUNK: This function tries to pull an accepted modality from a an input string based on an input pattern {{{
def get_modality_from_string(
    string: str,
    pattern: str,
    accepted_modalities: dict[str, str],
) -> str:
    modality = get_substring(string, pattern)
    if modality in accepted_modalities:
        return accepted_modalities[modality]
    else:
        raise ValueError(f"Unexpected modality: `{modality}`")


# }}}

# %% CHUNK: SampleConfig def {{{
class Command(object):
    def __init__(self, backend, sample_config, config):
        self.backend = backend
        self.sample_config = sample_config
        self.config = config
        self.inner = None
        self.subdirs = {
            "fastq": {
                "sub1": {
                    "sub2": {},
                    },
                },
            "outs_archive": {},
            "cellranger_wd": {},
            }
        self.wd_tree = {self.sample_config.sample: self.subdirs}
    def _gen(self):
        csv = f"./{self.sample_config.sample}_config.csv"
        self.sample_config.write(output = csv)
        if self.backend == "local":
            self.inner = ["cellranger", "multi", "--id", self.sample_config.sample, "--csv", csv]
    def _make_wd_tree(self, wd_tree, leading = []):
        for k,v in wd_tree.items():
            os.mkdir('/'.join(leading+[k]))
            self._make_wd_tree(v, leading+[k])
    def run(self):
        pass
        



# }}}
# %% CHUNK: SampleConfig def {{{
class SampleConfig(object):
    def __init__(self, sample, modalities, config):
        self.sample = sample
        self.modalities = modalities

    def _modality_configs(self):
        seen = []
        sections = ""
        for m in self.modalities.itertuples():
            if m.toml_header not in seen:
                sections += f"""[{m.toml_header}]
                reference = {getattr(config.references, m.modality)}

                """
                seen.append(m.toml_header)
        return sections

    def _libraries(self):
        columns = ["fastq_id", "fastqs", "feature_types"]
        section = self.modalities[columns].to_csv(sep=",", index=False)
        return f"[libraries]\n{section}"

    def write(self, output=None):
        config = ""
        config += self._modality_configs()
        config += self._libraries()
        with open(output, 'w') as sc:
            sc.write("\n".join([l.strip() for l in config.splitlines()]))

# }}}

# %% CHUNK: Run def {{{
class Run(object):
    def __init__(self, handler: Handler, config: munch.Munch):
        self.handler = handler
        self.config = config
        self.samples = {}
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

        self._populate()

    def _populate(self):
        cr_config = (
            self.handler.table[["sample", "modality", "technical_rep", "read_dir"]]
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
        cr_config["sample_rep"] = cr_config["sample"] + "_" + cr_config["technical_rep"]
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
        print(cr_config)
        for sample_rep, modalities in cr_config.groupby(["sample_rep"]):
            sample_rep = sample_rep[0]
            sc = SampleConfig(sample=sample_rep, modalities=modalities, config=config)
            if sample_rep in self.samples:
                raise KeyError(f"Duplicate sample_rep names: {sample_rep}")
            self.samples[sample_rep] = Command(backend = "local", sample_config = sc, config = config)
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
        "-o",
        "--output",
        default="./cellranger_table.tsv",
        type=str,
        help="Path to write cellranger table out to.",
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
    handler = Handler(path = args.path, working_dir = config.run.wd, pattern = config.filename_patterns.include)

    run = Run(handler=handler, config=config)
    print(run.handler.table)
    os.chdir(run.config.run.wd)
    for samp,cmd in run.samples.items():
        # cmd._gen()
        cmd._make_wd_tree(cmd.wd_tree)
        break
