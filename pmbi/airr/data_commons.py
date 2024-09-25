# %%
import os
import awkward as ak
import anndata
import itertools
import json
import copy
import pandas as pd
import numpy as np
import pmbi.anndata.io

from requests_html import HTMLSession

# %%
class Filter(object):
    def __init__(self, op: str, content: list|dict):
        self.op = op
        self.content = content
        self._check_types()
        if isinstance(self.content, list):
            self.content = [f.to_json() for f in self.content]
        else:
            self.content = self.content
    def _check_types(self):
        if isinstance(self.content, list):
            assert all([isinstance(f, Filter) for f in self.content]), "Content lists for Filter must be made up of other Filter instances"
        elif isinstance(self.content, dict):
            pass
        else:
            raise ValueError(f"Invalid content type: {type(self.content)}")
    def to_json(self) -> dict:
            return {
                    "op": self.op,
                    "content": self.content
                    }

# %%
class PostData(object):
    def __init__(
        self, 
        filters: list[Filter]|None = None,
        first_result: int = 0,
        size: int|None = None,
        fields: list[str]|None = None,
        format: str = "json"
    ):
        self.filters = filters
        self.first_result = first_result
        self.size = size
        self.fields = fields
        self.format = format
    def to_json(self) -> dict:
        selfcp = copy.deepcopy(self)
        if selfcp.filters is not None:
            if len(selfcp.filters) > 1:
                selfcp.filters = [f.to_json() for f in selfcp.filters]
            else:
                selfcp.filters = selfcp.filters[0].to_json()
        json = {}
        for key,value in selfcp.__dict__.items():
            if value is not None:
                json[key] = value
        return json

# %%
class Response(object):
    def __init__(self, query, info, schema, text):
        self.query = query
        self.info = info
        self.schema = schema
        self.text = text

# %%
def request(url: str, data: PostData) -> Response:
    session = HTMLSession()
    response = session.post(url, json=data.to_json(), headers = {"content-type": "applications/json"})
    session.close()
    try: 
        text = json.loads(response.text)
        schema_key = os.path.split(url)[1].capitalize() #The last piece of the url, in sentence case should be the scheme key
        if schema_key in text:
            return Response(query = data, info = text["Info"], schema = schema_key, text = text)
        else: 
            raise ValueError(f"Expected schema_key not present in response text: {schema_key}")
    except json.JSONDecodeError as e:
        print("Error parsing json...")
        print(response.text)
        raise e

# %%
def is_null_column(column: pd.Series) -> bool:
    if pd.api.types.is_object_dtype(column):
        if any([isinstance(v, str) for v in column]):
            if all([pd.isnull(v) or len(v) == 0 for v in column]):
                return True
            else:
                return False
        elif all(pd.isnull(column)):
            return True
        else:
            return False
    else:
        return False

# %%
    keep = np.where(seqdf.apply(lambda col: any([c is not None and c != "" for c in col]), axis = 0)) 
########################################
# %% Get the repertoire_ids for this project
########################################
data = PostData(
    filters=[ Filter(op = "=", content = {"field": "study.study_id", "value": "IR-Prak-000001"})],
    first_result=0,
    size=None,
    fields=["repertoire_id", "sample.sample_id"],
)
url = "https://hpap.ireceptor.org/airr/v1/repertoire"
data.to_json()
response = request(url, data)
repids = [rec["repertoire_id"] for rec in response.text["Repertoire"]]

########################################
# %% Get the sequence ids for these repertoire_ids
########################################
url = "https://hpap.ireceptor.org/airr/v1/rearrangement"
data = PostData(filters = [Filter(op = "in", content = {"field": "repertoire_id", "value": repids})], 
                first_result = 0, 
                size = None
                )
data.to_json()
response = request(url, data)
response.text

seqids = [rec["sequence_id"] for rec in response.text["Rearrangement"]]
seqids[1:10]

########################################
# %% Get the sequence data for these sequence_ids
########################################
url = "https://hpap.ireceptor.org/airr/v1/rearrangement"
all_seqid_resp = []
schema = "Rearrangement"
seqids = np.array(seqids)
for chunk in np.array_split(range(0,len(seqids)), 10):
    print(chunk)
    data = PostData(filters = [Filter(op = "in", content = {"field": "sequence_id", "value": seqids[chunk].tolist()})], 
                    first_result = 0, 
                    size = None
                    )
    all_seqid_resp.append(request(url, data).text[schema])

all_seqid_resp = list(itertools.chain.from_iterable(all_seqid_resp))
len(all_seqid_resp)

########################################
# %% Convert to DataFrame and drop columns that are all None/Null
# ## since we pulled the entirety of fields from Rearrangements
########################################
seqdf = pd.DataFrame(all_seqid_resp)
keep = np.where(seqdf.apply(lambda col: any([c is not None and c != "" for c in col]), axis = 0)) 
seqdf_filt = seqdf.iloc[:,keep[0]]
seqdf_filt.to_csv("/home/ubuntu/projmnt/betts/coculture/hpap_irecptor_seqid_response.csv", sep = ",", index = False)

seqdf.columns[np.where([x in mixcr_airr.columns for x in seqdf.columns])].shape
mixcr_airr
########################################
# %% Convert to AnnData
########################################
import scirpy as ir
from scirpy.io._datastructures import AirrCell
from scirpy.io._convert_anndata import from_airr_cells

# %%
cells = []
for idx,row in enumerate(seqdf_filt.itertuples()):
    cell = AirrCell(cell_id = str(row.sequence_id))
    if row.locus == "TRB":
    # if "TRB" in row.v_gene:
        ac = AirrCell.empty_chain_dict()
        ac.update({
            "locus": row.locus,
            "junction": row.junction,
            "junction_aa": row.junction_aa,
            "v_call": row.v_call,
            "j_call": row.j_call,
            "consensus_count": 0,
            "productive": row.productive
            })
    elif row.locus == "TRA":
    # elif "TRA" in row.v_gene:
        raise ValueError(f"Found alpha chain at row {idx}")
    else:
        raise ValueError(f"Unknown chain at row {idx}")
    for field_idx in np.where([f not in ["locus", "junction_aa", "v_call", "j_call", "productive", "Index", "ir_substring"] for f in row._fields])[0]:
        cell[row._fields[field_idx]] = row.__getattribute__(row._fields[field_idx])
    cell.add_chain(ac)
    cells.append(cell)
adata = from_airr_cells(cells)

########################################
# %% Write out AnnData
########################################
# none2nan_obs = adata.obs.apply(lambda col: [np.nan if val is None else val for val in col])
adata.obs = adata.obs.drop("_33", axis = 1)
ir.pp.index_chains(adata)
adata.write_h5ad("/home/ubuntu/projmnt/betts/coculture/hpap_ireceptor.h5ad")

########################################
# %% Drop sequences that don't have any productive chains before running ir_dist
########################################
f = lambda x: ~ak.is_none(x["junction_aa"], axis=-1)
cif = lambda x: ~ak.is_none(x["VDJ"], axis=-1)
keep = np.where(cif(adata.obsm["chain_indices"])[:,0])[0].to_numpy()
drop = np.where(~cif(adata.obsm["chain_indices"])[:,0])[0].to_numpy()
adata_prod = adata[keep,:].copy()

########################################
# %% Run distance metrics
########################################
ir.pp.ir_dist(adata_prod, n_jobs = 14, sequence = "aa", key_added = "hpap_ireceptor_self_ident")
ir.pp.ir_dist(adata_prod, n_jobs = 14, sequence = "aa", key_added = "hpap_ireceptor_self_levenshtein", metric = "levenshtein")
adata_prod.write_h5ad("/home/ubuntu/projmnt/betts/coculture/hpap_irecpetor_filt_dist.h5ad")

# %%
adata_prod = anndata.read_h5ad("/home/ubuntu/projmnt/betts/coculture/hpap_ireceptor_filt_dist.h5ad")
adata_prod

vdjdb = ir.datasets.vdjdb()
vdjdb.obsm["chain_indices"][:,"VDJ"]

########################################
# %% Just making doubly sure that obs is in the same order as obsm
########################################
keep_seqid = adata.obs.iloc[keep,:].sequence_id
drop_seqid = adata.obs.iloc[drop,:].sequence_id
# %%
url = "https://hpap.ireceptor.org/airr/v1/rearrangement"
data = PostData(filters = [Filter(op = "in", content = {"field": "sequence_id", "value": keep_seqid.tolist()[0:500]})],
                first_result = 0, 
                size = None,
                fields = ["productive"]
                )
data.to_json()
response = request(url, data)
all([x["productive"] for x in response["Rearrangement"]])
# %%
url = "https://hpap.ireceptor.org/airr/v1/rearrangement"
data = PostData(filters = [Filter(op = "in", content = {"field": "sequence_id", "value": drop_seqid.tolist()[0:500]})],
                first_result = 0, 
                size = None,
                fields = ["productive"]
                )
data.to_json()
response = request(url, data)
all([not x["productive"] for x in response["Rearrangement"]])


########################################
# %% Pull Repertoire metadata based on the subjects involved
# ## Save to pkl becuase of the nested-ness (e.g., subject {subject.id, ...})
########################################
url = "https://hpap.ireceptor.org/airr/v1/repertoire"
data = PostData(filters = [Filter(op = "and", content = [
    Filter(op = "in", content = {"field": "subject.subject_id", "value": adata.obs.subject.unique().tolist()}),
    Filter(op = "=", content = {"field": "study.study_id", "value": "IR-Prak-000001"})]
                                  )],
                first_result = 0, 
                size = None
                )
data.to_json()
response = request(url, data)

subject_meta = pd.DataFrame(response["Repertoire"])
subject_meta.columns
subject_meta = subject_meta.iloc[:,np.where(~subject_meta.apply(is_null_column, axis = 0))[0]]
subject_meta.to_pickle(path = "/home/ubuntu/projmnt/betts/coculture/metadata/hpap_irecptor_repertoire_metadata.pkl")

