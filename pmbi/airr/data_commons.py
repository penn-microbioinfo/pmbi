# %%
import json

from requests_html import HTMLSession


# %%
class Filter(object):
    def __init__(self, op: str, content: list|dict):
        self.op = op
        self.content = content
        self._check_types()
        if isinstance(self.content, list):
            self.content = [f.to_dict() for f in self.content]
        else:
            self.content = self.content
    def _check_types(self):
        if isinstance(self.content, list):
            assert all([isinstance(f, Filter) for f in self.content]), "Content lists for Filter must be made up of other Filter instances"
        elif isinstance(self.content, dict):
            pass
        else:
            raise ValueError(f"Invalid content type: {type(self.content)}")
    def to_dict(self) -> dict:
            return {
                    "op": self.op,
                    "content": self.content
                    }

# %%
class PostData(object):
    def __init__(
        self, filters: list[Filter] = [], first_result: int = 0, size: int|None = None, fields: list[str] = [], format: str = "json"
    ):
        self.filters = filters
        self.first_result = first_result
        self.size = size
        self.fields = fields
        self.format = format
    def to_dict(self) -> dict:
        if len(self.filters) == 0:
            filters_repr = None
        elif len(self.filters) > 1:
            filters_repr = [f.to_dict() for f in self.filters]
        else:
            filters_repr = self.filters[0].to_dict()
        return {
            "filters": filters_repr,
            "from": self.first_result,
            "fields": self.fields,
            "size": self.size,
            "format": self.format
        }
# %%
def request(url: str, data: PostData):
    session = HTMLSession()
    response = session.post(url, json=data.to_dict(), headers = {"content-type": "applications/json"})
    session.close()
    return json.loads(response.text)

# %%
data = PostData(
    filters=[
        Filter(op = "and",
               content = [
                   Filter(op = "=", content = {"field": "sample.tissue.label", "value": "spleen"}),
                   Filter(op = "=", content = {"field": "study.study_id", "value": "IR-Prak-000001"}),
                   ])

    ],
    first_result=0,
    size=None,
    fields=["repertoire_id", "sample.sample_id"],
)
url = "https://hpap.ireceptor.org/airr/v1/repertoire"
response = request(url, data)
repids = [rec["repertoire_id"] for rec in response["Repertoire"]]

# %%
url = "https://hpap.ireceptor.org/airr/v1/rearrangement"
data = PostData(filters = [Filter(op = "in", content = {"field": "repertoire_id", "value": repids[1:25:1]})], 
                first_result = 0, 
                size = None,
                fields = ["sequence_id", "productive", "v_call", "j_call", "clone_id", "sequence", "v_cigar", "j_cigar", "sequence_alignment", "v_alignment_start", "v_alignment_end","j_alignment_start", "j_alignment_end"])
data.to_dict()
response = request(url, data)
response["Rearrangement"][1]
sa = response["Rearrangement"][15358]["sequence_alignment"]
s = response["Rearrangement"][15358]["sequence"]
[r["sequence_alignment"][::-1][1:50] for r in response["Rearrangement"][1:20]]
sa
s[140:(140+151)]
len(s[(140+151):])

# %%
url = "https://hpap.ireceptor.org/airr/v1/clone"
data = PostData(filters = [],
        first_result = 0, 
                size = 25,
                fields = ["clone_id"])
data.to_dict()
session = HTMLSession()
response = session.post(url, json={"from": 0, "size": 25, "fields": ["clone_id"]}, headers = {"content-type": "applications/json"})
session.close()
response.text
response = request(url, {"from": 0, "size": 25, "fields": ["clone_id"]}
response["Clone"]
# %%

