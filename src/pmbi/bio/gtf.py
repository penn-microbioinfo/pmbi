import re
import io

import copy
import pandas as pd


class FeatureFile(object):
    """
    A class to represent and manipulate feature data from a GTF file.

    Attributes:
    ----------
    _features : pd.DataFrame
        A DataFrame containing the features data.

    Methods:
    -------
    __setitem__(key, value)
        Sets the value of a feature in the DataFrame.
    
    __getitem__(key)
        Retrieves the value of a feature from the DataFrame.
    
    __repr__()
        Returns a string representation of the DataFrame.
    
    update(features: pd.DataFrame, inplace=False)
        Updates the features DataFrame with new data, optionally in place.
    
    attributes()
        Parses and returns the attributes GTF column as a DataFrame.
    
    get_attribute(key)
        Retrieves a specific attribute from the parsed attributes column.
    
    from_gtf(handle)
        Creates a FeatureFile object from a GTF file handle.
    """

    def __init__(self, features: pd.DataFrame):
        self._features = features

    def __setitem__(self, key, value):
        self._features[key] = value

    def __getitem__(self, key):
        return self._features[key]

    def __repr__(self):
        return self._features.__repr__()

    def update(self, features: pd.DataFrame, inplace=False):
        if inplace:
            self._features = features
        else:
            new = copy.deepcopy(self)
            new._features = features
            return new

    @staticmethod
    def _parse_attributes_string(attribute_str: str):
        attrs = {}
        for tag_value_pair in attribute_str.split(";"):
            tag_value_pair = tag_value_pair.strip()
            tag_value_spl = tag_value_pair.split(" ", maxsplit=1)
            if len(tag_value_spl) == 0:
                raise ValueError("Empty tag-value pair")
            else:
                if len(tag_value_spl) == 2:
                    tag,value = tuple(tag_value_spl)
                    attrs[tag] = re.sub("[\"']", "", value)

                # Add empty value in the event that there 
                # is only a key and str.strip gets rid of 
                # empty value
                elif len(tag_value_spl) == 1:
                    tag,value = (tag_value_spl[0], "")

                else:
                    # This should be impossible
                    raise ValueError("Length of tag-value split by '<space>' > 2")

        return attrs

    def attributes(self):
        return pd.DataFrame(self["attributes"].apply(self._parse_attributes_string).tolist())

    def get_attribute(self, key):
        attrs = self.attributes()
        return attrs.apply(lambda row: row[key] if key in row else pd.NA, axis=1)


    @staticmethod
    def from_gtf(handle):
        ldict = []
        for line in handle:
            if not line.startswith("#"):
                spl = [x.strip() for x in line.split("\t")]
                # if spl[2] == "gene":
                ann = [x.strip() for x in spl[8].split(";")]
                d = dict()
                d["seqname"] = spl[0]
                d["source"] = spl[1]
                d["feature"] = spl[2]
                d["start"] = int(spl[3])
                d["end"] = int(spl[4])
                d["score"] = spl[5]
                d["strand"] = spl[6]
                d["frame"] = spl[7]
                d["attributes"] = spl[8]
                ldict.append(d)
        return FeatureFile(features = pd.DataFrame(ldict))
        
def parse_gtf(handle: io.TextIOWrapper) -> pd.DataFrame:
    """
    Parses a GTF file from a given file handle and returns a DataFrame.

    Parameters
    ----------
    handle : io.TextIOWrapper
        A file handle to the GTF file to be parsed.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the parsed GTF data.
    """
    ldict = []
    for line in handle:
        if not line.startswith("#"):
            spl = [x.strip() for x in line.split("\t")]
            # if spl[2] == "gene":
            ann = [x.strip() for x in spl[8].split(";")]
            d = dict()
            d["seqname"] = spl[0]
            d["source"] = spl[1]
            d["feature"] = spl[2]
            d["start"] = int(spl[3])
            d["end"] = int(spl[4])
            d["score"] = spl[5]
            d["strand"] = spl[6]
            d["frame"] = spl[7]
            for ele in ann:
                ele_spl = ele.split(" ")
                if len(ele_spl) == 2:
                    k = ele_spl[0]
                    v = ele_spl[1]
                    d[k] = re.sub("[\"']", "", v)
            ldict.append(d)

    return pd.DataFrame(ldict)
# %%

