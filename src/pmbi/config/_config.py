from __future__ import annotations
from typing import Union, Iterable, Tuple
import copy

import toolz
from munch import Munch

from .item import *


class Config:
    def __init__(self, items: Union[list[Item], None] = None):

        self._index = {}
        if items is not None:
            for it in items:
                if it.toml_path in self._index:
                    raise ValueError(f"duplicate toml_path: {it.toml_path}")

                self._index[it.toml_path] = it

    def __setitem__(self, key, val):
        if key not in self._index:
            raise KeyError(f"path not found in Config: {key}")

        # Raises ValueError in case of mismatch
        self._index[key].set_value(val)

    def __str__(self):
        return "\n".join([f"{path}: {item.value}" for path,item in self.items()])

    def __getitem__(self, key):
        return self._index[key]

    def add_item(self, item: Item):
        if item.toml_path in self.as_dict():
            raise KeyError(f"Key alread exists in config: {item.toml_path}")
        self._index[item.toml_path] = item

    def items(self) -> Iterable[Tuple[str, Item]]:
        return self._index.items()

    def validate(self):
        for toml_path,item in self._index.items():
            if not item.is_set():
                if not item.is_optional():
                    raise ValueError(f"Required Item `{item.toml_path}` is None")

    def set_values_from_config(self, config: Config, inplace=True) -> Union[Config,None]:
        for toml_path,item in config.items():
            if toml_path in self._index:
                self[toml_path] = item.value
        if not inplace:
            return copy.deepcopy(self)

    def set_values_from_dict(self, d: dict, inplace=True) -> Union[Config,None]:
        for _p, item in self.items():
           item.set_value_from_dict(d) 

        if not inplace:
            return copy.deepcopy(self)

    def update_snakemake_config(self, smconfig) -> dict:
        for k,v in self.to_dict().items():
            if k in smconfig:
                raise ValueError(f"{k} already snakemake config dict with value: {smconfig[k]}")
            smconfig[k] = v

        return smconfig


    @staticmethod
    def from_items(items: list[Item]):
        return Config(items)

    def to_dict(self):
        def _nodewise(*toml_paths, config, preceding=None):
            parts = []
            for s in toml_paths:
                spl = s.split(".", maxsplit=1)
                parts.append({"first": spl[0], "rest": spl[1]})
            d = {}
            for key, grp in toolz.groupby(lambda parts: parts["first"], parts).items():
                d[key] = {}
                if preceding is None:
                    prec_key_repr = key
                    next_prec = [key]
                else: 
                    prec_key_repr = ".".join(preceding)
                    next_prec.append[key]
                next_strings = [parts["rest"] for parts in grp]
                more_levels = ["." in ns for ns in next_strings]
                if not any(more_levels):
                    for ns in next_strings:
                        k = f"{prec_key_repr}.{ns}"
                        d[key][ns] = config[k].value 
                elif all(more_levels):
                    d[key][ns] = _nodewise(*next_strings, config=config, preceding=next_prec)
                else:
                    raise ValueError("Invalid config spec with uneven levels")
            return d

        config_dict = _nodewise(*self._index, config=self)

        return config_dict

    def to_munch(self):
        return Munch.fromDict(self.to_dict())

    def from_strings(*strings: str, split_char="."):
        parts = []
        for s in strings:
            spl = s.split(split_char, maxsplit=1)
            parts.append({"first": spl[0], "rest": spl[1]})
        d = {}
        for key, grp in toolz.groupby(lambda parts: parts["first"], parts).items():
            next_strings = [parts["rest"] for parts in grp]
            more_levels = [split_char in ns for ns in next_strings]
            if not any(more_levels):
                d[key] = next_strings
            elif all(more_levels):
                d[key] = MultiNode.from_strings(*next_strings)
            else:
                raise ValueError("Invalid config spec with uneven levels")
        return d

    def into_toml(self):
        raise NotImplementedError
# %%

