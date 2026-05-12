import re
from pathlib import Path
from typing import Any, Union, _UnionGenericAlias, get_args

from munch import Munch


class Item:
    def __init__(self, toml_path: str, type_=str, optional: bool = False, default=None):
        self.toml_path = toml_path
        self._split_char = "."
        self.toml_path_spl = toml_path.split(self._split_char)
        self.type_ = type_
        self.optional = optional
        self.default = default

        self.value = default

    def check_value(self, val):
        if isinstance(self.type_, _UnionGenericAlias):
            types_ = get_args(self.type_)
            type_check_pass = False
            for t in types_:
                if isinstance(val, t):
                    type_check_pass = True
            if not type_check_pass:
                raise TypeError(
                    f"val type mismatch; expected one of {types_}, got {type(val)}"
                )
        else:
            if not isinstance(val, self.type_):
                raise TypeError(
                    f"val type mismatch; expected {self.type_}, got {type(val)}"
                )

    def is_set(self):
        if self.value is None:
            return False
        else:
            return True

    def is_optional(self):
        return self.optional

    def on_set_value(self, val):
        return val

    def set_value(self, val):
        self.check_value(val)
        self.value = self.on_set_value(val)

    @property
    def name(self):
        return self.toml_path_spl[::-1][0]

    def set_value_from_dict(self, d: dict) -> None:
        pos = d.copy()
        for key in self.traverse_path():
            if key in pos:
                pos = pos[key]
            else:
                return None

        self.set_value(pos)
        
    # def from_toml(self, config: dict[str, Any]):
    #     try:
    #         v = config
    #         for i in self.toml_path_spl:
    #             v = v[i]
    #
    #         if not isinstance(v, self.type_):
    #             raise TypeError(
    #                 f"Unexpected type for config option at `{self.toml_path}`: {type(v)}"
    #             )
    #
    #         return v
    #
    #     except KeyError:
    #         raise KeyError(f"Missing expected config value at: {self.toml_path}")

    def as_dict(self):
        key, val = tuple(self.toml_path.split(self._split_char, maxsplit=1))
        if self._split_char in val:
            return {key: Item(val).as_dict()}
        else:
            return {key: val}

    def traverse_path(self):
        for s in self.toml_path_spl:
            yield s

    def update(self, dict_to_update):
        spl = self.toml_path.split(self._split_char)
        first = spl[0]
        rest = spl[1:]
        if len(rest) > 1:
            dict_to_update[first] = {rest[0]: None}
            sub_i = Item(".".join(rest))
            sub_i.update(dict_to_update[first])
        else:
            if first in dict_to_update:
                dict_to_update[first].append(rest[0])
            else:
                dict_to_update[first] = [rest[0]]


class RegexItem(Item):
    def __init__(
        self,
        toml_path: str,
        type_=str,
        optional=False,
        default=None,
        capture_groups: int = 0,
    ):
        super().__init__(toml_path, type_, optional, default)
        self.capture_groups = capture_groups

    def check_value(self, val):
        super().check_value(val)
        val_c = re.compile(val)
        if val_c.groups != self.capture_groups:
            raise ValueError(
                f"Expected {self.capture_groups} capture groups, got {val_c.groups}, for config regex pattern `{self.toml_path}`"
            )

    def on_set_value(self, val):
        val = super().on_set_value(val)
        return val

    def set_value(self, val):
        super().set_value(val)

    def from_toml(self, config: dict[str, Any]):
        v = super().from_toml(config)

        return v


class PathItem(Item):
    def __init__(self, toml_path: str, optional=False, default=None, must_exist=False):
        self.type_ = Union[Path, str]
        super().__init__(toml_path, self.type_, optional, default)
        self.must_exist = must_exist

    def check_value(self, val):
        super().check_value(val)
        if self.must_exist:
            if not Path(val).exists():
                raise OSError(f"Path value for PathItem does not exist: {val}")

    def on_set_value(self, val):
        val = super().on_set_value(val)
        return Path(val)

    def set_value(self, val):
        super().set_value(val)

    def from_toml(self, config: dict[str, Any]):
        v = Path(super().from_toml(config))
        return v


class ChoiceItem(Item):
    def __init__(
        self,
        toml_path: str,
        choices: list[Any],
        type_=str,
        optional=False,
        default=None,
    ):
        super().__init__(toml_path, type_, optional, default)
        self.choices = choices

    def check_value(self, val):
        super().check_value(val)
        if val not in self.choices:
            raise ValueError(f"Value for ChoiceItem not in choices: found {val}, expected one of {self.choices}")

