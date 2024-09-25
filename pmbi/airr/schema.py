import logging

import numpy as np
import numpy.typing
import pandas as pd
import requests
import yaml


def get_airr_schema(schema: str) -> dict:
    r = requests.get(
        "https://raw.githubusercontent.com/airr-community/airr-standards/refs/heads/master/specs/airr-schema.yaml"
    )
    if r.ok:
        airr_schema = yaml.load(r.text, Loader=yaml.Loader)
        if schema in r.text:
            return airr_schema[schema]
        else:
            raise ValueError(f"Schema not present in response text: {schema}")
    else:
        raise requests.exceptions.HTTPError(
            "AIRR schema request returned NOT OK. Check the link."
        )


def get_imgt_locus_names() -> numpy.typing.ArrayLike:
    return np.array(
        [
            "TRA",
            "TRG",
            "IGL",
            "IGK",
            "TRB",
            "TRD",
            "IGH",
        ]
    )


def calls_to_imgt_locus_names(calls: list[str]) -> str:
    locus_names = get_imgt_locus_names()
    matches_all_calls = []
    calls = [c for c in calls if not pd.isnull(c)]
    for c in calls:
        matches = locus_names[[l in c for l in locus_names]]
        if len(matches) == 0:
            logging.warning(f"Call does not match any known IMGT locus names: {c}")
        elif len(matches) > 1:
            raise ValueError(f"Call matches multiple IMGT locus names: {c}")
        else:
            matches_all_calls.append(matches[0])

    if len(set(matches_all_calls)) > 1:
        raise ValueError(f"Calls do not match the same locus: {calls}")
    else:
        return matches_all_calls[0]

