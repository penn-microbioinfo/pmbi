r"""
zenodo.py

This module provides utilities for preparing Zenodo deposition metadata.
It currently supports converting a plain-text list of authors into the
Zenodo `creators` JSON record format, using a user-supplied regular
expression to split each line into given/family name.

Functions:
- authors_text_to_json(text, pattern): Converts a newline-delimited author
  list into a list of Zenodo-style `person_or_org` JSON records, with all
  fields other than the parsed name left blank.

Dependencies:
- re: For matching author lines against a user-supplied pattern.

Usage:
    import re
    from pmbi.data_submission.zenodo import authors_text_to_json

    pattern = re.compile(r"^(?P<family_name>[^,]+),\s*(?P<given_name>.+)$")
    records = authors_text_to_json(author_list_text, pattern)
"""

import re
from typing import Any

from pmbi.logging import streamLogger

logger = streamLogger(__name__)


def authors_text_to_json(text: str, pattern: re.Pattern) -> list[dict[str, Any]]:
    """
    Converts a plain-text list of authors (one per line) into Zenodo-style
    JSON records, with every field other than the parsed name left blank.

    Each output record has the form:

        {
          "person_or_org": {
            "type": "personal",
            "family_name": "",
            "given_name": "",
            "identifiers": []
          },
          "affiliations": []
        }

    Args:
        text (str): The plain-text author list, one author per line.
        pattern (re.Pattern): Compiled regular expression used to parse each
            line. It should define the named capture groups `given_name`
            and/or `family_name`; either may be omitted if a line only
            encodes a single name component.

    Returns:
        list[dict[str, Any]]: One JSON-serializable record per author line
        that matched `pattern`. Lines that fail to match are skipped and
        logged as warnings.
    """
    records: list[dict[str, Any]] = []

    for line_num, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue

        match = pattern.match(line)
        if match is None:
            logger.warning(f"Line {line_num} did not match pattern: {line!r}")
            continue

        groups = match.groupdict()
        family_name = (groups.get("family_name") or "").strip()
        given_name = (groups.get("given_name") or "").strip()

        record: dict[str, Any] = {
            "person_or_org": {
                "type": "personal",
                "family_name": family_name,
                "given_name": given_name,
                "identifiers": [],
            },
            "affiliations": [],
        }
        records.append(record)

    return records
