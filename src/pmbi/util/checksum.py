"""
checksum.py

This module provides utilities for computing checksums of files using the MD5 hashing algorithm. It supports computing checksums for individual files and multiple files in parallel. The script can be executed directly to compute checksums for files in a specified directory.

Functions:
- md5sum(stream): Computes the MD5 checksum of a given stream.
- md5sum_file(path): Computes the MD5 checksum for a file at a specified path.
- checksum_multi(l, f, n_jobs): Computes checksums for a list of files in parallel.

Execution:
When executed as a script, it computes the MD5 checksums for files in a given directory. The user can specify the checksum type, directory path, recursion option, and number of parallel jobs through command-line arguments.

Dependencies:
- hashlib: For computing MD5 checksums.
- pathlib: For handling file paths.
- typing: For type annotations.
- joblib: For parallel processing.

Usage:
Run the script from the command line with the required checksum type and path. Use the optional arguments for recursion and parallel processing.

Example:
python checksum.py md5 /path/to/directory --recursive --n_jobs 4
"""

import hashlib
import typing
from pathlib import Path

from joblib import Parallel, delayed


def md5sum(stream):
    """
    Computes the MD5 checksum of a given stream.

    Args:
        stream: A file-like object that supports reading in binary mode.

    Returns:
        str: The MD5 checksum as a hexadecimal string.
    """
    m = hashlib.new("md5")
    m.update(stream.read())
    return m.hexdigest()


def md5sum_file(path: Path) -> tuple[str, str]:
    """
    Computes the MD5 checksum for a file at the specified path.

    Args:
        path (Path): The path to the file for which the MD5 checksum is computed.

    Returns:
        tuple[str, str]: A tuple containing the MD5 checksum as a hexadecimal string and the string representation of the file path.
    """
    with open(path, "rb") as stream:
        hexdig = md5sum(stream)
    return (hexdig, str(path))


def checksum_multi(
    l: typing.Iterable[Path], f: typing.Callable, n_jobs=1
) -> typing.Any:
    """
    Computes checksums for a list of files in parallel.

    Args:
        l (Iterable[Path]): An iterable of file paths for which checksums are computed.
        f (Callable): A function that computes the checksum for a single file.
        n_jobs (int, optional): The number of parallel jobs to run. Defaults to 1.

    Returns:
        Any: A list of checksum results for the provided file paths, as computed by the function `f`.
    """
    hexdigs = Parallel(n_jobs=n_jobs)(delayed(f)(p) for p in l)
    return hexdigs


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("checksum_type", help="Type of checksum.")
    parser.add_argument("path", help="Path to directory containing files to checksum")
    parser.add_argument(
        "-r", "--recursive", action="store_true", help="Generate checksums recursively"
    )
    parser.add_argument(
        "--n_jobs", type=int, help="Number of checksumming jobs to spawn"
    )
    args = parser.parse_args()

    match args.checksum_type:
        case "md5":
            if args.recursive:
                l = [p for p in Path(args.path).rglob("*") if not p.is_dir()]
                out = checksum_multi(l, md5sum_file, n_jobs=args.n_jobs)
                print("\n".join([f"{hexdig}  {fpath}" for hexdig, fpath in out]))
            else:
                raise NotImplementedError
        case _:
            raise NotImplementedError
