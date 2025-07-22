import argparse
from pathlib import Path
from pmbi.util.misc import symlink_is_valid

parser = argparse.ArgumentParser()
parser.add_argument("directory", nargs = 1, help="Directory within which to check for broken symlinks.")
args = parser.parse_args()

d = Path(args.directory[0])
if not d.exists():
    raise ValueError(f"Path does not exist: {d}")

for p in d.iterdir():
    if not symlink_is_valid(p):
        print(f"Broken symlink: {p}")


