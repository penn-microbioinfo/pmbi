import argparse

parser = argparse.ArgumentParser(description="Compare two lists of checksums (in files). Requires that the checksums are in newline-separated lists where each line is a space-separated string with the checksum in the first column and the filename in the second column.")
parser.add_argument("first", action = "store", help = "First list of checksums.")
parser.add_argument("second", action = "store", help = "Second list of checksums.")
args = parser.parse_args()

first = {}
with open(args.first, 'r') as fs:
    for line in fs:
        spl = line.split()
        first[spl[1]] = spl[0]

second = {}
with open(args.second, 'r') as fs:
    for line in fs:
        spl = line.split()
        second[spl[1]] = spl[0]

for fk,fv in first.items():
    if fk not in second:
        print(f"Missing filename: {fk}")
    else:
        if fv != second[fk]:
            print(f"Mismatched md5 checksum for file {fk}: {fv} -- {second[fk]}")

