import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("json", action = "store", help = "Path to the json file to print pretty-like.")
args = parser.parse_args()

parsed = json.load(open(args.json, "r"))
print(json.dumps(parsed, indent = 2))



