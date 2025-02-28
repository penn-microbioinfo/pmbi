import palettable
from matplotlib.colors import rgb2hex
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--palette", action = "store", required=True, help = "Name of palette.")
parser.add_argument("-s", "--s", default = 6, action = "store", help = "Number of colors from palette - check palettable docs for maximums.")
parser.add_argument("-r", "--reverse", action = "store_true", help = "Reverse the palette.")
parser.add_argument("--outfmt", action = "store", default = "r", help = "Output format for color palette.")
args = parser.parse_args()

pal_parts = args.palette.split('.') 
obj = palettable
for part in pal_parts:
    obj = (getattr(obj, part))

hex = [rgb2hex(c) for c in obj.mpl_colors]
quoted = [f"'{c}'" for c in hex]
if args.outfmt == 'r':
    print(f"c({', '.join(quoted)})")
