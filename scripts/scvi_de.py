import pmbi.wrappers.scvi as pmbint
import argparse
import importlib
import os
importlib.reload(pmbint)

def parse_comparisons(string):
    comparisons = []
    for pair in string.split(';'):
        comparisons.append(tuple(pair.split(',')))
    return comparisons

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--data", help = "Name of input h5ad/h5mu file.")
parser.add_argument("-m", "--model", help = "Path to model directory.")
parser.add_argument("-g", "--groupby", help = "Obs column to group by.")
parser.add_argument("-f", "--factorial", action = "store_true", help = "Pass flag to run deg analysis for factorial combinations of groupby categories.")
parser.add_argument("-c", "--comparisons",  help = "Pass flag to run deg analysis for factorial combinations of groupby categories.")
parser.add_argument("--fdr_target",  default=0.05, type = float, help = "FDR cutoff for differential expression.")
args = parser.parse_args()

sm = pmbint.ScviModeler.load(args.data, args.model)
sm.differential_expression(groupby=args.groupby, factorial = args.factorial, fdr_target=args.fdr_target, comparisons = parse_comparisons(args.comparisons))
