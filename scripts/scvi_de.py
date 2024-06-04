import os
import argparse
import pmbi.wrappers.scvi as pmbscvi
import pathlib
import pdb
import scipy.sparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-m", "--model", help = "Path to trained, saved model directory.")
parser.add_argument("-g", "--groupby", help = "Variable in AnnData.obs to group cells by.")
parser.add_argument("-f", "--factorial", action = "store_true", help = "Pass flag to test factorial combinations of the unique values of groupby variable.")
parser.add_argument("-o", "--output_path", required = True, help = "Output path to write the table of differentially expressed genes.")
parser.add_argument("--fdr_target", type = float, required = False, default = 0.05, help = "FDR target to pass to SCVI.differential_expression.")
args = parser.parse_args()

output_path = args.output_path or f"./{os.path.basename(args.model)}_deg.csv"
modeler = pmbscvi.ScviModeler(pathlib.Path(args.adata))
modeler.load_model(pathlib.Path(args.model))
deg = modeler.differential_expression(
        groupby = args.groupby,
        factorial = args.factorial,
        comparisons = None,
        fdr_target = args.fdr_target,
        )
deg = deg.reset_index(names="gene")
deg.to_csv(output_path, sep = ",", index = False)
