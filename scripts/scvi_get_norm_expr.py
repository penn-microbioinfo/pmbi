import os
import argparse
import pmbi.wrappers.scvi as pmbscvi
import scipy.sparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-m", "--model", help = "Path to trained, saved model directory.")
parser.add_argument("-c", "--chunksize", default = None, help = "Number of cells to compute expression values for at a time. Default: all the cells")
parser.add_argument("-n", "--n_samples", default = 10,  help = "Number of posterior samples upon which to base predictions. Default: 10")
parser.add_argument("-l", "--library_size", default = "latent", help = "Library size to scale expression to. Default: 'latent'")
parser.add_argument("-o", "--output_path", required = True, help = "Output path to write the normalized count matrix.")
args = parser.parse_args()

output_path = args.output_path or f"./{os.path.basename(args.model)}_normcounts.npz"
modeler = pmbscvi.ScviModeler(args.data)
modeler.load_model(args.model)
normexpr = modeler.get_normalized_expression(
        chunksize = args.chunksize,
        n_samples = args.n_samples,
        library_size = args.library_size
        )
scipy.sparse.save_npz(file=output_path, matrix=normexpr)
