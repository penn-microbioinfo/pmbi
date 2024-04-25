import argparse
import numpy as np
import pmbi.wrappers.scvi as pmbscvi
import scvi
import anndata

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-m", "--model", help = "Path to trained, saved model dir.")
parser.add_argument("-c", "--ncells", type = int, help = "Number of cells to use for log likelihood estimation.")
args = parser.parse_args()


adata = anndata.read_h5ad(args.adata)
trainer = scvi.model.SCVI.load(args.model, adata = adata) 
cellidx = np.random.choice(np.arange(0,adata.shape[0]), size = args.ncells, replace = False)
margll = trainer.get_marginal_ll(indices = cellidx)
print(args.model, margll)
