import argparse
import pmbi.wrappers.scvi as pmbscvi

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--adata", help = "Name of input h5ad file.")
parser.add_argument("-o", "--output", help = "Name of output h5ad file.")
parser.add_argument("-m", "--model", help = "Path to trained, saved model dir.")
parser.add_argument("-n", "--n_latent", type = int, help = "Number of latents.")
args = parser.parse_args()

m = pmbscvi.ScviModeler()
m.read(args.adata)
m.load_model(args.model)
m.get_normalized_expr(n_latent = args.n_latent, n_samples = 5)
m.write(m.fp.name.replace(m.ext, f"normExpr.{m.ext}"))
