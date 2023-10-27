import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn.objects as so
import anndata
import mudata
import pandas as pd
import argparse
import torch
import copy

class Modeler(object):
    def __init__(self, batch_key = "orig_ident"):
        self.batch_key = batch_key
        self.data = None
        self._reader = None
        self._writer = None

    def _io_funs(self, fp):
        if self.ext == "h5ad":
            return (anndata.read_h5ad, anndata.write_h5ad)
        elif self.ext == "h5mu":
            return (anndata.read_h5mu, anndata.write_h5mu)
        else:
            raise ValueError(f"Unexpected file extension: {fp.suffix}")

    def read(self, fp):
        self.fp =pathlib.Path(fp)
        self.ext = fp.replace('.', '')
        self._reader, self._writer = self._io_funs(self.fp)
        self.data = self._reader(self.fp)

    def write(self):
        self._writer(self.fp.name.replace(self.ext, f"latentReprs.{self.ext}"))

class ScviModeler(Modeler):
    def __init__(self, batch_key = "orig_ident"):
        super().__init__(batch_key)
        self._type_ = scvi.model.SCVI

    def setup_data(self):
        if isinstance(self.data, anndata._core.anndata.AnnData):
            self._type_.setup_anndata(self.data, batch_key = self.batch_key)

    def train(self, n_latent_values: list):
        for nlv in n_latent_values:
            model = self._type_(self.data, n_latent = nlv)
            model.train()
            obsm_key = f"X_totalvi_n_latent_{nlv}"
            self.data.obsm[obsm_key] = model.get_latent_representation()
            model.save(f"model_{obsm_key}", overwrite=True)

def de_volcano_plot(csv_path, fdr_target = 0.05, title=""):
    sample_name = os.path.basename(csv_path).split('_')[0]
    de = pd.read_csv(csv_path)
    de = de.set_index(de.columns[0])
    fdr_key =f"is_de_fdr_{fdr_target}"

    sp = so.Plot(data=de, x="lfc_mean", y="proba_de", color=fdr_key)\
            .add(so.Dots(pointsize=1), legend=False)\
            .scale(color=["red", "black"])\
            .label(x = "Log2FoldChange", y = f"p(DE), FDR = {fdr_target}", title = title)\
            .theme({"axes.facecolor": "w", "axes.edgecolor": "#000000"})
    print(plt.legend())
    '''
    return sp
    scat = sns.scatterplot(data = de, x = "lfc_mean", y = "proba_de", hue = fdr_key, palette = ["black", "red"], s=2, linewidth=0)
    scat.text(2, 0.50, f"DE Genes: {de.loc[de[fdr_key],].shape[0]}",)
    scat.set(xlabel = "Log2FoldChange", ylabel = f"p(DE), FDR = {fdr_target}", title = f"{sample_name} - {de['comparison'][0]}")
    plt.legend([],[], frameon=False)
    '''
    return copy.deepcopy(sp)
    #return copy.deepcopy(scat.get_figure())

