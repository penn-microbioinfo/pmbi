import sys
import anndata
import mudata
import pandas as pd
import argparse
import torch

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
