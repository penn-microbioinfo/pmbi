import anndata
import logging
import pmbi.wrappers.scanpy as scp
import typing
import itertools
import scipy
import numpy as np
import mudata
import argparse
import torch
import scvi
import pathlib
import pandas as pd
import os
import re

class Modeler(object):
    def __init__(self, batch_key = "orig_ident"):
        self.batch_key = batch_key
        self.sample_name = None
        self.data = None
        self._reader = None
        self._writer = None
        self.models = {}

    def _io_funs(self, fp):
        if self.ext == "h5ad":
            return (anndata.read_h5ad, anndata._io.h5ad.write_h5ad)
        elif self.ext == "h5mu":
            return (mudata.read_h5mu, mudata._core.mudata.MuData.write_h5mu)
        else:
            raise ValueError(f"Unexpected file extension: {self.fp.suffix}")

    def read(self, fp):
        self.fp =pathlib.Path(fp)
        self.ext = self.fp.suffix.replace('.', '')
        self.sample_name = self.fp.name.replace(self.ext, '')
        self._reader, self._writer = self._io_funs(self.fp)
        self.data = self._reader(self.fp)

        if not scp.obs_names_unique(self.data):
            logging.warning("Makeing obs_names unique and saving.")
            self.data.obs_names_make_unique()
            self.data.write(self.fp)

    def write(self, outfp):
        if self._writer is not None:
            self._writer(outfp, self.data)
        else:
            raise IOError(f"No _writer method available for {self._type_}")

class ScviModeler(Modeler):
    def __init__(self, batch_key = "orig_ident"):
        super().__init__(batch_key)
        self._type_ = scvi.model.SCVI

    def setup_data(self, layer = None):
        if isinstance(self.data, anndata._core.anndata.AnnData):
            self._type_.setup_anndata(self.data, batch_key = self.batch_key, layer = layer)

    def train(self, n_latent: int):
        model = self._type_(self.data, n_latent = n_latent)
        model.train()
        self.models[n_latent] = model

    def get_latent_repr(self, n_latent: int, as_obsm=True):
        obsm_key = f"X_scvi_n_latent_{n_latent}"
        if as_obsm:
            self.data.obsm[obsm_key] = self.models[n_latent].get_latent_representation()
        else:
            return self.models[n_latent].get_latent_representation()

    def get_normalized_expr(self, n_latent: int, as_layer=True, n_samples = 1): 
        if as_layer:
            self.data.layers[f"normExp_scvi_n_latent_{n_latent}"] = self.models[n_latent].get_normalized_expression(n_samples = n_samples)
        else:
            return self.models[n_latent].get_normalized_expression()

    def differential_expression(self, groupby: str, factorial: bool = False, comparisons: typing.Optional[list[tuple[str]]] = None, fdr_target = 0.05) -> dict:
        n_latent = list(self.models.keys())[0]
        if factorial and comparisons is not None:
            raise ValueError("factorial and comparisons cannot be specified at the same time")
        if factorial:
            for c in itertools.combinations(self.data.obs[groupby].unique(), r=2):
                de = self.models[n_latent].differential_expression(adata=self.data, groupby=groupby, group1=c[0], group2=c[1], fdr_target = fdr_target)
                de.to_csv(f"{self.sample_name}_scvi_deg_{groupby}_{c[0]}-vs-{c[1]}_nlatent{n_latent}.csv")
        else:
            if comparisons is not None:
                for c in comparisons:
                    de = self.models[n_latent].differential_expression(adata=self.data, groupby=groupby, group1=c[0], group2=c[1], fdr_target = fdr_target)
                    de.to_csv(f"{self.sample_name}_scvi_deg_{groupby}_{c[0]}-vs-{c[1]}_nlatent{n_latent}.csv")
            else:
                de = self.models[n_latent].differential_expression(adata = self.data, groupby = groupby, fdr_target = fdr_target)
                de.to_csv(f"scvi_deg_{groupby}_all-vs-rest.csv")

        return de

    def load_model(self, model):
        n_latent = int(re.search("[_]n[_]latent[_]([0-9]+)[_]", os.path.basename(model)).group(1))
        self.models[n_latent] = self._type_.load(model, adata=self.data)
        print(self.models)

    @staticmethod
    def load(fp, model):
        sm = ScviModeler()
        sm.read(fp)
        n_latent = int(re.search("[_]n[_]latent[_]([0-9]+)[_]", os.path.basename(model)).group(1))
        sm.models[n_latent] = sm._type_.load(model, adata=sm.data)
        return sm

    def save_model(self, n_latent):
        outname=f"{self.sample_name}_n_latent_{n_latent}_model"
        self.models[n_latent].save(outname, overwrite=True)
'''
    def get_normalized_counts(self, model, n_chunks=1):
        if not isinstance(model, scvi.model.SCVI):
            model = scvi.model.SCVI.load(model, adata = self.data)
        chunks = np.array_split(np.arange(0, len(self.data.obs_names)), n_chunks)
        gex_normexp = scipy.sparse.csc_array((1,len(self.data.var_names)), dtype=np.float32)
        for chunk in chunks:
            normexp_this_chunk = model.get_normalized_expression(adata = self.data, library_size="latent", indices=chunk)
            gex_normexp = scipy.sparse.vstack( (gex_normexp, normexp_this_chunk.astype(pd.SparseDtype("float32",0)).sparse.to_coo().tocsc()) )
        return gex_normexp[1:,:]
'''
