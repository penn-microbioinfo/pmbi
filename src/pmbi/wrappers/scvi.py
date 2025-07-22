import argparse
import copy
import itertools
import logging
import os
import re
import typing
from functools import wraps
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import mudata
import numpy as np
import pandas as pd
import scipy
import scvi
import seaborn as sns
import seaborn.objects as so
import torch
from pmbi.anndata.util import obs_names_unique

import pmbi.wrappers.scanpy as scp


def _is_model_loaded(inner):
    @wraps(inner)
    def wrapper(self, *args, **kwargs):
        if self.model is None:
            raise ValueError(
                "Model has not been loaded yet, so cannot run method that requires a model."
            )
        else:
            return inner(self, *args, **kwargs)

    return wrapper


# %%
class Modeler(object):
    def __init__(self, datafile: Path, batch_key="orig_ident"):
        self.batch_key = batch_key
        self.datafile = datafile
        self.type_ = None
        self.ext = self.datafile.suffix.replace(".", "")
        self.sample_name = self.datafile.name.replace(self.ext, "")
        self._reader, self._writer = self._io_funs()
        self.model = None
        self.model_n_latent = None

        self.data = self._reader(self.datafile)
        if not obs_names_unique(self.data):
            raise ValueError("obs_names must be unique before adding to Modeler")

    def _io_funs(self):
        if self.ext == "h5ad":
            return (anndata.read_h5ad, anndata._io.h5ad.write_h5ad)
        elif self.ext == "h5mu": return (mudata.read_h5mu, mudata._core.mudata.MuData.write_h5mu)
        else:
            raise ValueError(f"Unexpected file extension: {self.datafile.suffix}")

    def train(self, n_latent: int):
        self.model_n_latent = n_latent
        self.model = self.type_(self.data, n_latent=n_latent)
        self.model.train()

    def write(self, path=None):
        if self.model is None:
            raise ValueError(
                "Cannot write a model that does not exist. Train a model first"
            )

        path = path or f"./{self.sample_name}_n_latent_{self.model.n_latent}_model"

        self.model.save(outname, overwrite=True)

    def load_model(self, model):
        self.model = self.type_.load(model, adata=self.data)

    def write_data(self, path):
        self._writer(path, self.data)


def dataframe_into_csc(df: pd.DataFrame) -> scipy.sparse.csc_array:
    """
    Convert a pandas DataFrame into a compressed sparse column matrix.

    Parameters:
        df (pd.DataFrame): The pandas DataFrame to convert into a CSC matrix.

    Returns:
        scipy.sparse.csc_matrix: The resulting CSC matrix.

    Example:
        # Assuming 'dataframe' is a pandas DataFrame
        csc_matrix = dataframe_into_csc(dataframe)
    """
    return df.astype(pd.SparseDtype("float32", 0)).sparse.to_coo().tocsc()


class ScviModeler(Modeler):
    def __init__(self, datafile, batch_key="orig_ident"):
        super().__init__(datafile, batch_key)
        self.type_ = scvi.model.SCVI

    def setup_data(self, layer=None):
        if isinstance(self.data, anndata.AnnData):
            self.type_.setup_anndata(self.data, batch_key=self.batch_key, layer=layer)

    def get_latent_repr(self, n_latent: int, as_obsm=True):
        obsm_key = f"X_scvi_n_latent_{n_latent}"
        if as_obsm:
            self.data.obsm[obsm_key] = self.model.get_latent_representation()
        else:
            return self.model.get_latent_representation()

    @_is_model_loaded
    def get_normalized_expression(
        self,
        chunksize: int | None = None,
        n_samples: int = 10,
        library_size: str | int = "latent",
        **kwargs,
    ):
        """
        Get the normalized gene expression data using the specified SCVI model.

        Parameters:
            self (pmbi.wrappers.scvi.ScviModeler)
            model (scvi.model.SCVI): The trained SCVI model used for normalization. If not provided, will default to loaded model, if it exists.
            chunksize (Optional[int]): The number of cells to process simultaneously. If None, process all cells at once. Defaults to None.
            n_samples (int): The number of samples to draw from the posterior predictive distribution. Defaults to 10.
            library_size (Union[str, int]): Method for library size normalization. Defaults to "latent".
                                             If "latent", uses the inferred library size from the latent space.
                                             If int, uses the provided value as the library size.
            **kwargs: Additional keyword arguments to pass to the SCVI model's `get_normalized_expression` method.

        Returns:
            scipy.sparse.csc_matrix: The normalized gene expression data in a sparse matrix format.

        Raises:
            ValueError:
                If model is None, loaded model is not trained, or  no model is currently loaded (ie self.model is also None)
                If the provided model is not an instance of `scvi.model.SCVI`.

        Example:
            # Assuming 'model' is a trained SCVI model
            norm_expression = obj.get_normalized_expression(model, chunksize=100, n_samples=5)
        """

        if not isinstance(model, scvi.model.SCVI):
            raise ValueError

        if not model.is_trained:
            raise ValueError

        cell_indices = np.arange(0, len(self.data.obs_names))
        if chunksize is None:
            normexpr = model.get_normalized_expression(
                adata=self.data,
                n_samples=n_samples,
                library_size=library_size,
                indices=None,
                **kwargs,
            )
            normexpr = dataframe_into_csc(normexpr)

        else:
            chunks = np.split(cell_indices, indices_or_sections=chunksize)
            normexpr = scipy.sparse.csc_array(
                (1, len(self.data.var_names)), dtype=np.float32
            )
            for chunk in chunks:
                normexpr_chunk = model.get_normalized_expression(
                    adata=self.data,
                    n_samples=n_samples,
                    library_size=library_size,
                    indices=chunk,
                    **kwargs,
                )
                normexpr = scipy.sparse.vstack(
                    normexpr, dataframe_into_csc(normexpr_chunk)
                )

        return normexpr[0:, :]

    ### Needs checking ###
    @_is_model_loaded
    def differential_expression(
        self,
        groupby: str,
        factorial: bool = False,
        comparisons: typing.Optional[list[tuple[str, str]]] = None,
        fdr_target=0.05,
    ) -> pd.DataFrame:
        if factorial and comparisons is not None:
            raise ValueError(
                "factorial and comparisons cannot be specified at the same time"
            )
        if factorial:
            comparisons = itertools.combinations(self.data.obs[groupby].unique(), r=2)
        else:
            if comparisons is not None:
                comparisons = comparisons

            else:
                comparisons = [(None, None)]

        alldeg = pd.DataFrame()
        for c in comparisons:
            matchup = f"{c[0]}__vs__{c[1]}"
            de = self.model.differential_expression(
                adata=self.data,
                groupby=groupby,
                group1=c[0],
                group2=c[1],
                fdr_target=fdr_target,
            )
            de["matchup"] = matchup
            alldeg = pd.concat([alldeg, de])

        return alldeg

    @staticmethod
    def load(fp, model):
        sm = ScviModeler()
        sm.read(fp)
        n_latent = int(
            re.search("[_]n[_]latent[_]([0-9]+)[_]", os.path.basename(model)).group(1)
        )
        sm.models[n_latent] = sm.type_.load(model, adata=sm.data)
        return sm

    def de_volcano_plot(csv_path, fdr_target=0.05, title=""):
        sample_name = os.path.basename(csv_path).split("_")[0]
        de = pd.read_csv(csv_path)
        de = de.set_index(de.columns[0])
        fdr_key = f"is_de_fdr_{fdr_target}"

        sp = (
            so.Plot(data=de, x="lfc_mean", y="proba_de", color=fdr_key)
            .add(so.Dots(pointsize=1), legend=False)
            .scale(color=["red", "black"])
            .label(x="Log2FoldChange", y=f"p(DE), FDR = {fdr_target}", title=title)
            .theme({"axes.facecolor": "w", "axes.edgecolor": "#000000"})
        )
        print(plt.legend())
        return copy.deepcopy(sp)
