from typing import Union
import numpy as np
import pandas as pd
from pmbi.geometry import Curve
from sklearn.neighbors import KernelDensity


class Kde:
    def __init__(self, X: Union[pd.Series, np.ndarray, list] , x_len: int = 1000, backend:str = "sklearn", **kwargs):
        if isinstance(X, pd.Series):
            self.X = X.to_numpy()[:,np.newaxis]
        elif isinstance(X, np.ndarray):
            pass
        elif isinstance(X, list):
            self.X = np.array(X)[:,np.newaxis]
        else:
            raise ValueError("invalid type of X")

        if backend=="sklearn":
            x_pts = np.linspace(min(self.X), max(self.X), x_len)
            self.kde = KernelDensity(kernel="gaussian", **kwargs).fit(self.X)
            y_pts = np.exp(self.kde.score_samples(x_pts))
            self.curve = Curve.from_xy(x_pts, y_pts)
        else:
            raise NotImplementedError
