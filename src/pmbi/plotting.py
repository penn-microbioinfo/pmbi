import gc
import os
import copy
from functools import wraps

import matplotlib as mpl
import matplotlib.axes
import matplotlib.axis
import matplotlib.pyplot as plt
import munch
import numpy as np
import palettable
import pandas as pd


def _check_output_prefix(inner):
    @wraps(inner)
    def wrapper(self, *args, **kwargs):
        if self.output_prefix is None:
            raise ValueError(
                "No output_prefix set for Paneler despite requiring advancing image."
            )
        else:
            return inner(self, *args, **kwargs)

    return wrapper


# %%
class Paneler(object):
    def __init__(
        self,
        nrow,
        ncol,
        figsize=(3, 3),
        format="tiff",
        output_prefix=None,  # Only meaningful if intending to save anything to file
        layout="tight",
        dpi=400,
        **kwargs,
    ):
        Paneler.default_theme()
        self.nrow = nrow
        self.ncol = ncol
        self.output_prefix = output_prefix
        self.figsize = figsize
        self.layout = layout
        self.format = format
        self.dpi = dpi
        self.subplots_kwargs = kwargs
        self.fig, self.axs = self._new_fig()
        self.panel_idx = 0
        self.image_num = 1
        self.current_ax = self._get_ax((0, 0))

    def _new_fig(self):
        return plt.subplots(
            nrows=self.nrow,
            ncols=self.ncol,
            figsize=self.figsize,
            layout=self.layout,
            **self.subplots_kwargs,
        )

    @staticmethod
    def default_theme():
        plt.rcdefaults()
        plt.rcParams.update(
            {
                # "font.sans-serif": "Arial",
                "font.size": 5,
                "figure.dpi": 300,
                "axes.titlesize": 6,
            }
        )
        plt.rcParams["legend.frameon"] = False

    @staticmethod
    def update_theme(update_dict):
        plt.rcParams.update(update_dict)

    @staticmethod
    def subplot_idx_to_pos(nrow, ncol, idx):
        r = int(np.floor(np.true_divide(idx, ncol)))
        rem = idx % ncol
        c = rem
        if r >= nrow:
            r = r % nrow
        return (r, c)

    @staticmethod
    def pos_to_subplot_idx(nrow: int, ncol: int, coords: tuple[int, int]) -> int:
        row = coords[0]
        col = coords[1]
        if row > nrow - 1 or col > ncol - 1:
            raise IndexError(
                f"coordinates {coords} outside of limits of panel with shape ({nrow}, {ncol})"
            )
        idx = (row * ncol) + (col)
        return idx

    def _get_ax(self, coords: tuple[int, int]):
        if self.panel_idx_out_of_bounds():
            raise IndexError("paneler.panel_idx >= panel_dimensions")
        else:
            if self.nrow * self.ncol == 1:
                ax = self.axs
            elif self.nrow == 1 or self.ncol == 1:
                if self.nrow == 1:
                    ax = self.axs[coords[1]]
                else:
                    ax = self.axs[coords[0]]
            else:
                r, c = coords
                ax = self.axs[r, c]
            return ax

    def panel_idx_out_of_bounds(self):
        if self.panel_idx > (self.nrow * self.ncol) - 1:
            return True
        else:
            return False

    def next_ax(self) -> matplotlib.axes.Axes:
        if self.panel_idx_out_of_bounds():
            self._advance_image()
        coords = Paneler.subplot_idx_to_pos(self.nrow, self.ncol, self.panel_idx)
        print(self.image_num, coords)
        self.current_ax = self._get_ax(coords)
        self.panel_idx += 1
        return self.current_ax

    """
    def insert(self, plt_function, **kwargs):
        if self.nrow*self.ncol == 1:
            plt_function(ax=self.ax, **kwargs)
        elif self.nrow or self.ncol == 1:
            plt_function(ax=self.ax[panel_idx], **kwargs)
        else:
            row,col = Paneler.subplot_idx_to_pos(self.nrow, self.ncol, self.panel_idx)
            plt_function(ax=self.ax[row,col], **kwargs)
    """

    def subplots_adjust(self, **kwargs):
        self.fig.subplots_adjust(**kwargs)

    def _savefig(self, filename: os.PathLike):
        self.subplots_adjust(hspace=0.5, wspace=0.5)
        self.fig.savefig(filename, bbox_inches="tight")
        plt.close(self.fig)
        plt.clf()
        gc.collect()

    @_check_output_prefix
    def _advance_image(self):
        self._savefig(filename=f"{self.output_prefix}_{self.image_num}.{self.format}")
        self.fig, self.axs = self._new_fig()
        self.image_num += 1
        self.panel_idx = 0

    @_check_output_prefix
    def close(self):
        if self.image_num == 1:
            self._savefig(filename=f"{self.output_prefix}.{self.format}")
        else:
            self._savefig(
                filename=f"{self.output_prefix}_{self.image_num}.{self.format}"
            )

# %%
THEME_SPEC = {
    "axislabels": {
        "fontsize": 10.0,
        "text": {
            "x": "x-axis",
            "y": "y-axis",
            },
        },
    "ticklabels": {"fontsize": 8.0},
    "title": {
        "text": "",
        "fontdict": {
            "fontsize": 10.0,
            },
        },
}

class Theme(object):
    def __init__(
        self,
        theme_specification = THEME_SPEC,
        **kwargs
    ):
        theme_specification = copy.deepcopy(theme_specification)
        theme_specification.update(theme_specification)

        THEME_SPEC_KWARGS_MAP = {
                "xlab": "axislabels.text.x",
                "ylab": "axislabels.text.y",
                "title": "title",
                }

        self.inner = munch.Munch.fromDict(theme_specification)
        
        # for kw,arg in kwargs.items():
        #     attr_path = THEME_SPEC_KWARGS_MAP[kw].split(".")
        #     print(kw, arg, attr_path)
        #     curr = self.inner
        #     for i,attr in enumerate(attr_path): 
        #         print(attr)
        #         if attr in dir(curr):
        #             if i < len(attr_path):
        #                 curr = getattr(curr, attr)
        #             else:
        #                 setattr(curr, attr, arg)
        #         else:
        #             raise AttributeError

    def __getattr__(
        self, name
    ):  # Only called if default attribute access fails with AttributeError
        return getattr(self.inner, name)

    def xlab(self, label):
        self.inner.axislabels.text.x = label
        return self

    def ylab(self, label):
        self.inner.axislabels.text.y = label
        return self

    def title(self, label):
        self.inner.title.text = label
        return self

    def apply_to(self, ax):
        ax.set_title(self.inner.title.text, **self.inner.title.fontdict)
        # xlabel, ylabel = (ax.get_xlabel(), ax.get_ylabel())
        ax.set_xlabel(self.inner.axislabels.text.x, fontsize=self.inner.axislabels.fontsize)
        ax.set_ylabel(self.inner.axislabels.text.y, fontsize=self.inner.axislabels.fontsize)
        xticklabels, yticklabels = (ax.get_xticklabels(), ax.get_yticklabels())
        ax.set_xticklabels(xticklabels, fontsize=self.inner.ticklabels.fontsize)
        ax.set_yticklabels(yticklabels, fontsize=self.inner.ticklabels.fontsize)
        return None


# %%
def heatmap(
    matrix,
    ax,
    xlab,
    ylab,
    annot=None,
    cmap=palettable.scientific.sequential.Hawaii_20_r.mpl_colormap,
    theme=Theme(),
):
    ax.imshow(matrix, cmap=cmap)
    if annot is not None:
        cell_text = pd.DataFrame(
            np.array(
                list(np.ndindex(annot.shape)),
                dtype=np.dtype([("y", float), ("x", float)]),
            )
        )
        cell_text["s"] = np.array(annot.to_numpy().ravel(), dtype=np.dtype(str))
        for r in cell_text.itertuples():
            ax.text(
                x=r.x,
                y=r.y,
                s=r.s,
                fontsize=12,
                horizontalalignment="center",
                verticalalignment="center",
            )
    nrow, ncol = matrix.shape
    ax.set_xlabel(xlab, fontsize=theme.axislabels.fontsize)
    ax.set_ylabel(ylab, fontsize=theme.axislabels.fontsize)
    ax.set_xticks(range(0, ncol))
    ax.set_xticklabels(matrix.columns, fontsize=theme.ticklabels.fontsize)
    ax.set_yticks(range(0, nrow))
    ax.set_yticklabels(matrix.index, fontsize=theme.ticklabels.fontsize)
    return None
