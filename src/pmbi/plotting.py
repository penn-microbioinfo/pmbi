# %%
import copy
import gc
import os
import re
from functools import wraps
from itertools import count
from typing import Iterable

import matplotlib as mpl
import matplotlib.axes
import matplotlib.axis
import matplotlib.pyplot as plt
import munch
import numpy as np
import palettable
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.colorizer import Colorizer
from matplotlib.colors import BoundaryNorm, ListedColormap, hex2color, rgb2hex
from matplotlib.markers import MarkerStyle
from matplotlib.lines import Line2D


# %%
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
class MosaicPaneler:
    def __init__(
        self,
        design=[["a"]],
        figsize=(3,3),
        format="png",
        layout="tight",
        dpi=400,
        **kwargs
    ):
        Paneler.default_theme()
        self._figure = Figure(figsize=figsize, dpi=dpi, layout="tight")
        self.axes_dict = self._figure.subplot_mosaic(mosaic=design)
        self.axes_list = self._figure.axes
        self.aidx = 0

    def next_ax(self):
        if self.aidx == len(self.axes_list):
            raise IndexError("Axes exhausted")
        ax = self.axes_list[self.aidx]
        self.aidx+=1
        return ax

    def current_ax(self):
        return self.axes_list[self.aidx-1]

    def get_ax(self, key: str):
        return self.axes_dict[key]

    

# %%
def axvlines(ax, xintercepts, **kwargs):
    for x in xintercepts:
        ax.axvline(x=x, **kwargs)


def axvspans(ax, xmins, xmaxes, **kwargs):
    for xmin, xmax in zip(xmins, xmaxes):
        ax.axvspan(xmin=xmin, xmax=xmax, **kwargs)


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
    def __init__(self, theme_specification=THEME_SPEC, **kwargs):
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

    # def __getattr__(
    #     self, name
    # ):  # Only called if default attribute access fails with AttributeError
    #     return getattr(self.inner, name)

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
        ax.set_xlabel(
            self.inner.axislabels.text.x, fontsize=self.inner.axislabels.fontsize
        )
        ax.set_ylabel(
            self.inner.axislabels.text.y, fontsize=self.inner.axislabels.fontsize
        )
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
    aspect="auto",
    cmap=palettable.scientific.sequential.Hawaii_20_r.mpl_colormap,
    theme=Theme(),
):
    im = ax.imshow(matrix, cmap=cmap, aspect=aspect)
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
    ax.set_xlabel(xlab, fontsize=theme.inner.axislabels.fontsize)
    ax.set_ylabel(ylab, fontsize=theme.inner.axislabels.fontsize)
    ax.set_xticks(range(0, ncol))
    ax.set_xticklabels(matrix.columns, fontsize=theme.inner.ticklabels.fontsize)
    ax.set_yticks(range(0, nrow))
    ax.set_yticklabels(matrix.index, fontsize=theme.inner.ticklabels.fontsize)
    return im


# %%
def scatter(
    x, y, ax, marker="o", facecolors="none", edgecolors="none", linewidths=0, **kwargs
):
    ax.scatter(
        x=x,
        y=y,
        marker=marker,
        facecolors=facecolors,
        edgecolors=edgecolors,
        linewidths=linewidths,
        **kwargs,
    )


# %%
class CategoricalColormap:
    def __init__(self, values: Iterable[str], colors: dict[str, str]):
        self._values = np.array(values)
        if not all([x in colors.keys() for x in self._values]):
            raise ValueError(
                "Some values do not have a corresponding color in supplied dict"
            )
        if not all([isinstance(x,str) for x in colors.values()]):
            self._hex_colors = {k:rgb2hex(c) for k,c in colors.items()}
            self._cmap = ListedColormap(list(colors.values()))
        else:
            if not all([re.match("[#][A-Za-z0-9]{6}", x) for x in colors.values()]):
                raise ValueError("str colors do not look like hex colors")
            self._hex_colors = colors
            self._cmap = ListedColormap([hex2color(c) for c in self._hex_colors.values()])

        self._counter = count()
        self._mapping = self._set_mapping()

    def _set_mapping(self) -> dict[str, int]:
        mapping = {}
        for c in self._hex_colors:
            if c not in mapping:
                mapping[c] = next(self._counter)
        return mapping

    def encoded(self):
        return np.array([self._mapping[c] for c in self._values])

    def cmap(self):
        return self._cmap

    def norm(self):
        boundaries = np.arange(-0.5, len(self._mapping))
        return BoundaryNorm(boundaries, self._cmap.N)

    def mpl_color_dict(self):
        mpl_colors = self._cmap.colors
        labels = self._mapping.keys()
        return {k: v for k, v in zip(labels, mpl_colors)}

    # def to_colorizer(self):
    #     """Return a Colorizer instance (matplotlib >= 3.10)."""
    #     cmap, norm = self.to_cmap_norm()
    #     return Colorizer(cmap=cmap, norm=norm)


class Legend:
    def __init__(self, handles: list):
        self.handles = handles

    @staticmethod
    def from_mpl_color_dict(d, marker=".", markersize=12, linestyle="none", **kwargs):
        handles=[]
        for label, color in d.items():
            artist=Line2D(
                xdata=[],
                ydata=[],
                marker=marker,
                markersize=markersize,
                linestyle=linestyle,
                label=label,
                color=color,
                **kwargs
            )
            handles.append(artist)

        return Legend(handles=handles)

    def on_ax(self, ax):
        ax.legend(handles=self.handles)


