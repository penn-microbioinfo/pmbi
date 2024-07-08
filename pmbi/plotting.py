import matplotlib as mpl
import matplotlib.axes
import matplotlib.axis
import matplotlib.pyplot as plt
import numpy as np
import gc
from functools import wraps

def _check_output_prefix(inner):
    @wraps(inner)
    def wrapper(self, *args, **kwargs):
        if self.output_prefix is None:
            raise ValueError("No output_prefix set for Paneler despite requiring advancing image.")
        else:
            return inner(self, *args, **kwargs)

    return wrapper
    

class Paneler(object):
    def __init__(
        self,
        nrow,
        ncol,
        output_prefix=None, # Only meaningful if intending to save anything to file
        figsize=(3, 3),
        layout="tight",
        format="tiff",
        dpi=400,
    ):
        Paneler.theme()
        self.nrow = nrow
        self.ncol = ncol
        self.output_prefix = output_prefix
        self.figsize = figsize
        self.layout = layout
        self.format = format
        self.dpi = dpi
        self.fig, self.axs = self._new_fig()
        self.panel_idx = 0
        self.image_num = 1
        self.current_ax = self._get_ax((0, 0))

    def _new_fig(self):
        return plt.subplots(
            nrows=self.nrow, ncols=self.ncol, figsize=self.figsize, layout=self.layout
        )

    @staticmethod
    def theme():
        plt.rcdefaults()
        plt.rcParams.keys()
        plt.rcParams.update(
            {  # "font.sans-serif": "Arial",
                "font.size": 5,
                "figure.dpi": 300,
                "axes.titlesize": 6,
            }
        )
        plt.rcParams["legend.frameon"] = False

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

    def _savefig(self, filename: str):
        self.subplots_adjust(hspace=0.5, wspace=0.5)
        self.fig.savefig(filename, bbox_inches="tight")
        plt.close(self.fig)
        plt.clf()
        gc.collect()

    @_check_output_prefix
    def _advance_image(self):
        self._savefig(filename = f"{self.output_prefix}_{self.image_num}.{self.format}")
        self.fig, self.axs = self._new_fig()
        self.image_num += 1
        self.panel_idx = 0

    @_check_output_prefix
    def close(self):
        if self.image_num == 1:
            self._savefig(filename = f"{self.output_prefix}.{self.format}")
        else:
            self._savefig(filename = f"{self.output_prefix}_{self.image_num}.{self.format}")
