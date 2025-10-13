import numpy as np
import pandas as pd

# %% FUNC: Class to handle printing tables to markdown and printing to multiple slides
class MultiSlideTable(object):
    def __init__(self, table: pd.DataFrame, rows_per_page=25):
        self.table = table
        self.rows_per_page=25
        self.slide_title = "Regions of Interest (ROIs)"
        self.nrows = self.table.shape[0]
        self.splits = self._splits()
    def _splits(self):
        splits = []
        for start in range(0, self.nrows, self.rows_per_page):
            if start + self.rows_per_page > self.nrows:
                stop = self.nrows
            else:
                stop = start + self.rows_per_page
            splits.append((start, stop))
        return splits
    def to_markdown(self):
        out = ""
        out += "```{python}\n"
        out += "from IPython.display import Markdown, display\n"
        out += "from tabulate import tabulate\n\n"
        out += "```\n"
        i=0
        for start,stop in self._splits():
            if i==0:
                out+="# Regions of Interest (ROIs)\n\n"
            else:
                out+="# Regions of Interest (ROIs) cont.\n\n"
            # out+="```{python}\n"
            out+="\\fontsize{7}{7}\\selectfont\n"
            # out+="```{=html}\n"
            table_chunk=self.table.iloc[start:stop,:]
            out+=f"{table_chunk.to_markdown(index=False)}\n\n"  
            # out+=(f'table_md="""\n{table_chunk.to_markdown()}\n"""\n\n')
            # out+="display(Markdown(table_md))\n\n"
            # out+="```\n\n"
            i=i+1
        return out
# %%

