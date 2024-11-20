import logging

import anndata


def counts(adata: anndata.AnnData, layer: str | None = None, copy: bool = False):
    if layer is None:
        counts = adata.X
    else:
        if layer in adata.layers:
            counts = adata.layers[layer]
        else:
            if layer == "raw":
                logging.warning( f"Trying to pull counts from AnnData.raw because of provided layer name: {layer}")
                if adata.raw is not None:
                    counts = adata.raw.X
                else:
                    raise ValueError(f"Layer dose not exist in AnnData: {layer}")
            else:
                raise ValueError(f"Layer dose not exist in AnnData: {layer}")

    if copy:
        return counts.copy()
    else:
        return counts
