import matplotlib

def dict_try_insert(d: dict, key, value):
    if not key in d:
        d[key] = value
    return d

def figListToPdf(figs: list, filename: str):
    pdf = matplotlib.backends.backend_pdf.PdfPages(filename)
    for fig in figs:
        pdf.savefig(fig, bbox_inches = "tight")
    pdf.close()
    return None
