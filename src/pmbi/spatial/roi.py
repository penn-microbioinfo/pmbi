from PIL import Image, ImageDraw
from pathlib import Path
import pandas as pd
import numpy as np

# Allow big images, because the scans are big
Image.MAX_IMAGE_PIXELS = 1500000000

META=Path("/home/ubuntu/projmnt/tcsl/metadata/BreastVax_human/")
IMG=Path("/home/ubuntu/projmnt/tcsl/images/BreastVax_human")

# %%
# metadata = pd.read_excel(META.joinpath("Breastvax_geomx_metadata_table_JT_AR.xlsx"))
metadata = pd.read_csv(META.joinpath("Breastvax_geomx_phenoData_for_JT_completed_FULL.tsv"), sep = "\t")

# %%
scan = Image.open(IMG.joinpath("2023.10.11_PC23_JUNE772_1_12-12-22_A.jpeg"))
scan.size
metadata.columns

metadata["scan_name"].unique()
d = ImageDraw.Draw(scan)
for _idx,roi in metadata.iterrows():
    if roi["scan_name"] == '2023.10.11_pc23_june772#1_12-12-22_a':
        x,y = (roi["roi_x"], roi["roi_y"])
        radius = np.sqrt(np.divide(roi["area"], np.pi))
        print(x,y,radius)
        d.circle((x,y), radius=radius, fill = (255,255,255))

d.circle((10000,20000), radius=5000, fill = (255,255,255))
d
# %%
scan.save("/srv/http/breastvax/slide02_BX-006_test.jpeg", dpi=(95.986,95.986))
# %%

