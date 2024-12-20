from PIL import Image, ImageDraw, ImageEnhance, ImagePalette, ImageOps
from pathlib import Path
import pandas as pd
import numpy as np
import tifffile
import imagecodecs
import matplotlib.pyplot as plt

# Allow big images, because the scans are big
Image.MAX_IMAGE_PIXELS = 5000000000

META=Path("/home/ubuntu/projmnt/tcsl/metadata/BreastVax_human/")
IMG=Path("/home/ubuntu/projmnt/tcsl/images/BreastVax_human")

# %%
# metadata = pd.read_excel(META.joinpath("Breastvax_geomx_metadata_table_JT_AR.xlsx"))
metadata = pd.read_csv(META.joinpath("Breastvax_geomx_phenoData_for_JT_completed_FULL.tsv"), sep = "\t")

# %%
scan = tifffile.imread(IMG.joinpath("2023.10.11_PC23_JUNE772_1_12-12-22_A.ome.tiff"), is_ome = True, level=0)
# fig, axes = plt.subplots(ncols=1,nrows=1,squeeze=False,figsize=(12,3))
# for i in range(0,1):
#     axes[0,i].imshow(scan[i])
# plt.savefig("/srv/http/breastvax/slide02_BX-006_test_crop.png", dpi=2000)
# %%

scan.shape
roi.keys()
im_rgb.size
roi.scan_width*2

tiff_dim = np.array([scan.shape[2], scan.shape[1]])
tiff_dim
scan_dim = np.array([roi.scan_width,roi.scan_height])
scan_dim
print(tiff_dim, scan_dim)
scan_offset = np.array([roi.scan_offset_x, roi.scan_offset_y])
scan_corr = scan_dim+scan_offset
ratio_wOff = np.divide(tiff_dim, (scan_corr) )
ratio_noOff = np.divide(tiff_dim, scan_dim)
scan_dim*ratio_noOff
metadata["roi_x"].max()
metadata["roi_y"].max()



# %%
def correct_roi_coords(x, y, area, resolution_level, x_offset=0, y_offset=0):
    aligned = np.array([x+x_offset, y+y_offset])
    xy_scaled = aligned*(1/(2**resolution_level))
    # %% BUG: Aread scaling still isn't right
    area_scaled = area*(1/(2**resolution_level))
    return tuple([round(x) for x in [xy_scaled[0],xy_scaled[1], area_scaled]])

# %%
def i16_i8_arr(i16: np.array):
    ret = (i16-0+1)/(255-0+1)
    return ret.round()

# %%
def i16_i8(i16: np.ndarray|int):
    ret = (i16-0+1)/(255-0+1)
    if isinstance(ret, np.ndarray):
        return np.round(ret).astype(np.uint8)
    elif isinstance(ret, float):
        return np.uint8(np.round(ret))
    else:
        raise ValueError(f"Expected int or np.ndarray.")

# %%
dsp_channels = [
        ( "blue", (0,0,255), 113, 1292 ),
        ( "green", (0,255,0), 146, 572 ),
        ( "yellow", (255,255,0), 143, 327 ),
        ( "red", (255,0,0), 364, 766 ),
        ]
# %%
reslev = 0
scan_name='2023.10.11_slide4_pc23_june772_tma1'
scan_d = tifffile.imread(IMG.joinpath("2023.10.11_Slide4_PC23_JUNE772_TMA1.ome.tiff"), is_ome = True, level=reslev)
for cidx, cintensity in enumerate(dsp_channels):
    im = Image.fromarray(i16_i8(scan_d[cidx]), mode="L")
    cname, color, lower, upper = cintensity
    print(cname)
    im = ImageOps.colorize(im, black=(0,0,0), white=color, blackpoint=i16_i8(lower), whitepoint=i16_i8(upper))
    d = ImageDraw.Draw(im)
    for _idx,roi in metadata.iterrows():
        if roi["scan_name"] == scan_name:
            x,y,area = correct_roi_coords(x=roi["roi_y"], y=roi["roi_x"], area=roi["area"], resolution_level=reslev)
            r = np.sqrt(np.divide(area, np.pi))
            # print(x,y,r)
            d.circle((x,y), radius=r, fill = (0,0,0)) #0 intensity
    im.save(f"/srv/http/breastvax/slide_4_down_{cname}.png")

# %%

# %%
# # d.circle(s-(1000,1000), radius=500, fill = (138,43,226)) #purple
# # d.circle(s-(s-1000), radius=500, fill = (255,99,71)) #red
# d.polygon( xy=(((0+1500)*0.25,(0+3267)*0.25), ((0+1500)*0.25,(42617+3267)*0.25), ((29767+1500)*0.25,(42617+3267)*0.25), ((29767+1500)*0.25,(0+3267)*0.25)), outline=(138,43,226), width = 20) 
# offsets
# scan_dimens
# tiff_dim[0]-1500-29767
# tiff_dim[1]-3267-42618
# scan_dim
# %%
scan.shape
im = Image.fromarray(scan[0])
e = ImageEnhance.Brightness(im)
im.mode
im_rgb = im.convert("RGB")
im.size
e.enhance(1.0)
im_rgb.save("/srv/http/breastvax/slide02_BX-006_test_full.png")

im.size

metadata["scan_name"].unique()
d = ImageDraw.Draw(scan)
for _idx,roi in metadata.iterrows():
    if roi["scan_name"] == '2023.10.11_pc23_june772#1_12-12-22_a':
        x,y = (roi["roi_x"], roi["roi_y"])
        radius = np.sqrt(np.divide(roi["area"], np.pi))
        print(x,y,radius)
        d.circle((x,y), radius=radius, fill = (1))

# d.circle((10000,20000), radius=5000, fill = (255,255,255))
d
# %%
scan.save("/srv/http/breastvax/slide02_BX-006_test.tiff", dpi=(50,50))
# %%

