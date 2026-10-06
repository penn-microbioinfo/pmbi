class CellrangerConfigTemplate:
    def __init__(self):
        self.string = """
### Illumina demultiplexing -------
[metadata.colnames]
sample_id = "Experiment_Name"
# donor_id = "DonorID"
modality = "Library_Type"
# chip_well = "Chip_Well"
# index_plate = "Index_Plate"
# well = "Well"
library_name = "Library_Name"
# index_name = "10X_Index_Name"
### -------

### Modalities -------
[[modalities]]
name = "GEX"
reference = "" # cellranger-formatted GEX reference (ie from `cellranger-mkref`
create-bam = true
no-secondary = true

[[modalities]]
name = "ADT"
reference = "" # in cellranger Antibody Capture reference format

[[modalities]]
name = "HTO"
reference = "" # in cellranger Antibody Capture reference format
# ## -------

[filename_patterns]
sample = "^(well[_][0-9]+)[_][A-Z]+[_]S[0-9]+(?:[_]L[0-9]{3})?_R[0-9][_][0-9]{3}[.]fastq[.]gz$"
modality = "^well[_][0-9]+[_]([A-Z]+)[_]S[0-9]+(?:[_]L[0-9]{3})?_R[0-9][_][0-9]{3}[.]fastq[.]gz$"
illumina_bits = "^well[_][0-9]+[_]A-Z]+[_](S[0-9]+(?:[_]L[0-9]{3})?_R[0-9][_][0-9]{3})[.]fastq[.]gz$"
read_number = "[_](R[0-9])[_]"
include = "[.]fastq[.]gz$"

[run]
wd = "/home/amsesk/super1/montse/cr_runs"
cellranger_flavor = "multi"
n_jobs = 5

[command_line_args]
localcores = 8
localmem = 32
"""
    def __str__(self):
        return self.string

# %%
