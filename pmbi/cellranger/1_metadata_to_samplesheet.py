import pandas as pd
import argparse

### Column names
donorid = "DonorID"
modality = "Library_Type"
chip_well = "Chip_Well"
index_plate = "Index_Plate"
well = "Well"

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--metadata", action = "store", required = True, help = "Path to metadata sheet.")
parser.add_argument("-o", "--output", action = "store", required = True, help = "Output filename.")
parser.add_argument("-d", "--donors", action = "store", required = False, help = "Donors to include in output. Default is to include all donors.")
args = parser.parse_args()

md = pd.read_csv(args.metadata, dtype = str)
if args.donors is not None:
    donors = args.donors.strip().split(',')
    md = md[md[donorid].isin(donors)] 

md["Lane"] = "*"
md["Sample"] = md[donorid] + "_" + md[modality] + "_" + md[chip_well]
md["Index"] = "SI" + "-" + md[index_plate] + "-" + md[well]
ss = md[["Lane", "Sample", "Index"]]

ss.to_csv(args.output, index = False, header = True)
