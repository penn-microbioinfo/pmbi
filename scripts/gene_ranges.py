import sys
import pysam
import json
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import gffutils
import os

class Chrom(object):
    def __init__(self, name):
        self.name=name
        self.ranges=[]
        self.curr_range=None

    def insert(self, start, end):
        if end < start:
            start, end = (end, start)

        if self.curr_range is not None and (start in self.ranges[self.curr_range] or end in self.ranges[self.curr_range]):
            if not start in self.ranges[self.curr_range]:
                self.curr_range = range(start,max(self.ranges[self.curr_range]))
            if not end in self.ranges[self.curr_range]:
                self.ranges[self.curr_range]=range(min(self.ranges[self.curr_range]), end)
        else:
            self.ranges.append(range(start, end+1))
            self.curr_range = len(self.ranges)-1
    
    # Dosen't consider from the end of the_last range to end of the chromosome
    def uncovered_slow(self):
        self.curr_range = 0
        ranges = []
        c = 0
        for pos in range(0,max(self.ranges[::-1][0])+1):
            if not pos in self.ranges[self.curr_range]:
                if len(ranges) == 0:
                    ranges.append(range(0,1))
                else:
                    if pos == max(ranges[c])+1:
                        ranges[c] = range(min(ranges[c]), pos+1)
                    else:
                        ranges.append(range(pos, pos+1))
                        c=c+1

            if pos == max(self.ranges[self.curr_range]):
                self.curr_range = self.curr_range+1

        return ranges

    def uncovered_fast(self):
        ranges=[]
        c=0
        for i in range(0,len(self.ranges)):
            if i == 0:
                r = range(0,min(self.ranges[i]))
                if len(r) > 1:
                    ranges.append(r)

            elif i > 0 and i < len(self.ranges)-1:
                r = range(max(self.ranges[i-1]), min(self.ranges[i]))
                if len(r) > 1:
                    ranges.append(range(max(self.ranges[i-1]), min(self.ranges[i])))
            else:
                pass
        return ranges

    def len_uncovered(self):
        return sum([len(x) for x in self.uncovered_fast()])

def coveragesCsvToHist(path):
    csv = pd.read_csv(path)
    csv_nonzero=csv[csv["count"] > 0.0]
    ax = sns.histplot(csv_nonzero, x = "count", bins = np.logspace(1, np.ceil(np.log10(max(csv_nonzero["count"])))))
    plt.yscale("log")
    plt.xscale("log")
    return ax

def transform(d):
    try:
        d['transcriptId'] = d["transcript_id"]
    except KeyError:
        pass
    return d 

def gtf_db_exists(gtfpath):
    if os.path.isfile(f"{gtfpath}.db"):
        return True
    else:
        return False

def gtf_load_db(gtfpath):
    if not gtf_db_exists(gtfpath):
       db = gffutils.create_db(gtfpath, dbfn=f"{gtfpath}.db", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
    else: 
        db = gffutils.FeatureDB(f"{gtfpath}.db", keep_order=True)
    return db

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", help = "Path to BAM file.")
    parser.add_argument("-g", "--gtf", help = "Path to GTF for the reference genome.")
    parser.add_argument("-s", "--sample", help = "Sample name for printing outputs.")
    args = parser.parse_args()

    # Load database, creating if needed
    db=gtf_load_db(args.gtf)
    print(json.loads(db.execute(f"SELECT attributes FROM features WHERE seqid == 'NC_000001.11'").fetchone()["attributes"])["gene_id"])
    
    chroms = []
    curr_chrom = None
    for transcript in db.features_of_type('transcript', order_by='seqid'):
        tattr = transcript.attributes
        if tattr["transcript_biotype"][0] in ["mRNA", "transcript", "primary_transcript"]:
            chrom = transcript.seqid
            start = int(transcript.start)
            end = int(transcript.stop)
            if curr_chrom is None or chrom != curr_chrom.name:
                if curr_chrom is not None:
                    chroms.append(curr_chrom)
                curr_chrom = Chrom(chrom)
            curr_chrom.insert(start, end)

    bam = pysam.AlignmentFile(args.bam)
    region_coverages = []
    features = []
    stats = bam.get_index_statistics()
    for chr in chroms:
        chr_stats = [s for s in stats if s.contig==chr.name][0]
        for (rid,ucr) in enumerate(chr.uncovered_fast()):
            record = {}
            record["chr_name"] = chr.name
            #print(chr.ranges[rid], ucr)
            region_str = f"{chr.name}:{min(ucr)+1}-{max(ucr)}"
            region_len = max(ucr)-min(ucr)
            #print(chr.ranges[rid], region_str,region_len)
            record["region_str"] = region_str
            record["region_len"] = region_len
            record["reads_mapped"] = bam.count(region=region_str)
            record["reads_mapped_by_region_size"] = np.true_divide(record["reads_mapped"], record["region_len"])
            record["total_reads"] = chr_stats.total
            if record["total_reads"] == 0:
                record["prop_total_reads"] = 0.0
            else:
                record["prop_total_reads"] = np.true_divide(record["reads_mapped"], record["total_reads"])
            for row in db.execute(f"SELECT id,start,end,strand,attributes FROM features WHERE seqid=='{chr.name}' AND start>={min(ucr)+1} AND end<={max(ucr)} AND featuretype=='gene'").fetchall():
                row_attr = json.loads(row["attributes"])
                feature_record={}
                feature_record["region_str"] = region_str
                feature_record["feature_region_str"] = f"{chr.name}:{row['start']}-{row['end']}"
                feature_record["feature_len"] = int(row["end"])-int(row["start"])
                feature_record["id"] = row["id"] 
                feature_record["stand"] = row["strand"]   
                for k in [x for x in ["gene_id", "description", "pseudo", "gene_biotype"] if x in row_attr]:
                    feature_record[k] = row_attr[k][0]
                feature_record["reads_mapped"] = bam.count(region=feature_record["feature_region_str"])
                feature_record["reads_mapped_by_feature_len"] = np.true_divide(feature_record["reads_mapped"], feature_record["feature_len"])
                features.append(feature_record)
                #print(feature_record)

            region_coverages.append(record)
        print(f"Done with contig: {chr.name}")

    frame = pd.DataFrame(region_coverages) 
    frame.to_csv(f"{args.sample}_coverages.csv", index = False)

    feature_frame = pd.DataFrame(features) 
    feature_frame.to_csv(f"{args.sample}_features.csv", index = False)
    
    """
    frame["prop_reads_by_len"] = frame["prop_total_reads"]/frame["region_len"]
    sums = frame[["chr_name", "prop_total_reads"]].groupby("chr_name").aggregate("sum")
    sums = frame[["chr_name", "prop_total_reads"]].groupby("chr_name").aggregate("sum")
    sums["prop_reads_by_len"] = sums
    print(sums)
    print(sums.sum())
    print(sumsp55)
    """








