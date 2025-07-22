import argparse
from pathlib import Path
from pmbi.illumina.demux import fastq_missing_for

def main():
    parser = argparse.ArgumentParser(description="Check for missing FASTQ files based on an Illumina sample sheet")
    parser.add_argument("sample_sheet", type=str, help="Path to Illumina sample sheet")
    parser.add_argument("fastq_dir", type=str, help="Path to directory containing FASTQ files")
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    
    args = parser.parse_args()
    
    # Convert paths to Path objects
    sample_sheet_path = Path(args.sample_sheet)
    fastq_dir_path = Path(args.fastq_dir)
    
    # Check that paths exist
    if not sample_sheet_path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {sample_sheet_path}")
    
    if not fastq_dir_path.exists():
        raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir_path}")
    
    # Open sample sheet and check for missing files
    with open(sample_sheet_path, 'r') as sample_sheet:
        _missing_samples = fastq_missing_for(fastq_dir_path, sample_sheet)

if __name__ == "__main__":
    main()

