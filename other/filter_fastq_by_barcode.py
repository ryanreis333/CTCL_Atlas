#!/usr/bin/env python3
import gzip
from Bio import SeqIO

def load_barcodes(filename):
    """
    Load barcodes from a file into a set.
    """
    with open(filename, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def filter_fastq(fastq_r1, fastq_r2, out_r1, out_r2, barcodes, barcode_length=16):
    """
    Filter paired-end FASTQ files based on whether the read in R1 starts with one of the provided barcodes.
    """
    with gzip.open(fastq_r1, "rt") as handle_r1, gzip.open(fastq_r2, "rt") as handle_r2, \
         open(out_r1, "w") as out_handle_r1, open(out_r2, "w") as out_handle_r2:

        records_r1 = SeqIO.parse(handle_r1, "fastq")
        records_r2 = SeqIO.parse(handle_r2, "fastq")
        count = 0
        for rec1, rec2 in zip(records_r1, records_r2):
            cell_barcode = str(rec1.seq)[:barcode_length]
            if cell_barcode in barcodes:
                SeqIO.write(rec1, out_handle_r1, "fastq")
                SeqIO.write(rec2, out_handle_r2, "fastq")
                count += 1
        print(f"Total reads extracted: {count}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Filter paired-end FASTQ files by cell barcode.")
    parser.add_argument("--r1", required=True, help="Path to Read 1 FASTQ.gz file")
    parser.add_argument("--r2", required=True, help="Path to Read 2 FASTQ.gz file")
    parser.add_argument("--barcodes", required=True, help="Text file with barcodes to keep")
    parser.add_argument("--out_r1", required=True, help="Output file for filtered Read 1 (FASTQ)")
    parser.add_argument("--out_r2", required=True, help="Output file for filtered Read 2 (FASTQ)")
    parser.add_argument("--barcode_length", type=int, default=16, help="Length of cell barcode in R1")
    args = parser.parse_args()

    barcode_set = load_barcodes(args.barcodes)
    filter_fastq(args.r1, args.r2, args.out_r1, args.out_r2, barcode_set, args.barcode_length)
