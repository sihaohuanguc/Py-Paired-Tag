#!/usr/bin/env python
# -*- coding: utf-8

import gzip, os, sys
from .read_processor import read_processor
from itertools import islice

def combine3(file1, file2, out_prefix):
    total = 0
    passed = 0
    dna = 0
    rna = 0
    und = 0

    type_report_file = f"{out_prefix}_type_report.xls"      
    combined_file = f"{out_prefix}_combined.fq"              

    with open(type_report_file, 'w') as outType, open(combined_file, 'w') as outfile:
        with gzip.open(file1, 'rt') as red1, gzip.open(file2, 'rt') as red2:
            while True:
                lines_r1 = list(islice(red1, 4))
                lines_r2 = list(islice(red2, 4))

                if not lines_r1 or not lines_r2:
                    break
                
                total += 1
                title1, seq1, _, qual1 = lines_r1
                title2, seq2, _, qual2 = lines_r2
                title1 = title1.split(' ')[0]
                title2 = title2.split(' ')[0]
                seq1 = seq1.strip()
                seq2 = seq2.strip()
                qual1 = qual1.strip()
                qual2 = qual2.strip()

                if title1 != title2:
                    print(f"Title mismatch: {title1} vs {title2}")
                    continue

                read_2 = read_processor(seq2)  # Initialize read_2 with read 2 sequence
                read_2.trim()

                umi = read_2.umi
                seq2 = read_2.sbc1 + read_2.sbc2 + read_2.sbc4
                qual2 = qual2[:len(seq2)]

                if len(seq2) != 20:
                    continue

                read_type = read_2.type

                if read_type == "d":
                    dna += 1
                elif read_type == "r":
                    rna += 1
                elif read_type == "n":
                    und += 1
                
                updated_readname = f"{title1}:{umi}:{seq1}:{qual1}"
                outfile.write(f"{updated_readname}\n{seq2}\n+\n{qual2}\n")

                outType.write(f"{title1}\t{read_type}\n")
                passed += 1

    ratio = round((passed / total) * 100, 2) if total > 0 else 0.0
    print(f"==================================================")
    print(f"Paired-seq/Tag Barcode Locator Report: {out_prefix}")
    print(f"# Total raw reads: {total}")
    print(f"# Full barcoded reads: {passed}")
    print(f"# DNA reads: {dna}")
    print(f"# RNA reads: {rna}")
    print(f"# Undetermined reads: {und}")
    print(f"% Full barcode reads: {ratio}%")
    print(f"==================================================")

    os.system(f"gzip -f {out_prefix}_type_report.xls")
    os.system(f"gzip -f {out_prefix}_combined.fq")

if __name__ == "__main__":
    combine3(sys.argv[1], sys.argv[2], sys.argv[3]) 