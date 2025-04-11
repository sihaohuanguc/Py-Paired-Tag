#!/usr/bin/env python
# -*- coding: utf-8

import sys, os

def convert2(in_name, out_prefix):
    total = 0
    passed = 0

    # Read the input file directly
    with open(in_name, 'r') as inbam, open(out_prefix + '_cov.fq', 'w') as fout:
        for line in inbam:
            line = line.strip()
            if line.startswith('@'):
                continue

            total += 1
            fields = line.split('\t')
            readname = fields[0]
            chr_name = fields[2]   # this is the CB ID

            if chr_name == '*':
                continue

            tmp = readname.split(':')   # get the read 1 seq and score
            new_readname = f"@{':'.join(tmp[:7])}:{chr_name}:{tmp[7]}"
            seq = tmp[8]
            qual = readname[-len(tmp[8]):]  # get a string of the same length as seq from backward, which is the quality score

            fout.write(f"{new_readname}\n{seq}\n+\n{qual}\n")

            passed += 1
    
    os.system(f"gzip -f {out_prefix}_cov.fq")  # Compress the output file

    print(f"{total} reads processed.")
    print(f"{passed} mapped reads.")
    print(f"{total - passed} unmapped reads.")
    print(f"{passed / total:.2%} mapped reads.")

if __name__ == "__main__":
    convert2(sys.argv[1], sys.argv[2])
