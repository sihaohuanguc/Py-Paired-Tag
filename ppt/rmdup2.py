#!/usr/bin/env python
# -*- coding: utf-8

import subprocess, sys
from collections import defaultdict

def rmdup2(bam, out_prefix):
    dup_file = f"{out_prefix}_dup.xls"
    rmdup_bam = f"{out_prefix}_rmdup.bam"

    print("Reading BAM file...")
    cmd = f"samtools view {bam}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    hash_map = defaultdict(lambda: defaultdict(dict))
    dupcounts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for line in process.stdout:
        fields = line.decode('utf-8').strip().split('\t')
        chr_name = fields[2]
        pos = int(fields[3])
        umi = fields[0][-19:]  # this include both CB ID and the UMI   # change this if format is different !!!
        
        hash_map[chr_name][pos][umi] = True
        dupcounts[chr_name][pos][umi] += 1

    process.stdout.close()
    process.wait()

    print("Writing duplicates to file...")
    with open(dup_file, 'w') as f:
        for chr_name, pos_dict in dupcounts.items():
            for pos, umi_dict in pos_dict.items():
                for umi, count in umi_dict.items():
                    f.write(f"{chr_name}\t{pos}\t{umi}\t{count}\n")

    # now remove duplicates from the BAM file
    print("Removing duplicates...")
    cmd = f"samtools view -h {bam}"  # -h with head line
    in_process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    cmd_out = f"samtools view -b - -o {rmdup_bam}"
    out_process = subprocess.Popen(cmd_out, shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


    for line in in_process.stdout:
        line = line.decode('utf-8')
        if line.startswith('@'):
            out_process.stdin.write(line)
        else:
            fields = line.strip().split('\t')
            chr_name = fields[2]
            pos = int(fields[3])
            umi = fields[0][-19:]        # change this if format is different !!!
            
            if not hash_map[chr_name][pos][umi]:  # it's not the first time you see this UMI
                continue
            out_process.stdin.write(line)   # it is the first time you see this UMI
            hash_map[chr_name][pos][umi] = False  # after you record the UMI for the first time, mark it as False

    in_process.stdout.close()
    out_process.stdin.close()
    in_process.wait()
    out_process.wait()

    print("Duplicate removal completed.")

if __name__ == "__main__":
    rmdup2(sys.argv[1], sys.argv[2])
