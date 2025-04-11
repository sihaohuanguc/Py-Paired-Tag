#!/usr/bin/env python
# -*- coding: utf-8

import re

def extract_gene_info(gtf_file, output_file):
    with open(gtf_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            if line.startswith('#'):
                continue
            columns = line.strip().split('\t')
            if columns[2] != 'gene':
                continue
            chrom = columns[0]
            start = columns[3]
            end = columns[4]
            
            gene_id_match = re.search(r'gene_id "(.*?)"', columns[8])
            gene_name_match = re.search(r'gene_name "(.*?)"', columns[8])
            gene_name = gene_name_match.group(1) if gene_name_match else 'NA'
            gene_id = gene_id_match.group(1).split('.')[0] if gene_id_match else 'NA'
            
            out.write(f"{chrom}\t{start}\t{end}\t{gene_id} {gene_name}\n")  # the space between gene_id and gene_name is important
    print("Gene information extracted successfully.")

if __name__ == "__main__":
    extract_gene_info("gencode.v47.annotation.gtf","hg38.big_with_chr.txt")
