#!/usr/bin/env python
# -*- coding: utf-8

import subprocess, os, sys

def bam2Mtx2(ref, bam, out_prefix):
    # process ref file
    try:
        with open(ref, 'r') as f:
            pass
    except FileNotFoundError:
        print("Cannot find ref file.")
        exit(-1)

    genes = {}

    with open(ref, 'r') as f:
        for line in f:
            line = line.strip()
            tmp = line.split('\t')

            chr_name = tmp[0]
            pss = int(tmp[1])
            pse = int(tmp[2])
            gene_id = tmp[3] # this will read in both gene name and ID

            if chr_name not in genes:
                genes[chr_name] = {}
            
            for pos in range(pss // 1000, pse // 1000 + 1):  
                gene_line = f"{tmp[1]}\t{tmp[2]}\t{gene_id}"   # start, end , gene name and ID, use as a key for a dict
                if pos not in genes[chr_name]:
                    genes[chr_name][pos] = {}
                genes[chr_name][pos][gene_line] = 1  # remember, each pos can have multiple genes

    print("Reading ref annotation...")
    # print(genes)

    # check if the bam file exist
    try:
        with open(bam, 'r') as f:
            pass
    except FileNotFoundError:
        print("Cannot find bam file.")
        exit(-1)

    cmd = f"samtools view {bam}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    hash_data = {}
    genelist = {}
    celllist = {}

    print("Processing sam/bam file...")

    for line in process.stdout:
        line = line.decode('utf-8').strip()
        # print(line)
        fields = line.split('\t')

        sl_chr = fields[2]
        sl_pos = int(fields[3])
        readname = fields[0]

        if sl_chr not in genes:
            continue

        spos = sl_pos // 1000
        if spos not in genes[sl_chr]:
            continue

        for gene_line in genes[sl_chr][spos]:   # for each gene it could be annotated, it will all count once!!
            tmp = gene_line.split('\t')
            pss = int(tmp[0])
            pse = int(tmp[1])
            gene_id = tmp[2]

            if sl_pos < pss or sl_pos > pse:
                continue

            cell_id = readname[-19:-11]     
            umi = readname[-10:]          

            if cell_id not in hash_data:
                hash_data[cell_id] = {}
            if gene_id not in hash_data[cell_id]:
                hash_data[cell_id][gene_id] = {}

            hash_data[cell_id][gene_id][umi] = True   # finally, the number of such True item will be the UMI count, will be a huge number!!!
            genelist[gene_id] = 0
            celllist[cell_id] = 0

    process.stdout.close()
    process.wait()

    print("Passing cell post filter to new hash...")

    cutoff = 3
    total = 0
    new_hash = {}
    new_genelist = {}
    new_celllist = {}

    for cell_id, gene_data in hash_data.items():
        n_cell_umi = 0
        cur_genelist = {}

        for gene_id, umi_data in gene_data.items():
            n_umi = len(umi_data)
            cur_genelist[gene_id] = n_umi
            n_cell_umi += n_umi

        if n_cell_umi >= cutoff:   # the cell have at least 3 UMI
            for gene_id, n_umi in cur_genelist.items():
                total += 1
                if gene_id not in new_hash:
                    new_hash[gene_id] = {}
                new_hash[gene_id][cell_id] = n_umi
                new_genelist[gene_id] = 1
            new_celllist[cell_id] = 1

    print("Cleaning memory...")
    hash_data.clear()
    celllist.clear()
    genelist.clear()

    print("Finished processing.")
    # print(new_hash)
    # print(new_genelist)
    # print(new_celllist)
    print(total)

    # three files!!!

    mtx2_dir = f"{out_prefix}_mtx2"
    if not os.path.exists(mtx2_dir):
        os.makedirs(mtx2_dir)

    # genes.tsv
    genes_file = os.path.join(mtx2_dir, "genes.tsv")
    with open(genes_file, 'w') as out_file:
        gene_order = 0
        for gene_id in new_genelist:
            gene_order += 1
            tmp = gene_id.split(" ")
            output = f"{tmp[0]}\t{tmp[1]}\n"
            out_file.write(output)
            new_genelist[gene_id] = gene_order  # now they are not all 1

    # barcodes.tsv
    barcodes_file = os.path.join(mtx2_dir, "barcodes.tsv")
    with open(barcodes_file, 'w') as out_file:
        cell_order = 0
        for cell_id in new_celllist:
            cell_order += 1
            output = f"{cell_id}\n"
            out_file.write(output)
            new_celllist[cell_id] = cell_order # now they are not all 1

    # matrix.mtx
    matrix_file = os.path.join(mtx2_dir, "matrix.mtx")
    with open(matrix_file, 'w') as out_file:
        out_file.write("%%MatrixMarket matrix coordinate real general\n")
        out_file.write("%\n")
        output = f"{len(new_genelist)} {len(new_celllist)} {total}\n"
        out_file.write(output)

        for gene_id, cell_data in new_hash.items():
            for cell_id, umi_count in cell_data.items():
                gene_index = new_genelist[gene_id]  
                cell_index = new_celllist[cell_id]  
                output = f"{gene_index} {cell_index} {umi_count}\n"
                out_file.write(output)
    print("Files saved successfully.")

    return None

if __name__ == "__main__":
    bam2Mtx2(sys.argv[1],sys.argv[2],sys.argv[3])
