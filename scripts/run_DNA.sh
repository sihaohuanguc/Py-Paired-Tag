#!/bin/bash

in_folder="in_fastq_path"  # the folder that contains fastq files
s="sample_prefix"    # sample name prefix in the read 1 and read 2 fastq files
bowtie_index="path_to_bowtie_index"  # bowtie index prefix
bowtie2_index="path_to_bowtie2_index"    # bowtie2 index prefix
special_annotation="path_to_annotation_file"  # special annotation file for the bam2Mtx2 step, examples could be found in the original Paired-Tag package

mkdir -p DNA_add_barcode_results
mkdir -p DNA_trimmed_results
mkdir -p DNA_alignment_results
mkdir -p DNA_sparse_matrix_results

ppt combine3 -1 ${in_folder}/${s}_R1.fq.gz -2 ${in_folder}/${s}_R2.fq.gz -p DNA_add_barcode_results/${s}
zcat DNA_add_barcode_results/${s}_combined.fq.gz | bowtie ${bowtie_index} - --norc -m 1 -v 1 -S DNA_add_barcode_results/${s}_BC.sam
ppt convert2 -i DNA_add_barcode_results/${s}_BC.sam -p DNA_add_barcode_results/${s}
trim_galore DNA_add_barcode_results/${s}_cov.fq.gz -o DNA_trimmed_results
bowtie2 -x ${bowtie2_index} -U DNA_trimmed_results/${s}_cov_trimmed.fq.gz --no-unal -p 8 -S DNA_alignment_results/${s}.sam
samtools sort DNA_alignment_results/${s}.sam -o DNA_alignment_results/${s}_sorted.bam
ppt rmdup2 -i DNA_alignment_results/${s}_sorted.bam -p DNA_alignment_results/${s}
ppt bam2Mtx2 -r ${special_annotation} -i DNA_alignment_results/${s}_rmdup.bam -p DNA_sparse_matrix_results/${s}

