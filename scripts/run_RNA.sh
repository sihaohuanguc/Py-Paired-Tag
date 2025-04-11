#!/bin/bash

in_folder="in_fastq_path"  # the folder that contains fastq files
s="sample_prefix"    # sample name prefix in the read 1 and read 2 fastq files
bowtie_index="path_to_bowtie_index"  # bowtie index prefix
star_index="path_to_star_index"    # STAR index prefix
special_annotation="path_to_annotation_file"  # special annotation file for the bam2Mtx2 step, examples could be found in the original Paired-Tag package

mkdir -p RNA_add_barcode_results
mkdir -p RNA_trimmed_results
mkdir -p RNA_alignment_results
mkdir -p RNA_sparse_matrix_results

ppt combine3 -1 ${in_folder}/${s}_R1.fq.gz -2 ${in_folder}/${s}_R2.fq.gz -p RNA_add_barcode_results/${s}
zcat RNA_add_barcode_results/${s}_combined.fq.gz | bowtie ${bowtie_index} - --norc -m 1 -v 1 -S RNA_add_barcode_results/${s}_BC.sam
ppt convert2 -i RNA_add_barcode_results/${s}_BC.sam -p RNA_add_barcode_results/${s}
trim_galore RNA_add_barcode_results/${s}_cov.fq.gz -o RNA_trimmed_results
trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN RNA_trimmed_results/${s}_cov_trimmed.fq.gz -o RNA_trimmed_results
trim_galore -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN RNA_trimmed_results/${s}_cov_trimmed_trimmed.fq.gz -o RNA_trimmed_results
STAR --runThreadN 8 --genomeDir ${star_index} --readFilesIn RNA_trimmed_results/${s}_cov_trimmed_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix RNA_alignment_results/${s}_ --outSAMtype BAM Unsorted
samtools view -h -F 256 RNA_alignment_results/${s}_Aligned.out.bam -b > RNA_alignment_results/${s}_clean.bam
samtools sort RNA_alignment_results/${s}_clean.bam -o RNA_alignment_results/${s}_sorted.bam
ppt rmdup2 -i RNA_alignment_results/${s}_sorted.bam -p RNA_alignment_results/${s}
ppt bam2Mtx2 -r ${special_annotation} -i RNA_alignment_results/${s}_rmdup.bam -p RNA_sparse_matrix_results/${s}
