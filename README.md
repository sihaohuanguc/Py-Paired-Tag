# Description
Py-Paired-Tag (PPT) is the python version of [Paired-Tag](https://github.com/cxzhu/Paired-Tag) pipeline. It is used for Paired-Tag sequencing data analysis. Paired-Tag is a single cell technology that could simultaneously map one type of histone modification (or chromatin accessibility) together with nuclear RNA expressions. The output matrix file from this package contains 100% same information as the original C++ version, as long as they are run in the same environment, although the order of genes and barcodes could be different. This is nearly a mirror of the original Paired-Tag pipeline, with slight change in input/output strategy and fixation of several small bugs not related to the final matrix output. This is mainly for those who are more familiar with python than C++ and want to customize the data processing pipeline for their own purposes.

# Requirements
The package list and versions could be found in `for_py_paired_tag.yml`. To install all the packages necessary for the pipeline, run the following command. If you want to customize the name of the environment, just open the YAML file and change the `name` in the first line. Then run the commands below.
```bash
conda env create -f for_py_paired_tag.yml
```
And then activate the environment by running
```bash
conda activate for_py_paired_tag  # or replace "for_py_paired_tag" if you customized the name of the environment
```
In this way, your environment will be ready to install and run PPT. 

# Download
You could download the package by the following command.
```bash
git clone https://github.com/sihaohuanguc/Py-Paired-Tag.git
```
Then go to the folder with the `setup.py` file. And run
```bash
pip3 install .
```
Now you've installed the package. You could use it at any place of your account. If you are not clear about auy command, you could find help by the command
```bash
ppt -h
```
**Please follow the guidance and obey the rules of your cluster, when running the commands.**

# Reference and citation
This package is derived from the original C++ version of [Paired-Tag](https://github.com/cxzhu/Paired-Tag). You are always welcome to use the package or edit it in your work. If you would like to apply this package in your work, please refer to [this repository](https://github.com/cxzhu/Paired-Tag) for citation information.

# Protocol
## Pre-written scripts
If you would like to run the package directly without any customization, please follow the instructions in this section.
### DNA
In this case, your sample is Paired-Tag V1 histone or ATAC-seq data. Your cell barcode and UMI are on read 2. Only read 1 is used for genome alignment (single end mode). Run the `run_DNA.sh` in the `scripts` folder. Change the variables in the script. If the Bowtie/Bowtie2 index are not ready, see below in the `Preparation steps` section. Examples for the special annotation file for the bam2Mtx2 step could be found in the `reference` folder of the original [Paired-Tag](https://github.com/cxzhu/Paired-Tag) repository.
### RNA
In this case, your sample is Paired-Tag V1 RNA-seq data. Your cell barcode and UMI are on read 2. Only read 1 is used for genome alignment (single end mode). Run the `run_RNA.sh` in the `scripts` folder. Change the variables in the script. If the Bowtie/STAR index are not ready, see below in the `Preparation steps` section. Examples for the special annotation file for the bam2Mtx2 step could be found in the `reference` folder of the original [Paired-Tag](https://github.com/cxzhu/Paired-Tag) repository.

## Preparation steps
### Make index for bowtie for cell barcode whitelist
Examples for the cell barcode ID fasta files could be found in the `reference` folder of the original [Paired-Tag](https://github.com/cxzhu/Paired-Tag) repository. An example command line to build the index is shown below.
```bash
bowtie-build ./cell_id_full.fa cell_id/cell_id
```
### Make index for bowtie2 for alignment
To make index files for DNA alignment to genome by bowtie2, check the following commands as an example. In this case, the species is mouse.
```bash
mkdir bowtie2_genome_ref_mouse
bowtie2-build GRCm39.genome.fa bowtie2_genome_ref_mouse/mouse
```
### Make index for STAR for alignment
To make index files for RNA alignment to genome by STAR, check the following commands as an example. In this case, the species is mouse. The read length is 100nt.
```bash
mkdir index_overhang_99_gtf_mouse
STAR --runMode genomeGenerate \
     --genomeDir index_overhang_99_gtf_mouse \
     --runThreadN 8 \
     --genomeFastaFiles GRCm39.genome.fa \
     --sjdbGTFfile gencode.vM36.annotation.gtf \
     --sjdbOverhang 99
```

## Step by step protocol
### 1. Find out cell barcode and UMI for each read
The first step is to extract barcode and UMI information from read 2 and store them as fastq format. At the same time, read 1 sequence and quality scores are added to the first line of a fastq record, separated by ":".
```bash
ppt combine3 -1 your_read1.fq.gz -2 your_read2.fq.gz -p out_prefix
```
`-1` and `-2` are compressed read 1 and 2 fastq files, they are supposed to be in `.gz` format. This step will generate a report file `{out_prefix}_type_report.xls.gz` and a fastq file `{out_prefix}_combined.fq.gz`. Then we use bowtie to find out the ID for the cell barcode of each read. Here the bowtie index is for cell barcodes and could be generated by the command in the `Preparation steps` section. The output is in `.sam` format.
```bash
zcat your_prefix_combined.fq.gz | bowtie your_bowtie_index - --norc -m 1 -v 1 -S your_prefix_BC.sam
```
Next we generate the read 1 fastq file with cell barcode and UMI in the first line of each record, separated by ":".
```bash
ppt convert2 -i your_prefix_BC.sam -p out_prefix
```
### 2. Trim reads
For DNA, we just do the standard trim for read 1.
```bash
trim_galore your_prefix_cov.fq.gz -o your_output_folder
```
For RNA, we do the standard trim and also Oligo-dT primer and N6 primer trim for read 1.
```bash
trim_galore your_prefix_cov.fq.gz -o your_output_folder
trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN your_prefix_cov_trimmed.fq.gz -o your_output_folder
trim_galore -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN your_prefix_cov_trimmed_trimmed.fq.gz -o your_output_folder
```
### 3. Align reads to genome and generate the sparse matrix
The first step in this part is to do the alignment. For DNA, we use bowtie2 to align read 1 to the genome. You could generate the bowtie2 index files following `Preparation steps` section. After alignment, we use samtools to sort the output `.sam` file.
```bash
bowtie2 -x your_bowtie2_index -U your_prefix_cov_trimmed.fq.gz --no-unal -p 8 -S your_prefix.sam
samtools sort your_prefix.sam -o your_prefix_sorted.bam
```
For RNA, we use STAR to align read 1 to the genome. You could generate the STAR index files following `Preparation steps` section. After alignment, we use samtools to remove the secondary alignment lines and then sort the output `.sam` file.
```bash
STAR --runThreadN 8 --genomeDir your_star_index --readFilesIn your_prefix_cov_trimmed_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix out_prefix_ --outSAMtype BAM Unsorted
samtools view -h -F 256 your_prefix_Aligned.out.bam -b > your_prefix_clean.bam
samtools sort your_prefix_clean.bam -o your_prefix_sorted.bam
```
Next, for both DNA and RNA samples, we remove PCR introduced dulicate reads. Only one read will be kept when they share the same mapped chromosome, mapping start position and UMI. This step will generate a report file `{out_prefix}_dup.xls` and a bam file `{out_prefix}_rmdup.bam`.
```bash
ppt rmdup2 -i your_prefix_sorted.bam -p out_prefix
```
Then we convert the bam file into the expression matrix. The special annotation file derived from `.gtf` annotation files for DNA and RNA could be found in the `reference` folder of the original [Paired-Tag](https://github.com/cxzhu/Paired-Tag) repository. If you would like to make your own annotation files, you could use `tools/make_special_annotation_human.py` in this repository as a template. This script is for `.gtf` files downloaded from [GENCODE](https://www.gencodegenes.org) database. For annotation files from other databases, please customize the script.
```bash
ppt bam2Mtx2 -r your_special_annotation -i your_prefix_rmdup.bam -p out_prefix
```
This step will generate a folder ending with `mtx2` containing `genes.tsv`, `barcodes.tsv` and `matrix.mtx`. These files could be used in the downstream analysis by software like seurat and scanpy.



