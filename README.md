# Description
Py-Paired-Tag (PPT) is the python version of [Paired-Tag](https://github.com/cxzhu/Paired-Tag) pipeline. It is used for Paired-Tag sequencing data analysis. Paired-Tag is a single cell technology that could simultaneously map one type of histone modification (or chromatin accessibility) together with nuclear RNA expressions. The output matrix file from this package contains 100% same information as the original C++ version, as long as they are run in the same environment, although the order of genes and barcodes could be different. This is nearly a mirror of the original Paired-Tag pipeline, with slight change in input/output strategy and fixation of several small bugs not related to the final matrix output. This is mainly for those who are more familiar with python than C++ and want to customize the data processing pipeline for their own purposes.

# Requirements
The package list and versions could be found in `for_py_paired_tag.yml`. To install all the packages necessary for the pipeline, run the following command
```bash
conda env create -f for_py_paired_tag.yml
```
And then activate the environment by running
```bash
conda activate for_py_paired_tag  # or replace "for_py_paired_tag" if you customized the name of the environment
```
In this way, your environment will be ready to install and run PPT. If you want to customize the name of the environment, just open the YAML file and change the `name` in the first line. Then run the commands above.

# Download
You could download the package by the following command.
```bash
git clone (Placeholder)()()()()()
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
This package is derived from the cpp version of [Paired-Tag](https://github.com/cxzhu/Paired-Tag). You are always welcome to use the package or edit it in your work. If you would like to apply this package in your work, please refer to [this repository](https://github.com/cxzhu/Paired-Tag) for citation.

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
STAR --runMode genomeGenerate --genomeDir index_overhang_99_gtf_mouse --runThreadN 8 --genomeFastaFiles GRCm39.genome.fa --sjdbGTFfile gencode.vM36.annotation.gtf --sjdbOverhang 99
```

## Step by step protocol


