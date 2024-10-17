#!/bin/bash

set -euo pipefail

# 01
mkdir -p ~/hpc_workshop_2024
/usr/bin/cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/* ~/hpc_workshop_2024
[[ -f ~/hpc_workshop_2024/metadata.tsv ]] || { cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/metadata.tsv ~/hpc_workshop_2024; }
[[ -f ~/hpc_workshop_2024/data_01_R1.fq ]] || { echo "~/hpc_workshop_2024/data_01_R1.fq does not exist."; exit 1; }
[[ -f ~/hpc_workshop_2024/data_54_R1.fq ]] || { echo "~/hpc_workshop_2024/data_54_R1.fq does not exist."; exit 1; }
rm ~/combined.fq ~/lines_with_13.tsv

# 02
section02_files="fastqs/SRR1039520_1.fastq.gz  fastqs/SRR1039520_2.fastq.gz  fastqs/SRR1039521_1.fastq.gz  fastqs/SRR1039521_2.fastq.gz  fastqs/md5sum.txt run_salmon.e run_salmon.o salmon/SRR1039520/quant.sf salmon/SRR1039521/quant.sf"
for file in $section02_files
do
    [[ -f $file ]] || { echo "${file} does not exist."; exit 1; }
done

rm -r fastqc

# 3
section03_files="SummarizedExperiment.rds"
for file in $section02_files
do
    [[ -f $file ]] || { echo "${file} does not exist."; exit 1; }
done


echo "Good to go."
# module load bbc2/R/alt/R-4.2.1-setR_LIBS_USER
# Rscript -e 'bookdown::render_book("index.Rmd")'
