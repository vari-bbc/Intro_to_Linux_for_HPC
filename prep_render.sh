#!/bin/bash

set -euo pipefail

# 01

[[ -f ~/hpc_workshop/test_01_R1.fq ]] || { echo "~/hpc_workshop/test_01_R1.fq does not exist."; exit 1; }
[[ -f test_01_R1.fq ]] || { echo "test_01_R1.fq does not exist."; exit 1; }
[[ -f test_54_R1.fq ]] || { echo "test_54_R1.fq does not exist."; exit 1; }
rm -f combined.fq

# 02
section02_files="fastqs/SRR1039520_1.fastq.gz  fastqs/SRR1039520_2.fastq.gz  fastqs/SRR1039521_1.fastq.gz  fastqs/SRR1039521_2.fastq.gz  fastqs/md5sum.txt run_salmon.e run_salmon.o salmon/SRR1039520/quant.sf salmon/SRR1039521/quant.sf"
for file in $section02_files
do
    [[ -f $file ]] || { echo "${file} does not exist."; exit 1; }
done

rm -fr fastqc









echo "Good to go."
# module load bbc2/R/alt/R-4.2.1-setR_LIBS_USER
# Rscript -e 'bookdown::render_book("index.Rmd")'
