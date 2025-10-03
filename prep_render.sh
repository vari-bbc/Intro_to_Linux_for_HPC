#!/bin/bash

set -euo pipefail

# 01
mkdir -p ~/hpc_workshop_2024
[[ -f ~/hpc_workshop_2024/metadata.tsv ]] || cp "resources/sec1/metadata.tsv" "~/hpc_workshop_2024/"
[[ -f ~/hpc_workshop_2024/data_01_R1.fq ]] || cp "resources/sec1/data_01_R1.fq" "~/hpc_workshop_2024/" 
[[ -f ~/hpc_workshop_2024/data_54_R1.fq ]] || cp "resources/sec1/data_54_R1.fq" "~/hpc_workshop_2024/" 
rm -f ~/combined.fq ~/lines_with_13.tsv

# 02
section02_files="fastqs/SRR1039520_1.fastq.gz  fastqs/SRR1039520_2.fastq.gz  fastqs/SRR1039521_1.fastq.gz  fastqs/SRR1039521_2.fastq.gz  fastqs/md5sum.txt run_salmon.e run_salmon.o salmon/SRR1039520/quant.sf salmon/SRR1039521/quant.sf"
for file in $section02_files
do
    dir=`dirname "${file}"`
    [[ -f $file ]] || { mkdir -p $dir && cp "resources/sec2/${file}"  "${file}"; }
done

rm -fr fastqc
rm -fr multiqc

# 3
section03_files="SummarizedExperiment.rds"
for file in $section03_files
do
    [[ -f $file ]] || cp "resources/sec3/${file}" "./"
done

# 4
section04_files="de_res.tsv"
for file in $section04_files
do
    [[ -f $file ]] || cp "resources/sec4/${file}" "./"
done

echo "Good to go."
# module load bbc2/R/alt/R-4.4.0-setR_LIBS_USER
# Rscript -e 'bookdown::render_book("index.Rmd")'
