

# A toy bioinformatics project

Here, we will put together some of the new skills that we have discussed today to perform a toy bioinformatic analysis. Four RNA-seq fastq files representing paired-end reads from two human samples have already been downloaded from [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778). To save time for the workshop, each file has been randomly subsetted to just a small fraction of the total reads. Here we will:

1. Copy these fastq files to a private project directory for each user.
2. Use basic Linux commands to learn some characteristics of these files, and compare the results to results from FastQC, which we will also run.
3. Use Salmon to align (pseudo-align) these reads to the hg38 reference transcriptome and get transcripts-per-million (TPMs) for each annotated gene.

## Start up an interactive job

Remember we try to avoid running computationally-intensive operations on the submit node directly. Where a formal script is not needed, we can start up an 'interactive job'. Once an interactive job is started, we can use the commandline as before but now we are using pre-allocated computational resources.


```bash
qsub -I -l nodes=1:ppn=1 -l mem=32gb -l walltime=2:00:00

# one way to verify that we have started an interactive successfully is to run qstat
# to show only jobs submitted by ourselves, we use the `-u` option followed by our username.
qstat -u username
```

## Create a project directory for yourself and a subdirectory for storing the raw fastq files

It is good practice to try to make a separate directory for each project to minimize the risk of accidentally editing other files. It is also good to store raw data in its own location and not make any modifications to them.


```bash
cd /varidata/researchtemp/hpctmp/HPC_mini_workshop/Part3

mkdir firstname.lastname

cd firstname.lastname

mkdir fastqs
```

## Copy fastqs to working directory


```bash
# copy the md5sum file which will be used to verify successful file transfers.
cp /varidata/researchtemp/hpctmp/BBC/hpc_workshop_fqs/md5sum.txt ./fastqs/

# copy the fastqs
cp /varidata/researchtemp/hpctmp/BBC/hpc_workshop_fqs/*fastq.gz ./fastqs/

```

## Check that the files transferred properly

Run the following command. The tool will compare the md5sums in the `md5sum.txt` file to the md5sums calculated for the files in your directory currently. If each line says 'OK', then the transfer was successful.


```bash
cd fastqs

md5sum -c md5sum.txt
```

```
## SRR1039520_1.fastq.gz: OK
## SRR1039520_2.fastq.gz: OK
## SRR1039521_1.fastq.gz: OK
## SRR1039521_2.fastq.gz: OK
```

For the next steps, go back to the project directory.


```bash
cd ..
```

## Use zcat to take a look into fastq.gz files

Note that a `.gz` file suffix means that the file is compressed and running the usual `cat` on it will not return human-readable results. The `zcat` command first decompresses the file then prints the results.

Note the matching fastq IDs in the R1 and R2 files. Valid paired fastq files always have matching read IDs throughout the entire file.


```bash
zcat fastqs/SRR1039520_1.fastq.gz | head
```

```
## @SRR1039520.19151550 19151550 length=63
## CCAGGACATCAAGAAGCCAGCTGAAGATGAGTGGGGTAAAACCCCAGACGCCATGAAAGCTGC
## +
## HJJJJJJJJJJJJJJJIJJJJJJJJIJJJJJGGJJJFHIJIIJJIGHHFFFDDDDDCDDDCDD
## @SRR1039520.3404519 3404519 length=63
## TGAGACATGGTTATAGATAAGAGAGTACAAAATGACTCTTTTTCCTGTCAATTGAAATTTAAA
## +
## HIGDFBGIIJBHGIIEHIIJIJIIIEGGIIIHGIJJIJJJIIJGIIIIGIEHIGDCHHIICG7
## @SRR1039520.16253787 16253787 length=63
## CAGGAGACCAAAGACACTGCAATTTGTGTGTTTTCTACAGGGTGCTTTAGATGACGTCTCATT
```

Now look at the R2 file.


```bash
zcat fastqs/SRR1039520_2.fastq.gz | head
```

```
## @SRR1039520.19151550 19151550 length=63
## GAGATGGGGGTCCGTGCGGGCAGAACCCAGGGCATGAAGATCCAAAAGGGCCTGGTTCAGCTT
## +
## HJJJJJJJJJHHHIGIJIJJJJJJJHHHHFFFDDEEDDDDDDDDDDDDDDDDDDDDDDDDDDD
## @SRR1039520.3404519 3404519 length=63
## TTTGACCCTAGTATTGGCAATAGCCCTTTGCTATTTATATAATTAAAACTTTTCTTTAAATTT
## +
## HIJIJJJJIIEHGIGIIIDGIJJJJIIJJJIJJJIJJJIIIJJJJIJJGCHIJJJJJIJIJII
## @SRR1039520.16253787 16253787 length=63
## TACAGTTTGCAAAAGATGTCCAGATGGGTTCTTCTCAAATGAGACGTCATCTAAAGCACCCTG
```

## Use `wc` to see how many reads are in a fastq file and how long they are.

Remember that we can use `man wc` to see what the different options such as `-l` and `-m` do.

This command counts the new line characters at the end of each line. Typically, this corresponds to the number of lines in a file. Recall that each read is represented by 4 lines in a fastq file, so we need to divide the results of `wc -l` by 4 to get the number of reads.


```bash
zcat fastqs/SRR1039520_1.fastq.gz | wc -l
```

```
## 2000000
```

This command counts the number of characters in the sequence line (second line) in the first read in the file. you may notice that the results of `wc -m` will be 1 higher than the actual read length. This is because the tool counts the newline character also.


```bash
zcat fastqs/SRR1039520_1.fastq.gz | head -n2 | tail -n1 | wc -m

```

```
## 64
```

Record these numbers so that you can compare these to the output from FastQC (a commonly used tool for checking the quality of fastq files) which we will run below.

## How to see what packages are installed on HPC.

If you simply type `module av`, you will see all the modules available. There are a lot of modules installed. To parse through these more easily, we can **pipe** the results of `module av` into tools such as `head` to print just the first 'n' lines and `grep` to search for modules with specific keywords. There is a trick, though. The results of `module av` go to stderr, so we need to redirect stderr to stdout using `2>&1` before the pipe.

Print the first 50 lines.


```bash
module av 2>&1 | head -n50
```

```
## 
## ---------------------------- /cm/local/modulefiles -----------------------------
## cluster-tools/7.3 freeipmi/1.5.2    module-git        openldap
## cmd               gcc/6.1.0         module-info       shared
## dot               ipmitool/1.8.17   null              singularity/2.4.2
## 
## ---------------------------- /cm/shared/modulefiles ----------------------------
## acml/gcc/64/5.3.1
## acml/gcc/fma4/5.3.1
## acml/gcc/mp/64/5.3.1
## acml/gcc/mp/fma4/5.3.1
## acml/gcc-int64/64/5.3.1
## acml/gcc-int64/fma4/5.3.1
## acml/gcc-int64/mp/64/5.3.1
## acml/gcc-int64/mp/fma4/5.3.1
## acml/open64/64/5.3.1
## acml/open64/fma4/5.3.1
## acml/open64/mp/64/5.3.1
## acml/open64/mp/fma4/5.3.1
## acml/open64-int64/64/5.3.1
## acml/open64-int64/fma4/5.3.1
## acml/open64-int64/mp/64/5.3.1
## acml/open64-int64/mp/fma4/5.3.1
## aspera/AsperaConnect
## bbc/10x_bamtofastq/bamtofastq-1.2.0
## bbc/10x_bamtofastq/bamtofastq-1.3.2
## bbc/10x_subset-bam/10x_subset-bam-1.0
## bbc/7zip/7zip-16.02
## bbc/abyss/abyss-2.2.3
## bbc/AWS/aws-cli
## bbc/bamtools/bamtools-2.5.1
## bbc/bbc_projects/bbc_projects
## bbc/bcftools/bcftools-1.10.2
## bbc/bcftools/bcftools-1.12
## bbc/bcl2fastq/fastq-multx
## bbc/bcl2fastq2/bcl2fastq2-2.20.0
## bbc/bedops/bedops-2.4.37
## bbc/bedtools/bedtools-2.29.2
## bbc/bedtools/bedtools-2.30.0
## bbc/bioawk/bioawk-git
## bbc/biscuit/biscuit_0_3_14
## bbc/biscuit/biscuit_0_3_16
## bbc/biscuit/biscuit_1_0_0
## bbc/biscuit/biscuit_1_0_1
## bbc/biscuit/biscuit_dev
## bbc/bismark/bismark-0.22.3
## bbc/bismark/bismark-0.23.0
## bbc/blast+/blast+-2.10.0
## bbc/bonito/bonito-0.0.5
## bbc/bowtie/bowtie-1.2.3
```

Print the first 50 lines after subsetting to just the BBC-installed modules.


```bash
module av bbc 2>&1 | head -n50
```

```
## 
## ---------------------------- /cm/shared/modulefiles ----------------------------
## bbc/10x_bamtofastq/bamtofastq-1.2.0
## bbc/10x_bamtofastq/bamtofastq-1.3.2
## bbc/10x_subset-bam/10x_subset-bam-1.0
## bbc/7zip/7zip-16.02
## bbc/abyss/abyss-2.2.3
## bbc/AWS/aws-cli
## bbc/bamtools/bamtools-2.5.1
## bbc/bbc_projects/bbc_projects
## bbc/bcftools/bcftools-1.10.2
## bbc/bcftools/bcftools-1.12
## bbc/bcl2fastq/fastq-multx
## bbc/bcl2fastq2/bcl2fastq2-2.20.0
## bbc/bedops/bedops-2.4.37
## bbc/bedtools/bedtools-2.29.2
## bbc/bedtools/bedtools-2.30.0
## bbc/bioawk/bioawk-git
## bbc/biscuit/biscuit_0_3_14
## bbc/biscuit/biscuit_0_3_16
## bbc/biscuit/biscuit_1_0_0
## bbc/biscuit/biscuit_1_0_1
## bbc/biscuit/biscuit_dev
## bbc/bismark/bismark-0.22.3
## bbc/bismark/bismark-0.23.0
## bbc/blast+/blast+-2.10.0
## bbc/bonito/bonito-0.0.5
## bbc/bowtie/bowtie-1.2.3
## bbc/bowtie2/bowtie2-2.3.5.1
## bbc/bowtie2/bowtie2-2.4.1
## bbc/brename/brename-v2.11.1
## bbc/bustools/bustools-0.39.3
## bbc/bustools/bustools-0.40.0
## bbc/bwa/bwa-0.7.17
## bbc/bwa-mem2/bwa-mem2-v2.2.1
## bbc/byacc/byacc-1.9
## bbc/canu/canu-1.9
## bbc/CaVEMan/CaVEMan-1.15.1
## bbc/cellranger/cellranger-3.0.2
## bbc/cellranger/cellranger-3.1.0
## bbc/cellranger/cellranger-4.0.0
## bbc/cellranger/cellranger-6.0.2
## bbc/cellranger/cellranger-6.1.2
## bbc/cellranger-atac/cellranger-atac-1.1.0
## bbc/changeo/changeo-1.1.0
## bbc/choose/choose-1.3.3
## bbc/chromap/chromap
## bbc/ChromHMM/chromHMM-1.22
## bbc/ChromHMM/chromHMM-1.23
## bbc/cmake/cmake-3.19.4
```

Print the modules with the keyword, 'fastqc', in their names.


```bash
module av 2>&1 | grep 'fastqc'
```

```
## bbc/fastqc/fastqc-0.11.8
## bbc/fastqc/fastqc-0.11.9
```

## Run FastQC on the fastq files

We run fastqc on our raw fastq files for essentially all of the projects that we work on.

The output of fastqc includes an html file for each fastq file. You can open up these in a regular browser to look at basic QC metrics for these fastq files. In the interest of time, we will simply run FastQC without viewing the results because we will aggregate all of our results into one html file at the end of this exercise.


```bash
module load bbc/fastqc/fastqc-0.11.9

# make a directory for the output
mkdir fastqc

# run FastQC. Remember you can check out what the options do by typing `fastqc -h`.
fastqc -o fastqc fastqs/*fastq.gz
```

```
## mkdir: cannot create directory ‘fastqc’: File exists
## Started analysis of SRR1039520_1.fastq.gz
## Approx 5% complete for SRR1039520_1.fastq.gz
## Approx 10% complete for SRR1039520_1.fastq.gz
## Approx 15% complete for SRR1039520_1.fastq.gz
## Approx 20% complete for SRR1039520_1.fastq.gz
## Approx 25% complete for SRR1039520_1.fastq.gz
## Approx 30% complete for SRR1039520_1.fastq.gz
## Approx 35% complete for SRR1039520_1.fastq.gz
## Approx 40% complete for SRR1039520_1.fastq.gz
## Approx 45% complete for SRR1039520_1.fastq.gz
## Approx 50% complete for SRR1039520_1.fastq.gz
## Approx 55% complete for SRR1039520_1.fastq.gz
## Approx 60% complete for SRR1039520_1.fastq.gz
## Approx 65% complete for SRR1039520_1.fastq.gz
## Approx 70% complete for SRR1039520_1.fastq.gz
## Approx 75% complete for SRR1039520_1.fastq.gz
## Approx 80% complete for SRR1039520_1.fastq.gz
## Approx 85% complete for SRR1039520_1.fastq.gz
## Approx 90% complete for SRR1039520_1.fastq.gz
## Approx 95% complete for SRR1039520_1.fastq.gz
## Approx 100% complete for SRR1039520_1.fastq.gz
## Analysis complete for SRR1039520_1.fastq.gz
## Started analysis of SRR1039520_2.fastq.gz
## Approx 5% complete for SRR1039520_2.fastq.gz
## Approx 10% complete for SRR1039520_2.fastq.gz
## Approx 15% complete for SRR1039520_2.fastq.gz
## Approx 20% complete for SRR1039520_2.fastq.gz
## Approx 25% complete for SRR1039520_2.fastq.gz
## Approx 30% complete for SRR1039520_2.fastq.gz
## Approx 35% complete for SRR1039520_2.fastq.gz
## Approx 40% complete for SRR1039520_2.fastq.gz
## Approx 45% complete for SRR1039520_2.fastq.gz
## Approx 50% complete for SRR1039520_2.fastq.gz
## Approx 55% complete for SRR1039520_2.fastq.gz
## Approx 60% complete for SRR1039520_2.fastq.gz
## Approx 65% complete for SRR1039520_2.fastq.gz
## Approx 70% complete for SRR1039520_2.fastq.gz
## Approx 75% complete for SRR1039520_2.fastq.gz
## Approx 80% complete for SRR1039520_2.fastq.gz
## Approx 85% complete for SRR1039520_2.fastq.gz
## Approx 90% complete for SRR1039520_2.fastq.gz
## Approx 95% complete for SRR1039520_2.fastq.gz
## Approx 100% complete for SRR1039520_2.fastq.gz
## Analysis complete for SRR1039520_2.fastq.gz
## Started analysis of SRR1039521_1.fastq.gz
## Approx 5% complete for SRR1039521_1.fastq.gz
## Approx 10% complete for SRR1039521_1.fastq.gz
## Approx 15% complete for SRR1039521_1.fastq.gz
## Approx 20% complete for SRR1039521_1.fastq.gz
## Approx 25% complete for SRR1039521_1.fastq.gz
## Approx 30% complete for SRR1039521_1.fastq.gz
## Approx 35% complete for SRR1039521_1.fastq.gz
## Approx 40% complete for SRR1039521_1.fastq.gz
## Approx 45% complete for SRR1039521_1.fastq.gz
## Approx 50% complete for SRR1039521_1.fastq.gz
## Approx 55% complete for SRR1039521_1.fastq.gz
## Approx 60% complete for SRR1039521_1.fastq.gz
## Approx 65% complete for SRR1039521_1.fastq.gz
## Approx 70% complete for SRR1039521_1.fastq.gz
## Approx 75% complete for SRR1039521_1.fastq.gz
## Approx 80% complete for SRR1039521_1.fastq.gz
## Approx 85% complete for SRR1039521_1.fastq.gz
## Approx 90% complete for SRR1039521_1.fastq.gz
## Approx 95% complete for SRR1039521_1.fastq.gz
## Approx 100% complete for SRR1039521_1.fastq.gz
## Analysis complete for SRR1039521_1.fastq.gz
## Started analysis of SRR1039521_2.fastq.gz
## Approx 5% complete for SRR1039521_2.fastq.gz
## Approx 10% complete for SRR1039521_2.fastq.gz
## Approx 15% complete for SRR1039521_2.fastq.gz
## Approx 20% complete for SRR1039521_2.fastq.gz
## Approx 25% complete for SRR1039521_2.fastq.gz
## Approx 30% complete for SRR1039521_2.fastq.gz
## Approx 35% complete for SRR1039521_2.fastq.gz
## Approx 40% complete for SRR1039521_2.fastq.gz
## Approx 45% complete for SRR1039521_2.fastq.gz
## Approx 50% complete for SRR1039521_2.fastq.gz
## Approx 55% complete for SRR1039521_2.fastq.gz
## Approx 60% complete for SRR1039521_2.fastq.gz
## Approx 65% complete for SRR1039521_2.fastq.gz
## Approx 70% complete for SRR1039521_2.fastq.gz
## Approx 75% complete for SRR1039521_2.fastq.gz
## Approx 80% complete for SRR1039521_2.fastq.gz
## Approx 85% complete for SRR1039521_2.fastq.gz
## Approx 90% complete for SRR1039521_2.fastq.gz
## Approx 95% complete for SRR1039521_2.fastq.gz
## Approx 100% complete for SRR1039521_2.fastq.gz
## Analysis complete for SRR1039521_2.fastq.gz
```

Look at the output from FastQC.


```bash
ls fastqc

```

```
## SRR1039520_1_fastqc.html
## SRR1039520_1_fastqc.zip
## SRR1039520_2_fastqc.html
## SRR1039520_2_fastqc.zip
## SRR1039521_1_fastqc.html
## SRR1039521_1_fastqc.zip
## SRR1039521_2_fastqc.html
## SRR1039521_2_fastqc.zip
```

## Set up a job to run Salmon

First, we `exit` from our interactive job because we want to get back on to the submit node to submit a non-interactive job to run Salmon.


```bash
# exit the interactive job
exit

# verify that the interactive job has ended
qstat -u username

# go back to the project directory
cd /varidata/researchtemp/hpctmp/HPC_mini_workshop/Part3/firstname.lastname

```

Next, we can set up our job script to run Salmon. The job script code is below. Don't worry about the specifics of this code for now. Simply, copy the code and paste into your text editor to make your job script. Save the script as `run_salmon.sh` and ensure that it is in our project directory containing the `fastqs/` subfolder.


```bash
#PBS -l nodes=1:ppn=3
#PBS -l mem=36gb
#PBS -l walltime=1:00:00
#PBS -N salmon
#PBS -o salmon.o
#PBS -e salmon.e

set -euo pipefail # see explanation at https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425?permalink_comment_id=3945021

start_time=$(date +"%T")

# You need to navigate to your project directory. Conveniently, the $PBS_O_WORKDIR variable stores the path for where the job was submitted.
cd ${PBS_O_WORKDIR}

module load bbc/salmon/salmon-1.5.2

# Typically, you would have to first build an index before doing the aligning, but we have done this for you already. Here, we store the path to the index file in a variable called 'salmon_idx'.
salmon_idx="/varidata/research/projects/bbc/versioned_references/2022-03-08_14.47.50_v9/data/hg38_gencode/indexes/salmon/hg38_gencode/"

# make output directory for salmon
mkdir -p salmon

# this is called a for loop. We use this to run salmon quant on all the samples, one at a time. It is more efficient to run salmon on each sample "in parallel" but we will not do not today.
for samp_id in SRR1039520 SRR1039521
do
    salmon quant -p ${PBS_NUM_PPN} -l A -i $salmon_idx -1 fastqs/${samp_id}_1.fastq.gz -2 fastqs/${samp_id}_2.fastq.gz -o salmon/${samp_id} --validateMappings

done

end_time=$(date +"%T")
echo "Start time: $start_time"
echo "End time: $end_time"
```

## Submit the job

Submit the job and keep monitoring the job using `qstat` as shown below. The job should complete in about 2 minutes.


```bash
qsub run_salmon.sh

# check on the status of your job
qstat -u firstname.lastname

```

## Check the job logs to see if job finished running

It is good practice to check the job logs after your job is done to ensure that the job completed successfully. If there is an error, the output files may not be reliable or could be incomplete.

Look at the stderr from the script.


```bash
tail salmon.e
```

```
## [2022-09-20 23:32:16.397] [jointLog] [info] Marked 0 weighted equivalence classes as degenerate
## [2022-09-20 23:32:16.426] [jointLog] [info] iteration = 0 | max rel diff. = 202.929
## [2022-09-20 23:32:18.786] [jointLog] [info] iteration = 100 | max rel diff. = 6.3715
## [2022-09-20 23:32:21.080] [jointLog] [info] iteration = 200 | max rel diff. = 5.04421
## [2022-09-20 23:32:23.395] [jointLog] [info] iteration = 300 | max rel diff. = 0.631388
## [2022-09-20 23:32:25.725] [jointLog] [info] iteration = 400 | max rel diff. = 2.25551
## [2022-09-20 23:32:27.268] [jointLog] [info] iteration = 467 | max rel diff. = 0.00956146
## [2022-09-20 23:32:27.307] [jointLog] [info] Finished optimizer
## [2022-09-20 23:32:27.307] [jointLog] [info] writing output
```

Look at the stdout from the script.


```bash
tail salmon.o

```

```
## Start time: 23:30:31
## End time: 23:32:30
```

## Use grep to find the TPMs for specific genes

As an example, let's try to extract out the TPMs for MUC1. These values can be found in the `quant.sf` file in each sample's folder. Let's take a look at one of these files to figure out the format of these files.


```bash
head salmon/SRR1039520/quant.sf 

```

```
## Name	Length	EffectiveLength	TPM	NumReads
## ENST00000456328.2	1657	1502.015	0.000000	0.000
## ENST00000450305.2	632	477.136	0.000000	0.000
## ENST00000488147.1	1351	1196.015	0.000000	0.000
## ENST00000619216.1	68	1.932	0.000000	0.000
## ENST00000473358.1	712	557.094	0.000000	0.000
## ENST00000469289.1	535	380.252	0.000000	0.000
## ENST00000607096.1	138	24.897	0.000000	0.000
## ENST00000417324.1	1187	1032.015	0.000000	0.000
## ENST00000461467.1	590	435.180	0.000000	0.000
```

We can see that the TPMs are in the 4th column of this file. We can also see that the transcript IDs are in ENSEMBL format. If we look up MUC1 on [Ensembl](http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000185499;r=1:155185824-155192916), we will find that the canonical transcript is `ENST00000620103`, so let's `grep` for this transcript ID. 

Look for this transript for one specific sample.


```bash
grep 'ENST00000620103' salmon/SRR1039520/quant.sf
```

```
## ENST00000620103.4	1811	1656.015	0.000000	0.000
```

Look across all the samples at the same time. We can see that 'SRR1039521' has a TPM of 13.1 for MUC1 compared to 0 for 'SRR1039520'. Recall that the fastq files for this exercise were subsetted to a very small number of reads so don't interpret these results seriously.


```bash
grep 'ENST00000620103' salmon/*/quant.sf
```

```
## salmon/SRR1039520/quant.sf:ENST00000620103.4	1811	1656.015	0.000000	0.000
## salmon/SRR1039521/quant.sf:ENST00000620103.4	1811	1655.481	13.145898	5.723
```

## Use an interactive job to run multiQC on the Salmon and FastQC output

As our final exercise, let's run a tool called multiQC to summarize the results from FastQC and Salmon. Using the output from multiQC, we will be able to assess both the quality of the fastq sequences themselves and how well they match up with the reference transcriptome (mapping rate), and other useful information. We will start up another interactive job to do this.


```bash
qsub -I -l nodes=1:ppn=1 -l mem=32gb -l walltime=1:00:00

# go back to your project directory
cd /varidata/researchtemp/hpctmp/HPC_mini_workshop/Part3/firstname.lastname/
```

Now actually run multiQC.


```bash
module load bbc/multiqc/multiqc-1.12

# multiqc creates the output directory automatically. You don't have to run mkdir manually.
multiqc --outdir multiqc .
```

```
## 
## ### Loaded BBC module
## 	Loading this module prepends to $PYTHONPATH
## 	Don't use with conda.
## ### End BBC module message.
## 
##   /// MultiQC 🔍 | v1.12
## 
## |           multiqc | Search path : /varidata/research/home/daisy.fu/Intro_to_Linux_for_HPC
## |         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 141/141  
## |            snippy | Found 1 reports
## |          bargraph | Tried to make bar plot, but had no data: snippy_variants
## |            salmon | Found 2 meta reports
## |            salmon | Found 2 fragment length distributions
## |            fastqc | Found 4 reports
## |           multiqc | Compressing plot data
## |           multiqc | Previous MultiQC output found! Adjusting filenames..
## |           multiqc | Use -f or --force to overwrite existing reports instead
## |           multiqc | Report      : multiqc/multiqc_report_2.html
## |           multiqc | Data        : multiqc/multiqc_data_2
## |           multiqc | MultiQC complete
```

Look at the output from multiQC.


```bash
ls multiqc
```

```
## multiqc_data
## multiqc_data_1
## multiqc_data_2
## multiqc_report_1.html
## multiqc_report_2.html
## multiqc_report.html
```

Note the newly created `multiqc_report.html` file. Your final task today is to view this file in your browser. If you have mounted the HPC file system to your computer, you can try to open up this file directly. Alternatively, you can copy this file to your computer's local storage first and then open it.
