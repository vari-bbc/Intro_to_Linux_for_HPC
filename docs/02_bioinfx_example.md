

# **Mini Project**


Here, we will combine some of the commands we have learned today to perform a toy bioinformatic analysis. Two human paired-end samples (four total fastq files) have already been downloaded from [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) and are in the directory `/varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/fastqs`. To save time for the workshop, each file has been randomly subsetted to just a small fraction of the total reads. Here each user will:

1. Copy these fastq files their own private project directory.
2. Use basic Linux commands to learn some characteristics of these files.
3. Use Salmon to align (pseudo-align) these reads to the hg38 reference transcriptome and get transcripts-per-million (TPMs) for each annotated gene.


## **Start an Interactive Job**

While simple file manipulations can be done on the submit node, computationally intensive operations should be performed using preallocated computational resources. An **interactive job** allows us to preallocate resources while still being able to run commands line by line.

We will start an interactive job, requesting one CPU core and 1 hour and 30 minutes of walltime.


``` bash
srun -p short --nodes=1 --ntasks-per-node=1 --time=01:30:00 --pty bash 
```

After a job has been submitted, users can check on the status (e.g. how long it has been running) of it using the following.


``` bash
squeue --me
```

## **Create Directory**

Typically, one would work in their specific lab's folder on the HPC. For this workshop, we will work in a common directory so that the instructors and check on your progress.


``` bash
cd /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/
```

Create a directory based on your username to separate your work from other users'.


``` bash
mkdir <username>
```

Navigate to the username directory

``` bash
cd <username>
```

Create a `fastqs` directory inside the username directory. It is good practice to keep your raw data separate from your analyses and never make changes to them directly.


``` bash
mkdir fastqs
```

Use `ls` to confirm the `fastqs` directory was created.


``` bash
ls
```

## **Copy Fastq Files To `fastqs` subdirectory**

Verify the current directory is in the username directory.


``` bash
pwd
```

Copy the 4 fastq.gz files and `md5sum.txt`, which you will use to check the validty of the files.


``` bash
cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/fastqs/*fastq.gz ./fastqs/

cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/fastqs/md5sum.txt ./fastqs/
```

Verify that the copied files are valid. The following code should return with an **OK** for each file.


``` bash
cd fastqs

md5sum -c md5sum.txt
```

```
## SRR1039520_1.fastq.gz: OK
## SRR1039520_2.fastq.gz: OK
## SRR1039521_1.fastq.gz: OK
## SRR1039521_2.fastq.gz: OK
```

Go back to the username directory.


``` bash
cd ..
```

## **Use basic Linux commands to explore the fastq files**

Notice all of the sequence data files end with a `.gz` extension, indicating they are gzip compressed. To view such files, one must use `zcat` instead of `cat` so that the files are decompressed into human-readable format.


``` bash
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


``` bash
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

Notice that `SRR1039520_1.fastq.gz` and `SRR1039520_2.fastq.gz` files have the same read IDs. Valid paired fastq files always have matching read IDs throughout the entire file.

Next, we can figure out how many reads are in a fastq file using `wc -l`. Recall that each read is represented by 4 lines in a fastq file, so we need to divide the results of `wc -l` by 4 to get the number of reads.


``` bash
zcat fastqs/SRR1039520_1.fastq.gz | wc -l
```

```
## 2000000
```

Finally, we can use `wc -m` to figure out the read length in a fastq file. Below, we use a combination of `head -n2` and `tail -n1` to get the second line of the first read. You may notice that the results of `wc -m` will be 1 higher than the actual read length. This is because the tool counts the newline character also.


``` bash
zcat fastqs/SRR1039520_1.fastq.gz | head -n2 | tail -n1 | wc -m

```

```
## 64
```

## **Working With Environment Modules**

Type `module av bbc2` (modules beginning with just `bbc` are from the old HPC) to see all the software installed by the BBC. There are a lot of modules installed; to parse through these more easily, we can **pipe** the results of `module av` into `grep` to search for modules with specific keywords. There is a trick, though. The results of `module av` go to standard error, so we need to redirect standard error to standard out using `2>&1` before the pipe.

This command will output any modules containing the `fastqc` keyword.


``` bash
module av -w 1 bbc2 2>&1 | grep 'fastqc'
```

```
## bbc2/fastqc/fastqc-0.12.1
```

## **Run FastQC**

Create a `fastqc` directory inside of your username directory and run FastQC.


``` bash
module load bbc2/fastqc/fastqc-0.12.1

mkdir fastqc

# Run FastQC. You can also run `fastqc -h` to see what different options do.
fastqc -o fastqc fastqs/*fastq.gz
```

```
## Loading bbc2/fastqc/fastqc-0.12.1
##   Loading requirement: VARI/java/1.8.0_202
## application/gzip
## application/gzip
## Started analysis of SRR1039520_1.fastq.gz
## application/gzip
## application/gzip
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

See what was produced by FastQC.


``` bash
ls fastqc/

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

## **Set up a job script to run Salmon to align the reads**

First, we `exit` from our interactive job because we want to get back on to the submit node to submit a non-interactive job to run Salmon.


``` bash
# You should see one job when you run this command, corresponding to your interactive job.
squeue --me

# Exit the interactive job
exit

# Now you should see no jobs when you run this command because the interactive job has ended.
squeue --me

# Go back to your project directory
cd /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/<username>

```

Below is a SLURM job script to run Salmon. For now, do not worry about how the code works. Copy the code and paste it into a new file. Save it as `run_salmon.sh` in your username directory. If you have issues with this task, you can download the job script directly using the command, `wget https://raw.githubusercontent.com/vari-bbc/Intro_to_Linux_for_HPC/refs/heads/main/resources/sec2/run_salmon.sh`.


``` bash
#!/bin/bash

#SBATCH --export=NONE
#SBATCH -J run_salmon
#SBATCH -o run_salmon.o
#SBATCH -e run_salmon.e
#SBATCH --ntasks 4
#SBATCH --time 1:00:00
#SBATCH --mem=31G

start_time=$(date +"%T")

# You need to navigate to your project directory. Conveniently, the $SLURM_SUBMIT_DIR variable stores the path for where the job was submitted.
cd ${SLURM_SUBMIT_DIR}

module load bbc2/salmon/salmon-1.10.0

# Typically, you would have to first build an index before doing the aligning, but we have done this for you already. Here, we store the path to the index file in a variable called 'salmon_idx'.
salmon_idx="/varidata/research/projects/bbc/versioned_references/2022-03-08_14.47.50_v9/data/hg38_gencode/indexes/salmon/hg38_gencode/"

# make output directory for salmon
mkdir -p salmon

# This is called a for loop. We use this to run salmon quant on all the samples, one at a time. It is more efficient to run salmon on each sample "in parallel" but we will not do that today.
for samp_id in SRR1039520 SRR1039521
do
    salmon quant -p ${SLURM_NTASKS} -l A -i $salmon_idx -1 fastqs/${samp_id}_1.fastq.gz -2 fastqs/${samp_id}_2.fastq.gz -o salmon/${samp_id} --validateMappings

done

end_time=$(date +"%T")
echo "Start time: $start_time"
echo "End time: $end_time"
```

## **Submit a Job**

Type `ls` to ensure `run_salmon.sh` file exist in the username directory.


``` bash
ls
```

Use the `sbatch` command to submit the job.


``` bash
sbatch -p short run_salmon.sh

```

Users can check if the job is running with the following command. The job should take about two minutes to complete.


``` bash
squeue --me

```

## **Examine the Standard Output (stdout) and Standard Error (stderr) log files**

It is good practice to check the job logs after your job is done to ensure that the job completed successfully. If there is an error, the output files may not be reliable or could be incomplete.


``` bash
tail run_salmon.e
```

```
## [2025-10-03 10:29:55.282] [jointLog] [info] Marked 0 weighted equivalence classes as degenerate
## [2025-10-03 10:29:55.309] [jointLog] [info] iteration = 0 | max rel diff. = 203.028
## [2025-10-03 10:29:57.974] [jointLog] [info] iteration = 100 | max rel diff. = 19.6442
## [2025-10-03 10:30:00.719] [jointLog] [info] iteration = 200 | max rel diff. = 0.462189
## [2025-10-03 10:30:03.319] [jointLog] [info] iteration = 300 | max rel diff. = 1.92755
## [2025-10-03 10:30:05.909] [jointLog] [info] iteration = 400 | max rel diff. = 0.0347574
## [2025-10-03 10:30:07.538] [jointLog] [info] iteration = 463 | max rel diff. = 0.00959217
## [2025-10-03 10:30:07.581] [jointLog] [info] Finished optimizer
## [2025-10-03 10:30:07.581] [jointLog] [info] writing output
```


``` bash
tail run_salmon.o
```

```
## Start time: 10:28:40
## End time: 10:30:09
```

## **Use grep to find the TPMs for a specific gene**

As an example, let's try to extract out the TPMs (Transcripts per Million) for MUC1. These values can be found in the `quant.sf` file in each sample's folder. 

First, let's take a look at one of these files to figure out the format of these files.


``` bash
head salmon/SRR1039520/quant.sf

```

```
## Name	Length	EffectiveLength	TPM	NumReads
## ENST00000456328.2	1657	1501.734	0.000000	0.000
## ENST00000450305.2	632	476.867	0.000000	0.000
## ENST00000488147.1	1351	1195.734	0.000000	0.000
## ENST00000619216.1	68	2.021	0.000000	0.000
## ENST00000473358.1	712	556.827	0.000000	0.000
## ENST00000469289.1	535	379.986	0.000000	0.000
## ENST00000607096.1	138	24.859	0.000000	0.000
## ENST00000417324.1	1187	1031.734	0.000000	0.000
## ENST00000461467.1	590	434.911	0.000000	0.000
```

From the output above, we can see that the TPMs are in the 4th column of this file. The canonical transcript for [MUC1](http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000185499;r=1:155185824-155192916) is ENST00000620103, so we will search for that using `grep`.


``` bash
grep 'ENST00000620103' salmon/SRR1039520/quant.sf

```

```
## ENST00000620103.4	1811	1655.734	0.000000	0.000
```

Look for MUC1 across all the samples at the same time. We can see that 'SRR1039521' has a TPM of 13.1 for MUC1 compared to 0 for 'SRR1039520'. Recall that the fastq files for this exercise were subsetted to a very small number of reads so don't interpret these results seriously.


``` bash
grep 'ENST00000620103' salmon/*/quant.sf
```

```
## salmon/SRR1039520/quant.sf:ENST00000620103.4	1811	1655.734	0.000000	0.000
## salmon/SRR1039521/quant.sf:ENST00000620103.4	1811	1655.512	13.182825	5.741
```

## **BONUS: Use an interactive job to run multiQC on the Salmon and FastQC output**

In this final step, users will start an interactive job to perform multiQC to collect and summarize the outcomes from FastQC and Salmon, then evaluate the quality of the fastq sequences with the reference transcriptome(mapping rate).

Start an interactive job.


``` bash
srun --nodes=1 --ntasks-per-node=1 --time=00:30:00 --pty bash

```

Navigate to the project directory.

``` bash
cd /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/<username>
```

Load the environment module for multiQC. Run multiQC, which will summarize the results from FastQC and Salmon and output the results into a new directory called `multiqc/`.


``` bash
module load bbc2/multiqc/multiqc-1.14

multiqc --outdir multiqc .
```

```
## 
##   /// MultiQC 🔍 | v1.14
## 
## |           multiqc | MultiQC Version v1.31 now available!
## |           multiqc | Search path : /varidata/research/projects/bbc/research/Intro_to_Linux_for_HPC_devel
## |         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 143/143  
## |            fastqc | Found 4 reports
## |           multiqc | Compressing plot data
## |           multiqc | Report      : multiqc/multiqc_report.html
## |           multiqc | Data        : multiqc/multiqc_data
## |           multiqc | MultiQC complete
```

List the contents of the `multiqc` directory.


``` bash
ls multiqc
```

```
## multiqc_data
## multiqc_report.html
```

Note the newly created `multiqc_report.html` file. Try to view this file in your preferred internet browser. If you have mounted the HPC file system to your computer, you can simply double-click on this file. Alternatively, you can copy this file to your computer's local storage first and then open it.

