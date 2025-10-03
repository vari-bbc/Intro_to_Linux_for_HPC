#!/bin/bash

#SBATCH --export=NONE
#SBATCH -J run_salmon
#SBATCH -o run_salmon.o
#SBATCH -e run_salmon.e
#SBATCH --ntasks 4
#SBATCH --time 1:00:00
#SBATCH --mem=36G

start_time=$(date +"%T")

# You need to navigate to your project directory. Conveniently, the $SLURM_SUBMIT_DIR variable stores the path for where the job was submitted.
cd ${SLURM_SUBMIT_DIR}

module load bbc2/salmon/salmon-1.10.0

# Typically, you would have to first build an index before doing the aligning, but we have done this for you already. Here, we store the path to the index file in a variable called 'salmon_idx'.
salmon_idx="/varidata/research/projects/bbc/versioned_references/2022-10-06_14.25.40_v10/data/hg38_gencode/indexes/salmon/hg38_gencode/"

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
