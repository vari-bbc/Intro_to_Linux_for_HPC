

# A toy bioinformatics project

## Navigate to your research group's folder and see what's there


```bash
cd /varidata/research/projects/labname

ls
```

## Create project directory and navgiate into it

Also make a directory to store your fastq files.


```bash
mkdir hpc_workshop_202209 

cd hpc_workshop_202209

mkdir fastqs
```

## Copy fastqs to working directory


```bash
cp /varidata/researchtemp/hpctmp/BBC/hpc_workshop_fqs/md5sum.txt ./fastqs/
cp /varidata/researchtemp/hpctmp/BBC/hpc_workshop_fqs/*fastq.gz ./fastqs/

```

## Check that the files transferred properly


```bash
ls fastqs

cd fastqs

md5sum -c md5sum.txt
```

```
## md5sum.txt
## SRR1039520_1.5Msample.fastq.gz
## SRR1039520_2.5Msample.fastq.gz
## SRR1039521_1.5Msample.fastq.gz
## SRR1039521_2.5Msample.fastq.gz
## SRR1039520_1.5Msample.fastq.gz: OK
## SRR1039520_2.5Msample.fastq.gz: OK
## SRR1039521_1.5Msample.fastq.gz: OK
## SRR1039521_2.5Msample.fastq.gz: OK
```


```bash
cd ..
```

## Use zcat to take a look into fastq.gz files (usually we do not do this)

Note the matching fastq IDs in the R1 and R2 files.


```bash
ls
ls fastqs
zcat fastqs/SRR1039520_1.5Msample.fastq.gz | head
zcat fastqs/SRR1039520_2.5Msample.fastq.gz | head
```

```
## 01_linux_basics.md
## 01_linux_basics.Rmd
## 02_hpc_intro.md
## 02_hpc_intro.Rmd
## 03_bioinfx_example.Rmd
## bbc_workshop_202209.rds
## _bookdown.yml
## docs
## fastqs
## full_fastqs
## index.md
## index.Rmd
## _output.yml
## packages.bib
## render3f7d925f9a22f.rds
## toc.css
## md5sum.txt
## SRR1039520_1.5Msample.fastq.gz
## SRR1039520_2.5Msample.fastq.gz
## SRR1039521_1.5Msample.fastq.gz
## SRR1039521_2.5Msample.fastq.gz
## 
## gzip: fastqs/SRR1039520_1.5Msample.fastq.gz: not in gzip format
## 
## gzip: fastqs/SRR1039520_2.5Msample.fastq.gz: not in gzip format
```

## Run fastqc or zcat {} | wc -l to count how many reads


```bash

```


## How to see what packages are installed on HPC. Use grep.


```bash
module av
module av bbc | head -n50

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
## bbc/CogentAP/CogentAP-1.0
## bbc/CogentAP/CogentAP-1.5
## bbc/csvtk/csvtk-v0.24.0
## bbc/curl/curl-7.65.1
## bbc/cutadapt/cutadapt-2.10
## bbc/cutadapt/cutadapt-2.3
## bbc/cutadapt/cutadapt-3.2
## bbc/cutadapt/default
## bbc/datamash/datamash-1.3
## bbc/deeptools/deeptools-3.3.1
## bbc/deeptools/deeptools-3.4.3
## bbc/deeptools/deeptools-3.5.1
## bbc/delta/delta-0.10.2
## bbc/dust/dust-0.7.5
## bbc/emboss/emboss-6.6.0
## bbc/fancplot/fancplot-0.9.1
## bbc/fastahack/fastahack-1.0
## bbc/fastp/fastp-0.21.0
## bbc/fastqc/fastqc-0.11.8
## bbc/fastqc/fastqc-0.11.9
## bbc/fastq_screen/fastq_screen-0.14.0
## bbc/fastx_toolkit/fastx_toolkit-0.0.13
## bbc/freebayes/freebayes-1.3.1
## bbc/GAG/GAG-2.0.1
## bbc/gatk/gatk-3.8-1-0-gf15c1c3ef
## bbc/gatk/gatk-4.1.4.1
## bbc/gatk/gatk-4.1.8.1
## bbc/gdc/gdc-1.6.0
## bbc/genometools/genometools-1.6.2
## bbc/gff3sort/gff3sort_dd881e1
## bbc/gffcompare/gffcompare-0.12.6
## bbc/gffread/gffread-0.12.7
## bbc/gnuplot/gnuplot-5.2.8
## bbc/gridss/gridss-2.7.3
## bbc/gsl/gsl-2.5
## bbc/gsutil/gsutil-4.52
## bbc/guppy/guppy-3.4.5
## bbc/hisat2/hisat2-2.1.0
## bbc/hisat2/hisat2-2.2.1
## bbc/HMMRATAC/HMMRATAC-1.2.10
## bbc/HMMRATAC/HMMRATAC-1.2.9
## bbc/homer/homer-4.11
## bbc/HOMER/HOMER-4.11.1
## bbc/htop/htop-2.2.0
## bbc/htseq/htseq-0.13.5
## bbc/htslib/htslib-1.10.2
## bbc/htslib/htslib-1.12
## bbc/htslib/htslib-1.14
## bbc/iaap-cli/iaap-cli
## bbc/idr/idr-2.0.4.2
## bbc/igblast/igblast-1.17.1
## bbc/jcvi/jcvi-1.1.18
## bbc/jellyfish/jellyfish-2.3.0
## bbc/kallisto/kallisto-0.46.1
## bbc/kb-python/kb-python-0.24.4
## bbc/last/last-1256
## bbc/libBigWig/libBigWig-0.4.6
## bbc/libgd/libgd-2.2.5
## bbc/libpng/png16
## bbc/lisa/lisa-1.0
## bbc/macs2/macs2-2.2.6
## bbc/macs2/macs2-2.2.7.1
## bbc/manta/manta-1.6.0
## bbc/meme/meme-5.1.1
## bbc/meme/meme-5.3.3
## bbc/minimap2/minimap2-2.17
## bbc/minimap2/minimap2-2.22
## bbc/mosdepth/mosdepth-0.2.6
## bbc/multiqc/multiqc-1.11
## bbc/multiqc/multiqc-1.12
## bbc/multiqc/multiqc-1.8
## bbc/multiqc/multiqc-1.9
## bbc/Muscle/muscle-5.1
## bbc/nanoplot/nanoplot-1.28.1
## bbc/nanoplot/nanoplot-1.28.2
## bbc/nanopolish/nanopolish-0.11.3
## bbc/ncdu/ncdu-1.15.1
## bbc/ncftp/ncftp-3.2.5
## bbc/nextflow/nextflow-21.10.6
## bbc/ont-fast5-api/ont-fast5-api-2.0.1
## bbc/openssl/openssl-1.1.1c
## bbc/orthofinder/orthofinder-2.5.4
## bbc/pandoc/pandoc-2.7.3
## bbc/parallel/parallel-20191122
## bbc/pcre2/pcre2-10.34
## bbc/PEPATAC/PEPATAC
## bbc/perl/perl-5.30.0
## bbc/picard/picard-2.21.4-SNAPSHOT
## bbc/picard/picard-2.23.3
## bbc/pigz/pigz-2.4
## bbc/plink/plink-v1.90b6.18
## bbc/preseq/preseq-2.0.3
## bbc/preseq/preseq-3.1.2
## bbc/pyGenomeTracks/pyGenomeTracks-3.2
## bbc/python2/python2.7.0
## bbc/python2/python-2.7.18
## bbc/python3/python-3.6.10
## bbc/python3/python-3.7.3
## bbc/python3/python-3.7.4
## bbc/python3/python-3.8.1
## bbc/python3/python-3.8.1-bz
## bbc/python3/python-3.8.1-bz-sqlite3
## bbc/qualimap/qualimap_v.2.2.2
## bbc/R/alt/R-4.2.1-setR_LIBS_USER
## bbc/R/R-3.6.0
## bbc/R/R-4.0.0-pcre2
## bbc/R/R-4.0.0-pcre2-setR_LIBS_USER
## bbc/R/R-4.0.0-USE_PCRE2
## bbc/R/R-4.0.2
## bbc/R/R-4.0.2-setR_LIBS_USER
## bbc/R/R-4.1.0
## bbc/R/R-4.1.0-setR_LIBS_USER
## bbc/R/R-4.2.1
## bbc/racon/racon-1.4.3
## bbc/readline/readline-8.0
## bbc/rmats_turbo/rmats_turbo-4.1.1
## bbc/rsat/rsat-2020.02.07
## bbc/RSeQC/RSeQC-4.0.0
## bbc/rust/rust-1.52.1
## bbc/salmon/salmon-1.2.0
## bbc/salmon/salmon-1.3.0
## bbc/salmon/salmon-1.4.0
## bbc/salmon/salmon-1.5.2
## bbc/sambamba/sambamba-0.7.1
## bbc/samblaster/samblaster-0.1.24
## bbc/samblaster/samblaster-0.1.26
## bbc/samclip/samclip-v0.4.0
## bbc/samtools/samtools-1.12
## bbc/samtools/samtools-1.14
## bbc/samtools/samtools-1.9
## bbc/seqkit/seqkit-v2.1.0
## bbc/seqtk/seqtk-1.3-r115-dirty
## bbc/shasta/shasta-0.4.0
## bbc/skewer/skewer
## bbc/snakemake/snakemake-5.10.0
## bbc/snakemake/snakemake-5.14.0
## bbc/snakemake/snakemake-5.15.0
## bbc/snakemake/snakemake-5.17.0
## bbc/snakemake/snakemake-5.19.0
## bbc/snakemake/snakemake-5.20.1
## bbc/snakemake/snakemake-5.23.0
## bbc/snakemake/snakemake-5.28.0
## bbc/snakemake/snakemake-5.32.2
## bbc/snakemake/snakemake-6.1.0
## bbc/snakemake/snakemake-6.13.1_test
## bbc/snakemake/snakemake-6.15.0
## bbc/snakemake/snakemake-7.8.5
## bbc/snakePipes/snakePipes-2.1.2
## bbc/sniffles/sniffles-1.0.11
## bbc/SnpEff/SnpEff-4.3t
## bbc/sortmerna/sortmerna-4.3.4
## bbc/spaceranger/spaceranger-1.2.1
## bbc/spaceranger/spaceranger-1.3.0
## bbc/sratoolkit/sratoolkit-2.10.0
## bbc/sratoolkit/sratoolkit-2.11.0
## bbc/STAR/STAR-2.7.10a
## bbc/STAR/STAR-2.7.1a
## bbc/STAR/STAR-2.7.3a
## bbc/STAR/STAR-2.7.8a
## bbc/stringtie/stringtie-2.1.6
## bbc/stringtie/stringtie-2.2.1
## bbc/subread/subread-2.0.0
## bbc/taxonkit/taxonkit-v0.9.0
## bbc/tbmate/1.6
## bbc/tealdeer/tealdeer-1.4.1
## bbc/TEtranscripts/TEtranscripts-2.0.4
## bbc/TEtranscripts/TEtranscripts-2.2.1
## bbc/texlive/texlive-20190625
## bbc/texlive/texlive-20210928
## bbc/transdecoder/transdecoder-5.5.0
## bbc/tree/tree-1.8.0
## bbc/trim_galore/trim_galore-0.6.0
## bbc/Trimmomatic/trimmomatic-0.39
## bbc/Trinity/Trinity-2.13.0
## bbc/TRUST4/trust4-1.0.2b
## bbc/ucsc/ucsc-2020.06.11
## bbc/UMI-Tools/UMI-Tools-1.1.1
## bbc/vcflib/vcflib-1.0.1
## bbc/vcftools/vcftools-0.1.16
## bbc/vsearch/vsearch-2.21.1
## bbc/vt/vt-0.1.16
## bbc/WiggleTools/WiggleTools-1.2.11
## bbc/xlsx2csv/xlsx2csv-0.7.8
## bbc/Zzz_deprecated/cairo-1.15.12
## bbc/Zzz_deprecated/cairo-1.16.0
## bbc/Zzz_deprecated/multiqc-1.8
## bbc/Zzz_deprecated/snakemake/snakemake-5.8.2
## bbc/Zzz_deprecated/zlib/zlib-1.2.11
## blacs/openmpi/gcc/64/1.1patch03
## blacs/openmpi/open64/64/1.1patch03
## blas/gcc/64/3.6.0
## blas/open64/64/3.6.0
## bonnie++/1.97.1
## cmgui/7.3
## cuda10.0/blas/10.0.130
## cuda10.0/fft/10.0.130
## cuda10.0/nsight/10.0.130
## cuda10.0/profiler/10.0.130
## cuda10.0/toolkit/10.0.130
## cuda70/blas/7.0.28
## cuda70/fft/7.0.28
## cuda70/gdk/346.46
## cuda70/nsight/7.0.28
## cuda70/profiler/7.0.28
## cuda70/toolkit/7.0.28
## cuda80/blas/8.0.44
## cuda80/fft/8.0.44
## cuda80/nsight/8.0.44
## cuda80/profiler/8.0.44
## cuda80/toolkit/8.0.44
## cuda91/blas/9.1.85
## cuda91/fft/9.1.85
## cuda91/nsight/9.1.85
## cuda91/profiler/9.1.85
## cuda91/toolkit/9.1.85
## default-environment
## fftw2/openmpi/gcc/64/double/2.1.5
## fftw2/openmpi/gcc/64/float/2.1.5
## fftw2/openmpi/open64/64/double/2.1.5
## fftw2/openmpi/open64/64/float/2.1.5
## fftw3/openmpi/gcc/64/3.3.4
## fftw3/openmpi/open64/64/3.3.4
## gcc/6.5.0
## gcc/7.5.0
## gdb/7.11
## globalarrays/openmpi/gcc/64/5.4
## globalarrays/openmpi/open64/64/5.4
## gpfs/4.0
## hdf5/1.6.10
## hdf5_18/1.8.17
## hpcx/2.0.0
## hpl/2.2
## hwloc/1.11.3
## intel-cluster-checker/2.2.2
## intel-cluster-runtime/ia32/3.8
## intel-cluster-runtime/intel64/3.8
## intel-cluster-runtime/mic/3.8
## intel-tbb-oss/ia32/2017_20161128oss
## intel-tbb-oss/intel64/2017_20161128oss
## iozone/3_434
## java/11.0.5
## java/13.0.1
## java/1.8.0_60
## lapack/gcc/64/3.6.0
## lapack/open64/64/3.6.0
## maui/3.3.1
## moab/latest
## mpich/ge/gcc/64/3.1.4
## mpich/ge/open64/64/3.1.4
## mpiexec/0.84_432
## mvapich/gcc/64/1.2rc1
## mvapich/open64/64/1.2rc1
## mvapich2/gcc/64/2.2rc1
## mvapich2/open64/64/2.2rc1
## nccl/1.3.4
## netcdf/gcc/64/4.4.0
## netcdf/open64/64/4.4.0
## netperf/2.7.0
## open64/4.5.2.1
## openblas/dynamic/0.2.18
## openlava/3.0
## openmpi/cuda/64/3.0.0
## openmpi/gcc/64/1.10.1
## openmpi/open64/64/1.10.1
## puppet/3.7.5
## scalapack/mvapich2/gcc/64/2.0.2
## scalapack/openmpi/gcc/64/2.0.2
## sge/2011.11p1
## slurm/16.05.8
## torque/6.0.2
## torque/6.1.3
## VARI/bisulfite_seq
## VARI/bsoft
## VARI/ccache
## VARI/cryoem
## VARI/emhp
## VARI/genomics
## VARI/git-2.9.5
## VARI/macs2
## VARI/maxquant
## VARI/motioncor2
## VARI/p7zip-16.02
## VARI/pbsPretty
## VARI/php-5.6.23
## VARI/pyem
## VARI/relion1.4
## VARI/relion2
## VARI/relion2.1
## VARI/relion3
## VARI/relion3.0.3
## VARI/relion3.1.0
## VARI/relion3.1.0-4.1.0
## VARI/relion3.1b
## VARI/relion4b
## VARI/singularity
## VARI/topaz
## VARI/topaz0.2.3
## VARI/tree-1.8.0
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
## bbc/CogentAP/CogentAP-1.0
## bbc/CogentAP/CogentAP-1.5
## bbc/csvtk/csvtk-v0.24.0
## bbc/curl/curl-7.65.1
## bbc/cutadapt/cutadapt-2.10
## bbc/cutadapt/cutadapt-2.3
## bbc/cutadapt/cutadapt-3.2
## bbc/cutadapt/default
## bbc/datamash/datamash-1.3
## bbc/deeptools/deeptools-3.3.1
## bbc/deeptools/deeptools-3.4.3
## bbc/deeptools/deeptools-3.5.1
## bbc/delta/delta-0.10.2
## bbc/dust/dust-0.7.5
## bbc/emboss/emboss-6.6.0
## bbc/fancplot/fancplot-0.9.1
## bbc/fastahack/fastahack-1.0
## bbc/fastp/fastp-0.21.0
## bbc/fastqc/fastqc-0.11.8
## bbc/fastqc/fastqc-0.11.9
## bbc/fastq_screen/fastq_screen-0.14.0
## bbc/fastx_toolkit/fastx_toolkit-0.0.13
## bbc/freebayes/freebayes-1.3.1
## bbc/GAG/GAG-2.0.1
## bbc/gatk/gatk-3.8-1-0-gf15c1c3ef
## bbc/gatk/gatk-4.1.4.1
## bbc/gatk/gatk-4.1.8.1
## bbc/gdc/gdc-1.6.0
## bbc/genometools/genometools-1.6.2
## bbc/gff3sort/gff3sort_dd881e1
## bbc/gffcompare/gffcompare-0.12.6
## bbc/gffread/gffread-0.12.7
## bbc/gnuplot/gnuplot-5.2.8
## bbc/gridss/gridss-2.7.3
## bbc/gsl/gsl-2.5
## bbc/gsutil/gsutil-4.52
## bbc/guppy/guppy-3.4.5
## bbc/hisat2/hisat2-2.1.0
## bbc/hisat2/hisat2-2.2.1
## bbc/HMMRATAC/HMMRATAC-1.2.10
## bbc/HMMRATAC/HMMRATAC-1.2.9
## bbc/homer/homer-4.11
## bbc/HOMER/HOMER-4.11.1
## bbc/htop/htop-2.2.0
## bbc/htseq/htseq-0.13.5
## bbc/htslib/htslib-1.10.2
## bbc/htslib/htslib-1.12
## bbc/htslib/htslib-1.14
## bbc/iaap-cli/iaap-cli
## bbc/idr/idr-2.0.4.2
## bbc/igblast/igblast-1.17.1
## bbc/jcvi/jcvi-1.1.18
## bbc/jellyfish/jellyfish-2.3.0
## bbc/kallisto/kallisto-0.46.1
## bbc/kb-python/kb-python-0.24.4
## bbc/last/last-1256
## bbc/libBigWig/libBigWig-0.4.6
## bbc/libgd/libgd-2.2.5
## bbc/libpng/png16
## bbc/lisa/lisa-1.0
## bbc/macs2/macs2-2.2.6
## bbc/macs2/macs2-2.2.7.1
## bbc/manta/manta-1.6.0
## bbc/meme/meme-5.1.1
## bbc/meme/meme-5.3.3
## bbc/minimap2/minimap2-2.17
## bbc/minimap2/minimap2-2.22
## bbc/mosdepth/mosdepth-0.2.6
## bbc/multiqc/multiqc-1.11
## bbc/multiqc/multiqc-1.12
## bbc/multiqc/multiqc-1.8
## bbc/multiqc/multiqc-1.9
## bbc/Muscle/muscle-5.1
## bbc/nanoplot/nanoplot-1.28.1
## bbc/nanoplot/nanoplot-1.28.2
## bbc/nanopolish/nanopolish-0.11.3
## bbc/ncdu/ncdu-1.15.1
## bbc/ncftp/ncftp-3.2.5
## bbc/nextflow/nextflow-21.10.6
## bbc/ont-fast5-api/ont-fast5-api-2.0.1
## bbc/openssl/openssl-1.1.1c
## bbc/orthofinder/orthofinder-2.5.4
## bbc/pandoc/pandoc-2.7.3
## bbc/parallel/parallel-20191122
## bbc/pcre2/pcre2-10.34
## bbc/PEPATAC/PEPATAC
## bbc/perl/perl-5.30.0
## bbc/picard/picard-2.21.4-SNAPSHOT
## bbc/picard/picard-2.23.3
## bbc/pigz/pigz-2.4
## bbc/plink/plink-v1.90b6.18
## bbc/preseq/preseq-2.0.3
## bbc/preseq/preseq-3.1.2
## bbc/pyGenomeTracks/pyGenomeTracks-3.2
## bbc/python2/python2.7.0
## bbc/python2/python-2.7.18
## bbc/python3/python-3.6.10
## bbc/python3/python-3.7.3
## bbc/python3/python-3.7.4
## bbc/python3/python-3.8.1
## bbc/python3/python-3.8.1-bz
## bbc/python3/python-3.8.1-bz-sqlite3
## bbc/qualimap/qualimap_v.2.2.2
## bbc/R/alt/R-4.2.1-setR_LIBS_USER
## bbc/R/R-3.6.0
## bbc/R/R-4.0.0-pcre2
## bbc/R/R-4.0.0-pcre2-setR_LIBS_USER
## bbc/R/R-4.0.0-USE_PCRE2
## bbc/R/R-4.0.2
## bbc/R/R-4.0.2-setR_LIBS_USER
## bbc/R/R-4.1.0
## bbc/R/R-4.1.0-setR_LIBS_USER
## bbc/R/R-4.2.1
## bbc/racon/racon-1.4.3
## bbc/readline/readline-8.0
## bbc/rmats_turbo/rmats_turbo-4.1.1
## bbc/rsat/rsat-2020.02.07
## bbc/RSeQC/RSeQC-4.0.0
## bbc/rust/rust-1.52.1
## bbc/salmon/salmon-1.2.0
## bbc/salmon/salmon-1.3.0
## bbc/salmon/salmon-1.4.0
## bbc/salmon/salmon-1.5.2
## bbc/sambamba/sambamba-0.7.1
## bbc/samblaster/samblaster-0.1.24
## bbc/samblaster/samblaster-0.1.26
## bbc/samclip/samclip-v0.4.0
## bbc/samtools/samtools-1.12
## bbc/samtools/samtools-1.14
## bbc/samtools/samtools-1.9
## bbc/seqkit/seqkit-v2.1.0
## bbc/seqtk/seqtk-1.3-r115-dirty
## bbc/shasta/shasta-0.4.0
## bbc/skewer/skewer
## bbc/snakemake/snakemake-5.10.0
## bbc/snakemake/snakemake-5.14.0
## bbc/snakemake/snakemake-5.15.0
## bbc/snakemake/snakemake-5.17.0
## bbc/snakemake/snakemake-5.19.0
## bbc/snakemake/snakemake-5.20.1
## bbc/snakemake/snakemake-5.23.0
## bbc/snakemake/snakemake-5.28.0
## bbc/snakemake/snakemake-5.32.2
## bbc/snakemake/snakemake-6.1.0
## bbc/snakemake/snakemake-6.13.1_test
## bbc/snakemake/snakemake-6.15.0
## bbc/snakemake/snakemake-7.8.5
## bbc/snakePipes/snakePipes-2.1.2
## bbc/sniffles/sniffles-1.0.11
## bbc/SnpEff/SnpEff-4.3t
## bbc/sortmerna/sortmerna-4.3.4
## bbc/spaceranger/spaceranger-1.2.1
## bbc/spaceranger/spaceranger-1.3.0
## bbc/sratoolkit/sratoolkit-2.10.0
## bbc/sratoolkit/sratoolkit-2.11.0
## bbc/STAR/STAR-2.7.10a
## bbc/STAR/STAR-2.7.1a
## bbc/STAR/STAR-2.7.3a
## bbc/STAR/STAR-2.7.8a
## bbc/stringtie/stringtie-2.1.6
## bbc/stringtie/stringtie-2.2.1
## bbc/subread/subread-2.0.0
## bbc/taxonkit/taxonkit-v0.9.0
## bbc/tbmate/1.6
## bbc/tealdeer/tealdeer-1.4.1
## bbc/TEtranscripts/TEtranscripts-2.0.4
## bbc/TEtranscripts/TEtranscripts-2.2.1
## bbc/texlive/texlive-20190625
## bbc/texlive/texlive-20210928
## bbc/transdecoder/transdecoder-5.5.0
## bbc/tree/tree-1.8.0
## bbc/trim_galore/trim_galore-0.6.0
## bbc/Trimmomatic/trimmomatic-0.39
## bbc/Trinity/Trinity-2.13.0
## bbc/TRUST4/trust4-1.0.2b
## bbc/ucsc/ucsc-2020.06.11
## bbc/UMI-Tools/UMI-Tools-1.1.1
## bbc/vcflib/vcflib-1.0.1
## bbc/vcftools/vcftools-0.1.16
## bbc/vsearch/vsearch-2.21.1
## bbc/vt/vt-0.1.16
## bbc/WiggleTools/WiggleTools-1.2.11
## bbc/xlsx2csv/xlsx2csv-0.7.8
## bbc/Zzz_deprecated/cairo-1.15.12
## bbc/Zzz_deprecated/cairo-1.16.0
## bbc/Zzz_deprecated/multiqc-1.8
## bbc/Zzz_deprecated/snakemake/snakemake-5.8.2
## bbc/Zzz_deprecated/zlib/zlib-1.2.11
```

## Set up a job to run Salmon

The job script code is below. Don't worry about the specifics of this code for now. Simply, copy the code and paste into your text editor to make your job script. You will need to edit certain fields.


```bash

```

## Use qstat to check on the job


```bash
qstat
qstat | grep 'firstname.lastname'

```

## Check the job logs to see if job finished running


```bash
cat job.log

```

## Use grep to find the TPMs for specific genes


```bash

grep 'SNC' salmon.results.txt
```

## Transfer read counts table to laptop/desktop to view in Excel


```bash

```


