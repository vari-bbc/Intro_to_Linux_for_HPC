---
title: "Basic command lines for HPC"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography: [packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: vari-bbc/Intro_to_Linux_for_HPC
theme: "yeti"
---



# **Basic command lines for HPC**


## **Access HPC**

Ensure that your laptop is connected to wifi "vai", not "vai-guest". 


### **Mac users**

Open your terminal, and type **ssh vai_username@access.hpc.vai.org**, and enter your VAI password.


### **Windows user**

If MobaXterm is not installed, please download the [free MobaXterm version](https://mobaxterm.mobatek.net/download.html), then set up, and agree with license agreements to proceed installation.

- Open MobaXterm
- Click on `Session` in the top left corner
- Click on `SSH`(Secure Shell)
- Enter `access.hpc.vai.org` in Remote host
- Click the little box on the left of Specify username, then type VAI `username` and enter VAI `password`
- Next, click `OK`


## **File navigation**

Practice directory: `/varidata/researchtemp/hpctmp/BBC_workshop_June2023_I`. 

### **Navigate to a directory**
`cd` - change directory, navigate to a directory. Tip, ./ means current working directory. ../ means the parent directory. / is the root directory, ~ is home directory. Pay attention to these differences.


```bash
cd /varidata/researchtemp/hpctmp/BBC_workshop_June2023_I
```


### **Current working directory**
`pwd` - displays the current working directory


```bash
pwd
```

```
## /varidata/research/projects/bbc/research/hpc_workshop_202209
```


### **List content in a directory**
`ls`  - list content; Without anything specified after ls, it will list the current directory. It can also list content of another directory. 
`ls -lht` – will show more details of content

```bash
ls
ls -lht
ls -lht /varidata/researchtemp/hpctmp/BBC_workshop_June2023_I
```


### **Navigate to home irectory**
Without anything after `cd`, you will change to your home directory, which is `/home/username`

```bash
cd
pwd
```


### **Create a directory**
`mkdir dir_name` - will create a directory "dir_name" in your current directory. You will be creating a directory in your home drectory. Make sure you are in your home directory first (use pwd), then create a directory called "hpc_workshop"

```bash
pwd
mkdir hpc_workshop
```


### **Copy a file**
`cp` - copy file/files
Copy a file from `/varidata/researchtemp/hpctmp/BBC_workshop_June2023_I` to the directory you just created.

```bash
pwd
cp /varidata/researchtemp/hpctmp/BBC_workshop_June2023_I/test_01_R1.fq hpc_workshop
ls hpc_workshop
```


## **View and manipulate files**

### **Display content**
`head` - display as standard output (stdout) the first 10 lines/rows of the a file.
`tail` - display as standard output (stdout) the last 10 lines/rows of the a file. You can use "-n" to control how many lines you want to see, default is 10

```bash
cd hpc_workshop
head test_01_R1.fq
tail test_01_R1.fq
head -2 test_01_R1.fq
```

```
## bash: line 1: cd: hpc_workshop: No such file or directory
## @A00426:207:H537LDMXY:1:1101:2067:1000 1:N:0:NCTAAGAT+NCGCGGTT
## NTAACAGTGACTTGCGGGGGAAGTCTACGCGCGTGTGCACGCGGCACTCTC
## +
## #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
## @A00426:207:H537LDMXY:1:1101:2139:1000 1:N:0:NCTAAGAT+NCGCGGTT
## NTGGCCATTCACAGTATGGTATTTCTGAATAACAATCTTATCCACAGAGTC
## +
## #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
## @A00426:207:H537LDMXY:1:1101:2230:1000 1:N:0:NCTAAGAT+NCGCGGTT
## NAGACTAATCATCAGATCTCCTCTCTCTATGCTACATCCACTCCATTCAA
## +
## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
## @A00426:207:H537LDMXY:1:1101:7148:1063 1:N:0:NCTAAGAT+NCGCGGTT
## TGGCCTTGCTCACAGAGCTGCGTGAGAAACAGACGGTGCTTGCGATCTCTG
## +
## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
## @A00426:207:H537LDMXY:1:1101:7256:1063 1:N:0:NCTAAGAT+NCGCGGTT
## GCACGGGCGAGGGCGGGAACGGCGGAGCGGGAAGAAGCCGCGAGCGCGGAT
## +
## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
## @A00426:207:H537LDMXY:1:1101:2067:1000 1:N:0:NCTAAGAT+NCGCGGTT
## NTAACAGTGACTTGCGGGGGAAGTCTACGCGCGTGTGCACGCGGCACTCTC
```


### **Pattern search**
`grep` - search a pattern by line

```bash
grep "AATTGG" test_01_R1.fq
```

```
## NTCTGAATTGGGTTATGAGGCCCGGGAGGTGCCTCACCTCAGCCATTGAAC
## TGCCATTCTGTGCTCTCAGGACCTCTAATTGGGGGCCGTGGCAAAGGAGTG
```


### **Display the number of words, lines, and characters**

```bash
wc test_01_R1.fq
wc -l test_01_R1.fq
```

```
##  1000  1250 42064 test_01_R1.fq
## 1000 test_01_R1.fq
```


### **List the content as stdout**
`cat` – list the contents of a file as stdout; cat (short for concatenate) is one of the most frequently used commands. To run this command, type cat followed by the file’s name: `cat file.txt`

Combine two files: `cat file1 file2`.
It is useful to pipe the content of the file into another command.


### **Pipe - redirection**
A pipe is a form of redirection (transfer of sdout to other destinations). You can send output from one command/program to another for further processing, such as `command 1| command 2 | command 3`.

```bash
cat test_01_R1.fq | grep "@" | wc -l 
```

```
## 250
```


## **Exercise**
Below we will copy another fowward reads into your home directory (~ means home directory), and combine them into one and count how many reads were there. You can direct the stdout of any command to a new file use >.

```bash
cd ~/hpc_workshop
cp /varidata/researchtemp/hpctmp/BBC_workshop_June2023_I/test_54_R1.fq .
cat test_01_R1.fq test_54_R1.fq | grep "@" | wc -l 
cat test_01_R1.fq test_54_R1.fq > combined.fq
ls
```

```
## bash: line 1: cd: /home/kin.lau/hpc_workshop: No such file or directory
## 500
## 01_linux_basics.Rmd
## 02_bioinfx_example.Rmd
## 999_glossary.Rmd
## _bookdown.yml
## _output.yml
## bbc_bioinfx_book.rds
## combined.fq
## docs
## fastqs
## full_fastqs
## index.Rmd
## index.md
## multiqc
## packages.bib
## render1571c22d523052.rds
## run_salmon.sh
## salmon
## salmon.e
## salmon.o
## test_01_R1.fq
## test_54_R1.fq
## toc.css
```
