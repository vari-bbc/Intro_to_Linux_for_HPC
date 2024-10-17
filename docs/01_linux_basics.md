---
title: "Basic Command Line Utilities on the HPC"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography: [packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: vari-bbc/Intro_to_Linux_for_HPC
theme: "yeti"
---




# **Basic Command Line Utilities on the HPC**


## **Access HPC**

Ensure that your laptop is connected to wifi "vai", not "vai-guest". 

<details>
<summary>**archived content:** previous ssh connection instructions</summary>

<blockquote>
**Note:**
The platform-specific instructions below are **hidden by default** because of the 
newer HPC OnDemand platform, which we can all access from our laptops using a 
web browser.

However, for archival purposes, we will still provide the platform-specific 
instructions to access the HPC on your laptop.

#### **Mac users**

1. Open your terminal (search for the app)
2. Type the following command to connect to the HPC.
*Replace `firstname.lastname` with your actual VAI username.*
```bash
ssh firstname.lastname@access.hpc.vai.org
```
3. Enter your VAI password. (Note: you will not see the password as you type it,
but it *is* working! Press your enter/return key when you're done.)


#### **Windows users**

If MobaXterm is not installed, please download the [free MobaXterm version](https://mobaxterm.mobatek.net/download.html), then set up, and agree with license agreements to proceed installation.

1. Open MobaXterm
2. Click on `Session` in the top left corner
3. Click on `SSH`(Secure Shell)
4. Enter `access.hpc.vai.org` in the Remote host field
5. Click the box on the left of Specify username, then enter your VAI `username` and `password`
6. Click `OK`
</blockquote>
</details>

### **OnDemand Instructions**

OnDemand is a web-based portal that lets you access the HPC resources right from your browser, without needing VPN or ssh tools on your local computer. 
You can read more about this on the [HPC OnDemand sharepoint page](https://vanandelinstitute.sharepoint.com/sites/SC/SitePages/HPC-OnDemand.aspx).

Below are the steps to access the HPC command line using OnDemand:

1. Open a web browser (Chrome, Firefox, Safari, etc.)
2. Go to the following URL: [https://ondemand3.vai.zone/](https://ondemand3.vai.zone/)
    - if you see an error, try this instead: [https://ondemandlocal.hpc.vai.org](https://ondemandlocal.hpc.vai.org/) 
    - if neither or those links work, please let us know!
3. Select "hpc Shell Access" from the list of several "pinned apps". ("hpc Shell Access" is similar to the terminal on your laptop)

Congrats -- you've successfully connected to the HPC! 


## **File navigation**

Practice directory: `/varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II`. 

### **Navigate to a directory (folder)**

`cd` - Change Directory. Allows you to navigate to a different directory (folder). 
Note the following special symbols:

- `.` refers to the folder you're currently in (your "working" directory).
- `..` refers to the folder before where you currently are - this is your "parent" directory. 
- `/` is the root directory.
- `~` is your personal home directory.

`pwd` - Print Working Directory. Displays the current directory (folder) you are in.





```bash
# Display the current directory
pwd

# Change to the practice directory
cd /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II

# Display the current directory again to confirm the change
pwd
```

```
## /home/kaitlyn.denhaan
## /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II
```




### **List content in a directory**

- `ls`  - list contents (files and folders). Without anything specified after `ls`, it will list what's in the current directory. 
  - `ls -lht` – list contents, with added options (`-l` more details in a *long* list, `-h` *h*uman-readable sizes, `t` sorted by modification *t*ime).  
  - `ls foldername` - list contents of the specified folder. 

Shows more details about the files and folders, 

```bash
ls
ls -lht
ls -lht /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II
ls -lht ~
```


### **Navigate back to home directory**

Without anything after `cd`, you will go to your home directory, which is `/home/username` (same as `~`).


```bash
# same as "cd ~" or "cd /home/username"
cd

# check: where are you now?
pwd
```


### **Create a directory**

`mkdir dir_name` - will create a directory "dir_name" in your current directory. 
You will be creating a directory in your home drectory (because we just navigated to it).

Make sure you are in your home directory first (use pwd), then create a directory called "hpc_workshop_2024"


```bash
pwd
mkdir hpc_workshop_2024
```


### **Copy a file**

`cp` - copy file/files

Copy a file from `/varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II` to the directory you just created.


```bash
pwd
cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/metadata.tsv hpc_workshop_2024
ls hpc_workshop_2024
```


<!-- What if we want to copy *multiple* files? -->
<!-- There are a few different ways of doing this, but we'll use `*`.  -->
<!-- This symbol (an asterisk) can be thought of as a "select all" symbol. -->

<!-- To copy all the files in this folder, we can use the following command: -->

<!-- ```{bash, eval=FALSE, engine="sh"} -->
<!-- cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/* hpc_workshop_2024 -->
<!-- ls hpc_workshop_2024 -->
<!-- ``` -->


## **View and manipulate files**

### **Display content**

`cat` - display the entire contents of a file.  
`head` - display the *first* 10 lines/rows of the a file.  
`tail` - display the *last* 10 lines/rows of the a file. You can use "-n" to control how many lines you want to see, default is 10.


```bash
cd ~/hpc_workshop_2024

cat metadata.tsv
head metadata.tsv
tail metadata.tsv
head -n 2 metadata.tsv
```

The `cat` command can also be used to combine files. (This is where the command's name comes from: con*cat*enate)
Note that we won't run the command below, but keep this in mind for future reference.


```bash
cat file1.txt file2.txt > combined.txt
```



### **Pattern search**

`grep` - search a pattern by line


```bash
pwd
grep "13" metadata.tsv 
```

```
## /home/kaitlyn.denhaan/hpc_workshop_2024
## M	13	sample 13	replicate a
## Z	26	sample 13	replicate b
```


### **Display the number of words, lines, and characters**

`wc` - word count. It counts the number of words, lines, and characters in a file.

Adding the `-l` option afterwards just gives the number of ***l**ines* in the file.


```bash
wc metadata.tsv
wc -l metadata.tsv
```

```
##  26 156 675 metadata.tsv
## 26 metadata.tsv
```


### **Pipe - redirection**

A pipe is a form of redirection (instead of printing output to the screen, 
it sends it to other destinations). You can send output from one command/program
to another for further processing, such as: 
`command 1 | command 2 | command 3`.

- Above, the output from `command 1` is used as input for `command 2`, and the 
output from `command 2` is used as input for `command 3`.

In the example below, we will use 3 commands subsequently to count the number of lines that contain the number "13" in the file.


```bash
cat metadata.tsv | grep "13" | wc -l
```

```
## 2
```

### **Output redirection**

Instead of printing output to your screen (typical command output) or 
another command (pipe), you can redirect the output to a file.

`>` - redirect output to a file. Note that it will **overwrite** the file if it already exists -- be careful!  
`>>` - append output to a file. It will **add the output to the end** of the file.


```bash
# note that no output will be displayed on the screen -- it's saved in the file instead
grep "13" metadata.tsv > lines_with_13.tsv

# check: do we see the file we just created?
ls

# display the content of the file
cat lines_with_13.tsv
```

```
## combined.fq
## data_01_R1.fq
## data_54_R1.fq
## lines_with_13.tsv
## metadata.tsv
## M	13	sample 13	replicate a
## Z	26	sample 13	replicate b
```

## **Exercise**

Below, we'll be doing something similar to the exercises above, but with real genomic data!
We'll be using a fastq file (format containing raw sequencing & quality information).

We will:

1. copy fastq files from the practice directory to the folder you created
2. combine the two fastq files into one
3. count the number of reads in the combined file

See if you can do each step on your own -- if you get stuck, don't worry! 
Try to remember the commands we've learned so far, and you can always **refer back to the examples above.**
To check your work, you can view the commands below.

### **Step 1: Copy the files**

Once you're in your `hpc_workshop_2024` directory, copy the files `data_01_R1.fq` and `data_54_R1.fq`
from the practice directory (`/varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II`) into your folder.

<details>
<summary>Click here to see a solution</summary>


```bash
cd ~/hpc_workshop_2024
cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/data_01_R1.fq .
cp /varidata/researchtemp/hpctmp/BBC_workshop_Oct2024_II/data_54_R1.fq .
```

</details>

### **Step 2: Combine the files**

Combine the two fastq files into one file called `combined.fq`.

<details>
<summary>Click here to see a solution</summary>


```bash
cat data_01_R1.fq data_54_R1.fq > combined.fq
ls
```

```
## combined.fq
## data_01_R1.fq
## data_54_R1.fq
## lines_with_13.tsv
## metadata.tsv
```

</details>

### **Step 3: Count the number of reads**

Count the number of reads in the combined file.

<details>
<summary>Hint #1 -- what is a "read"?</summary>
In a fastq file, reads start with the "@" symbol.
</details>


<details>
<summary>Hint #2 -- how to approach this?</summary>
See if you can **combine** a few commands together! Plan out each of the steps you need to do, then put them all together. 
</details>
<br>

<details>
<summary>Click here to see an example solution</summary>


```bash
# use that combined file
grep "@" combined.fq | wc -l

# a different way to do the same thing
cat combined.fq | grep "@" | wc -l

# or, combine the commands all together
cat data_01_R1.fq data_54_R1.fq | grep "@" | wc -l
```

```
## 500
## 500
## 500
```

</details>
<br>

Can you explain why your solution works? Which solution above does yours look most similar to?

Congratulations! You've successfully completed the exercise.

**BONUS:** If you have extra time, use the skills you've learned to explore the fastq file format further:

1. Count the total number of **lines** in the combined file.
2. Using your answer from the exercise above (count the number of **reads**) -- how are the number of reads related to the total number of lines in the file?
3. Look at some of the contents of the combined fastq file. What do you notice about the structure of the file? Can you find any patterns? What do you think each line represents?  
4. What are some questions you have about the fastq file format?  

## **Summary**

In this section, we've covered the basic commands for navigating directories, viewing and manipulating files, using pipes, and redirecting output.

You've learned how to:

- Navigate to a directory
- List content in a directory
- Create a directory
- Copy a file
- Display content
- Search for patterns
- Display the number of words, lines, and characters
- Use pipes for redirection
- Use output redirection

These are the fundamental commands you'll need to work with files and directories on the HPC.

In the next section, we'll work through a real bioinformatics "mini-project" that builds on these commands, so you can see how they can be used together to solve a problem.

[Next: Bioinformatics Mini Project](mini-project.html)
