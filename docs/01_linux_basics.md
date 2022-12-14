

# Basics of using the HPC and the Linux commandline

Here we will learn the basics of using the Linux commandline on the HPC.

## Ex 1 Log into the HPC

**Step 1.** Make sure you are on the VAI network, not VAI guest. 

**Step 2.** 
- For Mac user: open terminal, type `ssh you_vai_user_name@submit.hpc.vai.org`, and type your password. 
- For Windows user, open MobarXterm, click new session under “Session”.
Click on “SSH”, type “submit.hpc.vai.org” as Remote host name, your user name for  VAI as username, click “OK”. You need to type your VAI password. 


## Commands you need to learn for Ex2:

**Cmd1.** `pwd` - displays the current working directory

It is useful when directory changes are made frequently.

**Cmd2.** `ls` - list content 

This command lists directory contents. Without anything specified after ls, it will list the current directory. It can also list content of another directory. You might not have any file in your home directory at this point, you can type `touch file.txt` to create your first file.

`ls /varidata/`

`ls -lht`  -- show more details of the files. 

**Cmd3.**  `mkdir` -- create a new directory

Syntax is mkdir new_directory name. You can create a new directory in the current directory or another directory. For example, `mkdir projects` will create a new directory called "projects" in your home directory.

**Cmd4.** `mv` -- move and rename file

The arguments in `mv` are similar to the `cp` command. 

To move a file, syntax is mv file_name destination’s directory. For example: `mv file.txt /home/username/projects`, will move `file.txt` from current directory to `/home/username/projects`

To rename a file, the syntax is `mv old_file_name.txt new_file_name.txt`

**Cmd5.** `cd` -- change directory

Let’s say you’re in `/home/username/` and you want to go to `projects`, a subdirectory of documents, type: `cd projects`

Another scenario is if you want to switch to a completely new directory, for example, `/home/username/research`, you need to type `cd` followed by the directory’s absolute path: `cd /home/username/research`

_Tip_, `./` means current working directory. `../` means the parent directory. `/` is the root directory, `~` is home directory. Pay attention to these differences.  

**Cmd6.** `rm` -- remove a file 

`rm file_name` will delete the file. 

`rm -r dir` will delete the whole directory. 

Note, `rm` is permanent. Unlike Window and Mac where deleted files can be recovered from the Trash. Think carefully before you type `rm`. 

**Cmd7.** `cp` -- copy a file or a dir 

Syntax: `cp old_file new_file` for file copy. 
`cp -r old_dir new_dir` for directory copy. 

**Cmd8.** `find` -- look for a file based on its name 

Typical usage: `find /home/username -name notes.txt` will search for a file called `notes.txt` within your home directory and its subdirectories.

_Tip_, if you want to find all pdf files under the current directory, do `find . -name “*.pdf”`. Here, `*` is called wildcard, it will match any name. 


## Ex2: file navigation

1. Use the commands you just learned, go to

/varidata/researchtemp/hpctmp/HPC_mini_workshop/Ex2


2. List files in this the directory, how many files and directories did you see? Which file was created the earliest?


3. Make a new directory in your home directory called “hpc_mini_workshop”


4. Try to copy (not to move) one of the fastsq files to hpc_mini_workshop you just created in your home directory. 


5. Use “find” command to list all the fastq files (file ends with fq) in /varidata/researchtemp/hpctmp/HPC_mini_workshop/Ex2


## Commands you need to learn for Ex 3

**Cmd9.**  `head` & `tail` -- view the first and last a few rows/lines of a file

_Usage_: `head -n 5 file_name` will display as standard output (stdout) the first 5 lines/rows of the a file. Without -n, it will default display 10 lines. 

**Cmd10.** `less` & `more` -- view content of a file one page at a time. 

_Tip_, you can use `q` key to quit the viewing. 

**Cmd11.** `cat` -- list the contents of a file as stdout

`cat` (short for concatenate) is one of the most frequently used commands. 

To run this command, type `cat` followed by the file’s name: `cat file.txt`

Combine two files: `cat file1 file2`

It is useful to **pipe** the content of the file into another command. 


### So, what is pipe?

A pipe is a form of redirection (transfer of sdout to other destinations). You can send output from one command/program/process to another command/program/process for further processing. You can do so by using the pipe character `|`.

Lastly, you can direct the stdout of any command to a new file use `>`. 

For example:

`bwa mem genome.fa reads.fastq | samtools sort -o output.bam`

`cat protein.fa | head -5 > first_5_protein.fa`


**Cmd 12.** `grep` -- search a pattern by line 

`grep “GCGGA” sequence_file.fastq` will output lines with GCGGA in it. 

**Cmd 13.** `wc` - counting words, lines, and bytes in files

Without any specification, wc will output three numbers. 

The most useful is wc -l, it only outputs the line number. 

## Exercise 3:  File viewing and manipulation

6. Count how many lines in the fastq file you just copied to your home directory. 


7. Create a file that only has protein sequence name from file:

`/varidata/researchtemp/hpctmp/HPC_mini_workshop/Ex2/protein_seq.fa`

(please copy this file to your home directory first). 

_Hint_, use `grep` to look for pattern that is shared among all protein names; then use redirection to create a new file. 


8. Copy a different fastq file from `/varidata/researchtemp/hpctmp/HPC_mini_workshop/Ex2/` to your home directory, and combine two fq files into a new file called `combined.fq`.


