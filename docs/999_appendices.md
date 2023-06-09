

# **Appendices**

## **Workshop Powerpoint files**

[June 8, 2023; Summer workshop I](https://vanandelinstitute-my.sharepoint.com/personal/daisy_fu_vai_org/Documents/Workshop_Jun_2023/Workshop1_BBC_merged.pdf)

## **Bash cheatsheet**

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Name </th>
   <th style="text-align:left;"> Command Line </th>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:left;"> Example </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Print Working Directory </td>
   <td style="text-align:left;"> `pwd` </td>
   <td style="text-align:left;"> Displays the current working directory </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ pwd`**

Result:/home/username </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> List </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `ls` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Lists the files and directories in the current directory </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ ls`**

Result: It returns empty after the $
  symbol since nothing has been created. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> List More Detail </td>
   <td style="text-align:left;"> `ls -lht` </td>
   <td style="text-align:left;"> Display more details about the file </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ ls -lht`**

Result: Display file detail in the current director </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Make Directory </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `mkdir` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Creates a new directory </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ mkdir hpc_mini_workshop`**

Result: A hpc_mini_workshop folder is created. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Move </td>
   <td style="text-align:left;"> `mv` </td>
   <td style="text-align:left;"> Moves or renames files or directories </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ mv hpc_mini_workshop workshopTraining`**

Result: Now the hpc_mini_workshop directory is called workshopTraining </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Change Directory </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `cd` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Change to an existing directory </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ cd workshopTraining`**

Result: [username\@submit002 workshopTraining]$ Notice ~ was in the home directory, now in the workshopTraining directory. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Remove </td>
   <td style="text-align:left;"> `rm` </td>
   <td style="text-align:left;"> Deletes files and directories </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ rm -r TaskProject`**

*Note: -r means directory

Result: TaskProject is deleted </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Copy </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `cp` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Copies files or directories </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ cp -r Task1 Project`**

Result: Task1 directory has moved to the Project directory. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Search for File </td>
   <td style="text-align:left;"> `find` </td>
   <td style="text-align:left;"> Search for files and directories based on various criteria like name, size, and modification time </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ find Project`**

Result: Task1 will appear within the Project directory

Project
Project/Task1 </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Head </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `head` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Display at the beginning of a file </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ head -n5 file_name`**

Result: It will display the first 5 lines from the beginning of the file_name. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tail </td>
   <td style="text-align:left;"> `tail` </td>
   <td style="text-align:left;"> Display at the end of a file </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ tail -n5 file_name`**

Result: It will display the last 5 lines from the end of the file_name </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Less </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `less` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Load the necessary portion of a file </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ less file.txt`**

Result: The user is able to view a portion of the file.txt. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> More </td>
   <td style="text-align:left;"> `more` </td>
   <td style="text-align:left;"> Load the entire file </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ more file.txt`**

Result: The user is able to view the entire file.txt. </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Quit </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `q` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Stop viewing the current file </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> quit viewing the current file </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Concatenate </td>
   <td style="text-align:left;"> `cat` </td>
   <td style="text-align:left;"> Display the contents of a file </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ cat file.txt`**

Result: Display the contents of file.txt </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Search for Text </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `grep` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Search for a specific pattern of text within files </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ grep “GCGGA” sequence_file.fastq`**

Result: Display any GCGGA pattern in sequence_file.fastq </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Word, Line, and Character Count </td>
   <td style="text-align:left;"> `wc` </td>
   <td style="text-align:left;"> Display the number of words, lines, and characters in a file </td>
   <td style="text-align:left;"> **`[username\@submit002 ~]$ wc file.txt`**

Result: Display the number of words, lines, and characters in the file.txt. </td>
  </tr>
  <tr>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Touch </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> `touch` </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> Create a new empty file </td>
   <td style="text-align:left;background-color: #e9f5f8 !important;"> **`[username\@submit002 ~]$ touch exampleFile.txt`**

Result: The command line will display file details in
  the current directory </td>
  </tr>
</tbody>
</table>
