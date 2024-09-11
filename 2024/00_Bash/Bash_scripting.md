# Bash scripting

## Background Processes
Programms can take a while and, in the meantime, we do not have access to the command line until they have finished or we kill the program (Ctrl+c). To avoid that we work in backgroung.

```
screen #open new screen
Ctrl+a+d #close screen
Ctrl+d #kill screen
screen -ls #list of all screen attached
screen -r <screenidnumber> #open screen with that idnumber 
```
---
## Concatenate commands and programs
### Pipe | 
| connects the standard output of one process to the standard input of another

```
grep "word" file1 | sed 's/ /\t/g' | program1 > file2 #sed use as input the output of grep; when using pipe the "original" input has to be specified only in the first command
```
### semicolon ;
; concatenate different commands or programs sequentially
```
grep ">" file1.fasta >output1 ; grep "_" file2.fasta output2
```
### &&
Concatenate two programs so that program2 run only if program1 completed successfully

```
program1 input.txt > intermediate-results.txt && program2 intermediate-results.txt > results.txt
```
---
## Standard output and standard error

Many programs use a "standard output" for outputting data, we usually redirected the standard output in an output file using ">". A separate file, called "standard error" is needed for errors, warnings, and messages. We can redirect the standard error with "2>"
```
program1 file 2> program1.stderr > results.txt
```
---
## Variables
A variable in bash can be anything, a number, a character, a string of characters, a file, a folder.
Variables are in bash are indicated with $
```
var1="ciao"
echo $var1 #print ciao
echo "$var1" #print ciao
echo '$var1' #print var1
```
---

## Command substitution 
Command substitution runs a Unix command inline and returns the output as a string that can be used in another command.
```
echo "$(cat file.fasta)"
echo "There are $(grep -c '^>' input.fasta) entries in my FASTA file." # show the string "There are 416 entries in my FASTA file."
```
---
## For loop
```
for i in *.fasta; do echo $i; done
for i in *.fasta; do mv $i ${i::-5}”_for_trimming”; done
for i in *.fasta; do mv $i ${i:1:3}”.fasta” ; done
for i in *fasta; do sed ‘s/>Locus/>/’ > $i”_editname” ; done
for i in *fasta; do grep –c”>” $i ; done > counts
for i in *fasta; do program1 $i > “output_”$i; done
for i in */ ; do cd $i; cp *.fasta ../; cd ..; done
```
---
## Bash script

Bash scripts are indicated with the .sh extention (python scripts with .py, perl scripts with .pl). 

Create bash script with vi:

```
#!/bin/bash
mkdir test
for i in *fasta; do mv $i test; done
cd test
for i in *fasta; do program1 $i > $i"_output"; done
for i in *output; do grep ">" $i; done > list_of_sequences
```
We need now to make the .sh file executable
```
chmod 777 namescript.sh
```
To execute the bash script
```
bash namescript.sh
```
## Bash script with variables

fake_script.sh
```
#!/bin/bash
#$1=fasta file

cp "$1" ../data
```
To execute it
```
bash fake_script.sh file1.fasta
```

This is the script get_sequences_from_list_of_loci.sh
```
#!/bin/bash
#retrieve sequences from a list of loci 
#$1=list of loci to grep in a fasta file
#$2=fasta file


file="$1" ; name=$(cat $file) ; for locus in $name ; do grep -w -A1 $locus "$2" ; done > "$1".fasta
```
Execute the script
```
bash get_sequences_from_list_of_loci.sh list_file M_musculus.fasta
```

## Conda
```
conda init bash # initialize conda
conda env list # see list of environments
conda activate "environment_name" # activate environment_name
conda list # see packages installed in that environment_name
conda deactivate # close environment 
```

## Useful Path:
/var/local/ # where are databases in our server

## Tips

> * Work in a robust and reproducible way
> * Document each step
> * Check everything between computational steps, errors can be silent
> * Code should be readable and organized in a logical way
> * Files, file names and folders organized in a logical way
> * Humans doing rote activities tend to make many mistakes, have your computer do as much of this rote work
as possible