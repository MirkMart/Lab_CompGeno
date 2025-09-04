# Bash scripting

## Background Processes

Programms can take a while and, in the meantime, we do not have access to the command line until they have finished or we kill the program (Ctrl+c). To avoid that we work in backgroung. We can you different tool, but the two most common are `screen` and `tmux`. Personally, I prefer the second beacuse it manages with more agility mutiple terminal sessions, allowing to switch between them easily.

with **tmux**. There are many shortcuts for tmux and an helpful [cheat sheet](https://tmuxcheatsheet.com/):

```bash
tmux new -s <session_name> #open a new session names
tmux a -t <session_name> #attach to a screen named
tmux ls #list all existing sessions
Ctrl+a+d #close screen
Ctrl+d #kill screen
Ctrl+c #create a new window in the current session
Ctrl+n #move to the window on the right (next) in the current session
Ctrl+p #move to the window on the left (previous) in the current session
Ctrl+| #split vertically the current window
Ctrl+_ #split horizontally the current window
Ctrl+[ #enter in read mode
```

with **screen**:

```bash
screen #open new screen
screen -R <screen_name> #open a new screen name screen_name
Ctrl+a+d #close screen
Ctrl+d #kill screen
screen -ls #list of all screen attached
screen -r <screenidnumber> #open screen with that idnumber 
```

---

## Concatenate commands and programs

### Pipe |

`|` connects the standard output of one process to the standard input of another

```bash
grep "word" file1 | sed 's/ /\t/g' | program1 > file2 #sed use as input the output of grep; when using pipe the "original" input has to be specified only in the first command
```

### semicolon ';'

`;` concatenate different commands or programs sequentially

```bash
grep ">" file1.fasta >output1 ; grep "_" file2.fasta output2
```

### &&

Concatenate two programs so that program2 run only if program1 completed successfully

```bash
program1 input.txt > intermediate-results.txt && program2 intermediate-results.txt > results.txt
```

---

## Standard output and standard error

Many programs use a "standard output" for outputting data, we usually redirected the standard output in an output file using ">". A separate file, called "standard error" is needed for errors, warnings, and messages. We can redirect the standard error with "2>"

```bash
program1 file 2> program1.stderr > results.txt
```

---

## Variables

A variable in bash can be anything, a number, a character, a string of characters, a file, a folder. Variables are in bash are indicated with $

```bash
var1="ciao"
echo $var1 #print ciao
echo "$var1" #print ciao
echo '$var1' #print var1
```

---

## Command and process substitution

**Command substitution** runs a Unix command inline and returns the output as a string that can be used in another command. It is used when we need a **text** as the input of the following command.

```bash
echo "$(cat file.fasta)"
echo "There are $(grep -c '^>' input.fasta) entries in my FASTA file." # show the string "There are 416 entries in my FASTA file."
```

**Process substitution** runs a Unix command and return a FIFO pipeline, a file that does not store data permanently. It is used when we need a **file** as the input of the following command.

```bash
join <(sort -k1,1 input1.txt) <(sort -k1,1 input2.txt) > output.txt
```

Process substitution can be used also in **output** while we want to sort in different way something produced by a command.

```bash
bash script.sh >(grep "warning" > warninigs.log) >(grep "error" > errors.log) > script.log
```

---

## For loop

```bash
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

Bash scripts are indicated with the '.sh' extention (python scripts with '.py', perl scripts with '.pl').

Create bash script with nano:

```bash
#!/bin/bash
mkdir test
for i in *fasta; do mv $i test; done
cd test
for i in *fasta; do program1 $i > $i"_output"; done
for i in *output; do grep ">" $i; done > list_of_sequences
```

We need now to make the .sh file executable. Normally it is already executable. This could be needed if some program uses this script and encounters some problems while trying to run it.

```bash
chmod 777 namescript.sh
```

To execute the bash script

```bash
bash namescript.sh
```

## Bash script with variables

```bash
nano fake_script.sh
```

```bash
#!/bin/bash
#$1=fasta file

cp "$1" ../data
```

To execute it

```bash
bash fake_script.sh file1.fasta
```

We can also name our variable so they will be more clear when called.

```bash
#!/bin/bash

list_OG=$1
aln_dir=$2
nucleo_dir=$3

mkdir -p "$aln_dir"/old

for OG in $(cat "$list_OG"); do
    cp "$aln_dir"/"$OG"_aligned.faa "$aln_dir"/old/"$OG"_aligned_old.faa
    sed -i -f <(awk '{print "/" $0 "/{N;d;}" }' <(grep -v -w -Ff <(grep ">" "$nucleo_dir"/"$OG".fna) <(grep ">" "$aln_dir"/"$OG"_aligned.faa))) "$aln_dir"/"$OG"_aligned.faa
done
```
