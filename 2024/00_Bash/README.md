# Bash

## Server didattica (login using Guacamole)

[https://rlab.unibo.it](https://rlab.unibo.it)

login using your unibo email address and password

* To copy text (with Windows): select text to copy and paste it on command line using right click
* To copy text (with Mac): select text and press control+Shift+C to copy, control+Shift+V to paste
* If these methods do not work: Ctrl+Shift+Alt and use the text space provided as a bridge between Guacamole and any exterior windows

## Common ways to login to a Server

```bash
ssh username@ip # from Ubunt or Mac terminal
```

Using [putty](https://www.putty.org/) for Windows

## Paths, files and directories

Open linux or mac terminal; if you have windows download [git bash](https://gitforwindows.org/)

```bash
ssh username@ip #login to a server via ssh, it will ask for your password
Ctrl+d #close session
Ctrl+c #clean the command line
top # display Linux processes
whoami # username
who # which users are logged in
df -h # disk free
du -h # disk usage
pwd #print working directory
cd #change directory (relative or absolute path)
cd .. # bring you back one directory up
cd ~ #cd home
ls #print list of elements in directory
ls –lh # print sizes in human readable format
man list # (or any command, gives the manual for that command)
history #dispay last commands
ls *.fasta #list all fasta files
```

> **special characters:**

* ~ home directory
* . Current Directory
* / Path Directory Separator
* Space
* \* any character
* ; Shell Command Separator
* \> redirect output to new file
* \t tab
* \n end of line in bash
* \r end of line in mac
* \n\r end of line in windows
* many others (& | < ? []{} etc..)

---

## Permissions

```bash
sudo command #allows users to run programs if he/she has sudoer priviledges
chmod #change permissions of your files (or directory with chmod –r)  
```

![chmod](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/chmod.png)

---

## Edit files

`vi nomefile` #create new empty file

> **avoid special characters**; once created new file press “i" to write, after editing press Ctrl+c or esc and type `:wq` to save and exit from file or ':q!' to exit without saving

```bash
cp filename pathwheretocopy #copy file somewhere using absolute or relative path of where to copy
mv filename pathwheretocopy #mv file somewhere using absolute or relative path of where to copy
mv filename new_filename #rename file
less filename #see file
cat filename #similar to less
head filename #see first 10 rows of the file
tail filename #see last 10 rows of the file
head –n5 filename #(or tail –n 5) see only first 5 (or last) 5 rows
wc filename #outputs the number of words, lines, and characters
wc –l filename #outputs the number of lines
rm file #remove file
rm * #remove every file in that directory
mkdir newfoldername# make new directory
mkdir -p zmays-snps/{data/seqs,scripts,analysis} #create directory and subdirectories with one command
cp –r foldername #copy folder
mv foldername pathwheretomove #move folder
rm –r foldername #remove folder
```

> tip: be **VERY careful with `rm`**, once you removed something there is no way to undo it; remember bash is case sensitive, the file, folder or scritp "Data" is different from "data".
---

## Download and transfer data

`wget` can handle HTTP and FTP links

```bash
wget https://github.com/jacopoM28/CompOmics_2022/archive/refs/heads/main.zip
```

For public data we don’t need any authentication to gain access to this file, otherwise use flags *--user=*  and *--ask-password*

```bash
curl –O link #-O saves the file with its original filename
```

`curl` can transfer files using more protocols than wget. It supports FTP, FTPS, HTTP, HTTPS, SCP, SFTP, TFTP, TELNET, DICT, LDAP, LDAPS, FILE, POP3, IMAP, SMTP, RTMP and RTSP

`scp` (secure copy): transfer data from local computer to remote host, or from two remote hosts. scp works just like cp, except we need to specify both host and path

## Common way to transfer data

using Ubuntu or Mac terminal

```bash
scp username@ip:path/to/file/to/copy /where/to/paste/it # copy a file from remote host to local host
scp path/to/file/to/copy username@ip:/where/to/paste/it # copy a file from local host to remote host
scp -r username@ip:path/to/file/to/copy /where/to/paste/it # copy a directory
```

using Windows:
[Cyberduck](https://cyberduck.io/)

## Data integrity

```bash
md5sum nomefile
```

---

## Compress and decompress data

```bash
gzip file #compress file file.gz
bzip2 file #slower than gzip but higher compression ratio
gzip –k file #keep also the not compressed file
gunzip file.gz #uncompress file
zless file.gz #less compressed file
zgrep "word" file.gz #use grep in compressed file
```

With **gzip** you don't get compression across files, each file is compressed independent of the others in the archive, advantage: you can access any of the files contained within.

With **tar** the only gain you can expect using tar alone would be by avoiding the space wasted by the file system as most of them allocate space at some granularity.

In **tar.gz** compression: create an archieve and extra step that compresses the entire archive, you cannot access single files without decompressing them.

```bash
tar -zcvf myfolder.tar.gz myfolder # c create archive; z gzip archive; f specify new output; v verbose
tar xvfz ./nome_archivio.tgz #decompress archive
```

---

## Merge and sort files

```bash
cat file1 file2 file3 … #merge multiple files in 1  
cat file1 file2 file3 > newfilename #redirect output to new file
sort file #sort the file, careful to computational sorting of file
sort –h file #human numeric sort
sort -t file #specify field separator (e.g., -t",")
sort -k1,1 -k2,2n file #sort by first column adn then numerically by second column
sort -k1,1 -k2,2nr file #sort by first column adn then numerically by second column in reversed order
sort -k1,1V -k2,2n file #as before but human sorted
join -1 1 -2 1 sorted_file1 sorted_file2 #join two files according to first column of each file
join -1 1 -2 1 -a 1 sorted_file1 sorted_file2 #keep also non joined rows
paste file1 file2 #merge lines of files
```

---

## Compare two sorted files

```bash
diff -y file1 file2 #Compare FILES line by line and show side by side
comm file1 file2 #compare two sorted files line by line
comm -1 file1 file2 #lines unique to file1
comm -2 file1 file2 #lines unique to file2
comm -12 file1 file2 #print only lines present in both file1 and file2
comm -3 file1 file2 #print lines in file1 not in file2, and vice versa
```

---

## Grep

```bash
grep "word" file #print all rows that contains "word"
grep -w "word" file #print all rows that contains exactly the pattern "word"
grep -V "word" file #inverted match, print all rows that not contain the patter "word"
grep -i "word" file #ignore case distinctions, grep both "word" and "WORD"
grep -c "word" file #count how many rows contain the patter "word"
grep –A10 "word" file # print rows containing pattern "word" and the 10 rows after
grep –B10 "word" file # print rows containing pattern "word" and the 10 rows before
grep –C10 "word" file # print rows containing pattern "word" and the 10 rows after and before
grep "Locus10[12]" file #print Locus101 and Locus102 
greo -E "(Locus101|Locus102)" file #print Locus101 and Locus102 
```

special characters:

* \^ starting with ; grep "^>" file #print lines starting with ">"
* \$ ending with ; grep ">\$" file #print lines ending with ">"
* Regular Expressions: a sequence of characters that specifies a pattern

---

## Awk

```bash
awk '{print $1}' file #print first column
awk '{print $0}' file #print all columns
awk '{print $NF}' file #print last column
awk '{print $4"\t"$1}' file #change orders of column and use tab as field separator in the output
awk -F";" '{print $3,$4}' file #fiels separator is ";"
awk '$1==$2 {print}' file #if first column = second column, print all columns
awk '$1 ~ /chr2|chr3/ { print $0 "\t" $3 - $2 }' file #if first column contain "chr2" or "chr3", print all columns, and a column with the difference between $3 and $2
awk '$1 ~ /chr1/ && $3 - $2 > 10 '{print}' file #if both the first column  contain "chr1" AND $3-$2>0 , print all columns
awk '{if ($1~">") print $0; else print length}' "fasta_file #print length instead of sequence in fasta file
```

![awk](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/awk.png)

---

## Sed

```bash
sed 's/Locus/Transcript/' file #for each line subtitute "Locus" with "Transcripts" at first occurrance
sed 's/Locus/Transcript/g' file #for each line subtitute "Locus" with "Transcripts" at first occurrance
sed -i 's/Locus/Transcript/g' file # overwrite input with the output
sed '/Locus/d' file #delete any row containing "Locus"
```

> $ echo "chr1:28427874-28425431" | \
sed -E 's/^(chr[^:]+):([0-9]+)-([0-9]+)/\1\t\2\t\3/'
> chr1 28427874 28425431

The first component of this regular expression is `^\(chr[^:]+\):`. This matches the text that begins at the start of the line (the anchor ^ enforces this), and then captures everything between \( and \). The pattern used for capturing begins with “chr” and matches one or more characters that are not “:”, or delimiter. We match until the first “:” through a character class defined by everything that’s not “:”, [^:]+. The second and third components of this regular expression are the same: match and capture more than one number. Finally, our replacement is these three captured groups, interspersed with tabs, \t.

> $ echo "chr1:abRsZtjf-dhThdbUdj" | \
sed -E 's/^(chr[^:]+):([a-zA-Z]+)-([a-zA-Z]+)/\1\t\2\t\3/'
> chr1 abRsZtjf dhThdbUdj

---

## Conda

```bash
conda init bash # initialize conda
conda env list # see list of environments
conda activate "environment_name" # activate environment_name
conda list # see packages installed in that environment_name
conda deactivate # close environment 
```

---

## Git and GitHub

```bash
git add <file> #add file to the stagin area
git commit -m "commit message" #creates a new a commit 
git push #push any new commit
git pull #pull any changes from online repositories
git rm <file> #remove any file automatically adding this change to the staging area
git mv <file> #mv any file automatically adding this change to the staging area
git status #display the current state of the git repository
git rm --cached <FILENAME> #remove a file from those that are synchronised
```

---

## Bash scripting

All these commands can be combined to create small (or not) customed program that can help us resolving daily tasks in bioinformatics. More is present in [Bash scripting](./Bash_scripting.md)
