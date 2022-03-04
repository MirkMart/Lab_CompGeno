# *Omics* data and related usefull bash commands

## Fastq

FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character for brevity.

FASTQ file: four lines per sequence. 
* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence

![header](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/illumina_seq_id.png)

![ascii](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/ascii_2.png)

![ascii_2](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/ascii.png)

![ascii_3](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/ascii33.gif)

---

## Fasta

FASTA format is a text-based format for storing biological sequences.

FASTA file: usually two lines per sequence.
* Line 1 begins with a '>' character and is followed by a sequence identifier and an optional description
*	Line 2 is the sequence.

**NB** Fasta files are usually multi - line (*i.e.* after sequence name there are more line with a fixed length for each sequence), however in bash is more straightforward to work with oneliners. A usefull command to conver a fasta file from multi - line to oneliner :

```
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < infile.fa > outfile.fa
``` 

---

## Exercises and usefull commands

1. Count the total number of reads in ech file (N.B you canf find them in /home/PERSONALE/jacopo.martelossi2/Data/Reads_Ex)  

#### IMPORTANT :
Do NOT copy the file, get used to [SYMLINKS](https://linuxize.com/post/how-to-create-symbolic-links-in-linux-using-the-ln-command/), don't waste space with garbage!

  Tips :
  
  * Files are comrpessed, use the appropriate command.
  * Remember how is structured a fastq file (How many lines are presents for each sequence ?).
  * See [here](https://linuxhint.com/bash_arithmetic_operations/) for Bash Arithmetic Operations
  
2. Which is the mean read length per file?

  Tips :
  
  * Think about the sequencing strategy.
  * Usefull ```awk``` to calculate mean value of a column, ([here](https://stackoverflow.com/questions/19149731/use-awk-to-find-average-of-a-column) a full explanation) :
  
```
awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }'
```
**NB** Remember to change accordly the column number

3. Calculate mean expected coverage for each library given a genome size of 12Mb. What should you change if you had both PE reads?

  Tips :
  
  * Usefull ```awk``` to sum values of a column :
  
```
awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }'
```

4. Count number of sequences in the fasta file ```Example.fa```

5. Extract the sequence ```XP_001647772.23``` and store it in a new file. How was the protein annotated?

  Tips :
  
  * Watch out for "nested" patterns! 
