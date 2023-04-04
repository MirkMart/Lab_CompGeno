# *Omics* data and related usefull bash commands

## Fastq

FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character for brevity.

FASTQ file: four lines per sequence. 
* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence


![Phred score](https://raw.githubusercontent.com/jacopoM28/CompOmics_2022/main/Figures/Screenshot_2022-03-04%20Phred%20quality%20score%20-%20Wikipedia.png)

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

## Annotation files (GFF and BED)

**GFF**

```
Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.  
    source - name of the program that generated this feature, or the data source (database or project name).  
    feature - feature type name, e.g. Gene, Variation, Similarity.  
    start - Start position* of the feature, with sequence numbering starting at 1.  
    end - End position* of the feature, with sequence numbering starting at 1.  
    score - A floating point value.  
    strand - defined as + (forward) or - (reverse).  
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..  
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.  
```
For more informations see [here](https://www.ensembl.org/info/website/upload/gff.html).  

**BED**

```
The first three required BED fields are:

    - chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).  
    - chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.  
    - chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position notation chr1:1-100.  

The 9 additional optional BED fields are:

    - name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.  
    - score - A score between 0 and 1000.  
    - strand - Defines the strand. Either "." (=no strand) or "+" or "-".  
```

For more informations see [here](http://genome.ucsc.edu/FAQ/FAQformat#format1).  

[Example](https://github.com/jacopoM28/CompOmics_Tutorship/blob/main/2023/1_FastaFastq/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff.gz) of GFF3 file. [Example](https://github.com/jacopoM28/CompOmics_Tutorship/blob/main/2023/1_FastaFastq/Cgig_Genes.bed.gz) of BED file.

[BEDTOOLS](https://bedtools.readthedocs.io/en/latest/) is a comprehensive toolset for working with BED/GFF3 files. If you want to perform arithmetic operations in a genome (e.g merge, compare, extend intervals)

## Exercises and usefull commands

Necessary files are stored in ```/home/PERSONALE/jacopo.martelossi2/Data/Reads_Ex```
#### IMPORTANT :
Do NOT copy the file, get used to [SYMLINKS](https://linuxize.com/post/how-to-create-symbolic-links-in-linux-using-the-ln-command/), don't waste space with garbage!

1. Count the total number of reads in ech file  

  Tips :
  
  * Files are compressed, use the appropriate command.
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
awk '{s+=$1}END{print s}' file'
```

4. Count number of sequences in the fasta file ```Example.fa```

5. Extract the sequence ```XP_001647772.2``` and store it in a new file. How was the protein annotated?

  Tips :
  
  * Watch out for "nested" patterns! 
