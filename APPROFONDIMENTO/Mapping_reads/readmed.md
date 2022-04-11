# Mapping reads

## [BWA](http://bio-bwa.sourceforge.net/)

Software package for mapping low-divergent sequences against a large reference. Three different algorithms:
+ **BWA-backtrack** designed for Illumina sequence reads up to 100bp
+ **BWA-SW**  for longer sequences (from 70bp to 1Mbp)
+ **BWA-MEM**  recommended, faster and more accurated; for longer sequences (from 70bp to 1Mbp)

1. Index the referece

    ``` bwa index ref.fa```

2. Mapping reads against the indexed reference

    ```bwa mem ref.fa read1.fq read2.fq > aln-pe.sam```

## [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

Ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s 

1. Index the referece

    ``` bowtie2-build <reference_in> <bt2_base>```
    
    bt2_base: The basename of the index files to write

2. Mapping reads against the indexed reference

    ```bowtie2 -x <bt2-idx> -1 <mate1> -2 <mate2> -S <sam> > alignment_rate```
    
    --no-discordant By default, bowtie2 looks for discordant alignments if it cannot find any concordant alignments. This option disables that behavior.

    ... take  a look at the options in the manual
---
# SAM and BAM files

SAM: Sequence Alignmetn Map. SAM files are TAB-delimited text files, they are the format for mapping data. In SAM files, ach sequencing read has an alignment entry --> alignment data files are massive and require space-efficient complex binary file formats

The binary analog of SAM is BAM, which stores the same data in a compressed binary representation (if BAM file = 1M, SAM file = 4G).

SAM format consist of a header section and an alignment section.

[SAMTOOLS](http://www.htslib.org/doc/samtools.html): tool to extract the information in sam and bam files and explore data

```
samtools view -H samfile #look at the header
samtools view samfile # look at samfile without header
samtools view -h samfile # look at the whole samfile
samtools view -b input.sam > output.bam #from sam to bam
samtools view -h input.bam > output.sam  #from bam to sam
```

**Header:**

+ @SQ information about the reference sequences
+ @RG important read group and sample metadata
+ @PG info about the programs used to create and processa the SAM/BAM file

Example @RG:     ID:1    LB:lib1 PL:Illumina     SM:Sample_01

**Alignment section:**

Contains read alignments. Each alignment entry is composed of 11 required fields (and optional fields after this).

1. **QNAME** query name (read name)

2. **FLAG** the bitwise flag, which contains information about the alignment

3. **RNAME** reference name (chromosome name, transcript name..)

4.	**POS** the position on the reference sequence

5. **MAPQ** mapping quality
6. **CIGAR** the CIGAR string
7. **RNEXT**  reference name of a paired-end read’s partner ("*" info unavailable; "=" NNEXT is the same of QNAME; if != discondart mapping between forward and reverse read) 
8. **PNEXT** position of a paired-end read’s partner
9. **TLEN** (Template LENgth) distance between the mapped end of the template and the
mapped start of the template
10. **SEQ** the original read sequence, in the orientation it is aligned
11. **QUAL** the original read base quality 


## FLAGS
![Flags samtools](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/flag_samtools.png)

```
samtools flags 147 # 0x93 147 PAIRED,PROPER_PAIR,REVERSE,READ2
samtools flags 0x93` # 0x93 147 PAIRED,PROPER_PAIR,REVERSE,READ2
samtools flags paired,proper_pair # 0x3 3 PAIRED,PROPER_PAIR
```

## CIGAR
![Cigar](https://raw.githubusercontent.com/MariangelaIannello/didattica/main/images/cigar.png)


## Mapping Quality

Q = -10 log10(P) 

P = probabliy of incorrect mapping position

mapping quality=20 --> 1% of probability that the alignment is incorrect



## Filter SAM files and index BAM files

```
samtools view #see options
samtools -f INT # keep reads with the flags in INT
samtools -F INT # exclude reads with the flasgs in INT
samtools view -h -f 0x2 -F 256 -q 30 -Sb input.sam > file.bam # keep only proper paired mapping reads (-f 0x2), with minimum mapping quality of 30 (-q 30), exclude secondary alignments (-F 256, or -F 0x200) and save the output in bam format (-Sb)
samtools sort file.bam > file_sorted.bam # can be very computationally intensive
samtools index file_sorted.bam # index the file so that we can access to specific regions
samtools view file_sorted.bam M_assembly49:10000-12000 #show reads that align to the "M_assembly49" chromosome, from position 10000 to 12000

samtools idxstats input.sortedbam #stats about mapping reads for each chromosome (prints reference name, reference length, number of mapped reads and number of unmapped reads)
```
## View alignment

(Bam file must be sorted and indexed)

```samtools tview -p assembly3:1-500 Sample_23_merged.sorted.bam MF.fasta```

Use arrows to move, Crtrl+h to move 1000b forward, ctrl+l 1000b backward


> . = base that matched the reference on the forward strand 

> , = base that matched the reference on the reverse strand

> AGTCN = base that did not match the reference on the forward strand

> agtcn = base that did not match the reference on the reverse strand

Alternative: [igv](http://software.broadinstitute.org/software/igv/home)
