# K-mer-based genome survey

## What is a *k-mer*?

It’s a mathematical concept that represents substrings (*mers*) of a specific length (*k*) in a biological sequence. *K-mers* are not related to real, physical molecules.

![k-mer](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/Figures/kmers.png)

### K-mer-based genome survey

1. **Genome size estimation**
2. **Repeat Content estimation**
3. **Heterozigosity estimation**

For a given sequence of length L, and a k-mer size of k, the total number of possible k-mers (*n*) will be given by :

```math
( L – k ) + 1
```

Since we sequenced a genome multiple times, the total length (N) will be :

```math
n = [( L – k ) + 1] * C
```

```math
N = n/C
```

In real case scenarios we have a **non-uniform** coverage of the genome ( **Repetitive content**; **Sequencing errors**; **heterozigosity**; ecc…)

-----

## Preliminary step: QC of short-read data

Necessary softwares :

1. **[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**
2. **[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**

### FastQC (20 min each)

It is a quality control tool for high-throughput sequencing data.

```bash
fastqc seqfile1 seqfile2 .. seqfileN
# Sum all results into a single report
multiqc .
```

Reads are here:

> /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads/SRR11672503_1.fastq.gz
> /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads/SRR11672503_2.fastq.gz

#### Read quality results

[Analysis modules](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

- Per base sequence quality: overview of the range of quality values across all bases at each position in the FastQ file.
- Per Sequence Quality Scores: allows you to see if a subset of your sequences have universally low quality values. A subset of sequences may have universally poor quality, often because they are poorly imaged (e.g., at the edge of the field of view); however, these should represent only a small percentage of the total sequences.
- Per Base Sequence Content: proportion of each base position in a file for which each of the four normal DNA bases has been called
- Per Sequence GC Content: measures the GC content across the whole length of each sequence in a file and compares it to a modelled normal distribution of GC content
- Per Base N Content: percentage of base calls at each position for which an N was called
- Sequence Length Distribution: graph showing the distribution of fragment sizes in the file which was analysed (warning if all sequences are not the same length).
- Sequence Duplication Levels: counts the degree of duplication for every sequence in a library and creates a plot showing the relative number of sequences with different degrees of duplication
- Overrepresented Sequences: lists all of the sequences which make up more than 0.1% of the total
- Adapter Content: generic analysis of all of the Kmers in your library to find those which do not have even coverage through the length of your reads. This can find a number of different sources of bias in the library which can include the presence of read-through adapter sequences building up on t

##### RNA

- per-base sequence content (there's "random priming" performed that isn't completely random, so the first ~8 bases show a typical bias pattern).
- per-sequence GC content (particularly in sRNA-seq you're enriching for a subset of genes, they won't have some nice GC distribution and that's perfectly OK).
- Sequence duplication level (if this doesn't fail, then everything is uniformly expressed...which would be quite unlikely)
- Overrepresented sequences (as above, though this will be worse for sRNA-seq since reads will often cover the entire transcript, so there's even less sequence diversity)

### Trim reads in a sliding window approach: **TRIMMOMATIC** (20 mins)

```bash
trimmomatic PE -threads 20 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```

where :

- 'PE': specifies that Trimmomatic should run in paired-end mode
- '-threads 20': use 20 threads for parallel processing
- '-phred33': type of quality score format
- the next six are input and output names that MUST be given in a specific order:
  - first input (1)
  - second input (2)
  - first output (first input) for trimmed and paired forward reads (1)
  - second output (first input) for the trimmed and unpaired forward reads (those that cannot be paired) (1)
  - first output (second input) for trimmed and paired forward reads (2)
  - second output (second input) for the trimmed and unpaired forward reads (those that cannot be paired) (2)
- 'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10'
  - path to adapters FASTA;
  - seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed;
  - palindromeClipThreshold: specifies how accurate the match between the two ‘adapter ligated’ reads must be for PE palindrome read alignment;
  - simple clip threshold: the minimum score of the full alignment between adapter andcd ../ read for the clipping to take place
- 'LEADING:3' #Remove leading low quality or N bases (below quality 3)
- 'TRAILING:3' #Remove trailing low quality or N bases (below quality 3)
- 'SLIDINGWINDOW:4:15'
  - windowSize: specifies the number of bases to average across;
  - requiredQuality: specifies the average quality required.
- 'MINLEN:36'
- '2> stats_trimmomatic': standard error redirection to a file named stats_trimmomatic

See the [manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for a deeper explanation of parameters.

Other possibilities :

- **bbduck** from [BBMap](https://sourceforge.net/projects/bbmap/) suite
- **[TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)**

-----

## Compute *k-mer* frequency (~ 7 mins for 20 threads)

Analysing k-mers can be done with many different tools ([KAT](https://kat.readthedocs.io/en/latest/index.html), [Jellyfish](https://github.com/gmarcais/Jellyfish), or [BBmap](https://sourceforge.net/projects/bbmap/)). We will use KAT.

KAT (K-mer Analysis Toolkit) is a bioinformatics toolset for analyzing k-mer distributions from sequence data, primarily for genome assembly and evaluation. K-mers help assess the structure, quality, and completeness of genome assemblies and compare datasets such as sequencing reads and assembled genomes.

```bash
kat hist -t <THREADS> -m <MER LENGTH> -o <OUTPUT PREFIX> <paired_1> <paired_2>
```

In particular:

- '- t ': informs the program about how many threads (parallelisation) it should use (4).
- '-m': is the dimension of the k-mer we will use (27).
- '-o': determines the prefix of the output.
- Finally, append the two input files.

## Genome size, heterozygosity and repetitive content estimation

We will use the file without any extension. We can rename it using the suffix '.hist'. To use it, delete all the comments at the beginning, then upload it to [GenomeScope](http://genomescope.org/genomescope2.0/).

A useful description of how to [interpret a GenomeScope graph](https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/genomescope.html#gsc.tab=0)

[This](http://genomescope.org/genomescope2.0/analysis.php?code=vWJpnAj5PQRPhLKl3lTZ) is what our result should look like.

**Watch out for lower and upper count boundaries**.

-----

### Some other example

[*Ruditapes philippinarum*](./Data/Rphil_kmer27.png)
[*Reticulitermes lucifugus*](./Data/Rluc.kmc_30_Genomescope.png)
