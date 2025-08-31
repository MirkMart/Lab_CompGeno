# Kmer - based genome survey

## What is a *k-mer*?

It’s a metematical concept representing substrings (*mers*) of a certain length (*k*) contained within a biological sequence. *K-mers* are not related to real, physical molecules.

![k-mer](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/Figures/kmers.png)

### Kmer based genome survey

1. **Genome size estimation**
2. **Repeat Content estimation**
3. **Heterozigosity estimation**

For a given sequence of length L, and a k-mer size of k, the total k-mer’s possible (*n*) will be given by :

```math
( L – k ) + 1
```

Since we sequenced a genome multiple time the total length (N) will be :

```math
n = [( L – k ) + 1] * C
```

```math
N = n/C
```

In real case scenarios we have a **non-uniform** coverage of the genome ( **Repetitive content**; **Sequencing errors**; **heterozigosity**; ecc…)

-----

## Preliminary step: QC of short reads data

Necessary softwares :

1. **[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**
2. **[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**

### Fastqc (20 mins cad)

It is a quality control tool for high throughput sequence data.

```bash
fastqc seqfile1 seqfile2 .. seqfileN
```

### Trim reads in a sliding window approach: **TRIMMOMATIC** (20 mins)

```bash
trimmomatic PE -threads 20 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```

where :

- 'PE' : specifies that trimmomatic should run in paired-end mode
- '-threads 20' : use 20 threads for parallel processing
- '-phred33' : type of quality score format
- the next six are input and output names that MUST be given in a specific order:
  - first input (1)
  - second imput (2)
  - first output (first input) for trimmed and paired forward reads reads (1)
  - second output (first input) for the trimmed and unpaired forward reads (those that cannot be paired) (1)
  - first output (second input) for trimmed and paired forward reads reads (2)
  - second output (second input) for the trimmed and unpaired forward reads (those that cannot be paired) (2)
- 'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10'
  - path to adapters fasta;
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

See the [manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for deeper explanation of parameters.

Other possibilities :

- **bbduck** from [BBMap](https://sourceforge.net/projects/bbmap/) suite
- **[TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)**

-----

## Compute *k-mer* frequency (~ 7 mins for 20 threads)

Analysing k-mers can be done with many different tools ([KAT](https://kat.readthedocs.io/en/latest/index.html), [Jellyfish](https://github.com/gmarcais/Jellyfish), or [BBmap](https://sourceforge.net/projects/bbmap/)). We will use KAT.

KAT (K-mer Analysis Toolkit) is a bioinformatics toolset used to perform various analyses on k-mer distributions from sequence data, primarily in genome assembly and evaluation. K-mers are used to analyze the structure, quality, and completeness of genome assemblies, and to compare datasets such as sequencing reads and assembled genomes.

```bash
kat hist -t <THREADS> -m <MER LENGTH> -o <OUTPUT PREFIX> <paired_1> <paired_2>
```

In particular:

- '-t' : inform the program about how many thread (parallelisation) it should use (4).
- '-m' : is the dimension of the k-mer we will use (27).
- '-o' : determines the prefix of the output.
- finally, append the two input files.

## Genome size, heterozigosity and repetitive content estimation

The file we will use is the one without any extension. We can rename it using the suffix '.hist'. To use it, we have to delete all the comments in  its beginning, then we can upload it into [genomescope](http://genomescope.org/genomescope2.0/).

A useful description of how [interpret a genomescope graph](https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/genomescope.html#gsc.tab=0)

[This](http://genomescope.org/genomescope2.0/analysis.php?code=vWJpnAj5PQRPhLKl3lTZ) is how our result should like.

**Watch out for lower and upper count boundaries**.

-----

### Some other example

[*Ruditapes philippinarum*](./Data/Rphil_kmer27.png)
[*Reticulitermes lucifugus*](./Data/Rluc.kmc_30_Genomescope.png)
