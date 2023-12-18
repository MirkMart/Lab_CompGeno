Genome assembly
================

# Some theory

## What is a genome assembly ?

[Ensembl
definition](https://www.ensembl.org/info/genome/genebuild/assembly.html)
“A genome assembly is a computational representation of a genome
sequence. Because we are not able to sequence along the complete length
of a chromosome, each chromosome assembly is made up of short stretches
of sequenced DNA pasted together.”

Basically, assemblers needs to solve giant jigsaw puzzles with millions
of pieces and in which some pecies are missing altogether and some
pieces contain errors. **They perform some of the most complex
computations in all biology**

-----

## Main steps in basic genome assembly :

  - 1.  Contig assembly + error correction (depending on the type of
        reads)

  - 2.  Scaffolding

  - 3.  Gap
filling

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/Genome_assembly.png)

#### Some defitinitions :

**Contig** : “a contiguous sequence generated from determining the
non-redundant path along an order set of component sequences. A contig
should contain no gaps but often the terms contig and scaffold are used
interchangeably”

**Scaffold** : “an ordered and oriented set of contigs. A scaffold will
contain gaps, but there is typically some evidence to support the contig
order, orientation and gap size estimates.”

**Assembly** : “a set of chromosomes, unlocalized and unplaced (random)
sequences and alternate loci used to represent an organism’s genome.
Most current assemblies are a haploid representation of an organism’s
genome, although some loci may be represented more than once (see
Alternate locus, above). This representation may be obtained from a
single individual (e.g. chimp or mouse) or multiple individuals
(e.g. human reference assembly). Except in the case of organisms which
have been bred to homozygosity, the haploid assembly does not typically
represent a single haplotype, but rather a mixture of haplotypes. As
sequencing technology evolves, it is anticipated that diploid sequences
representing an individual’s genome will become available.”

**Diploid assembly** : “A genome assembly for which a Chromosome
Assembly is available for both sets of an individual’s chromosomes, as
defined by the NCBI Assembly model. It is anticipated that a diploid
genome assembly is representing the genome of an individual. Therefore
it is not anticipated that alternate loci will be defined for this
assembly, although it is possible that unlocalized or unplaced sequences
could be part of the assembly.”

**Primary Assembly** : “Relevant for haploid assemblies only. The
primary assemblies represents the collection of assembled chromosomes,
unlocalized and unplaced sequences that, when combined, should represent
a non-redundant haploid genome. This excludes any of the alternate locus
groups.”

**Unlocalized sequence/scaffold** : “A sequence found in an assembly
that is associated with a specific chromosome but cannot be ordered or
oriented on that chromosome.”

**Unplaced Sequence/scaffold** : “A sequence found in an assembly that
is not associated with any chromosome.”

Source : [Genome Reference
Consortium](https://www.ncbi.nlm.nih.gov/grc/help/definitions/#:~:text=The%20primary%20assemblies%20represents%20the,a%20non%2Dredundant%20haploid%20genome.)

### 1\. Contig - level assembly

Main assembly algorithms :

  - 1.  **Overlap Layout Consensus
(OLC)**

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/Overlap-Layout-Consenus.png)

  - 2.  ***De Bruijn Graph*
approach**

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/DBG.png)

### 2\. Scaffolding

  - 1.  **De - novo scaffolding**
    
    <!-- end list -->
    
      - Genetic mapping
      - Bacterial Artificial Chromosome–End Sequencing
      - Linked-Read Sequencing
      - Optical Maps
      - Proximity Ligation : Regions of the genome that are close
        together in sequence generally have more frequent physical
        contact than parts of the genome that are far apart in sequence.
        Nevertheless, regions of the same chromosome, even those
        megabases away, contact each other more often than they contact
        other chromosomes.

  - 2.  **Reference - based/Sinteny - based
scaffolding**

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/Scaffolding.png)

### 3\. Gap - filling

-----

# What can influence our genome assembly?

  - 1.  Nature and characteristic of the genome being sequenced (Ploidy,
        repetitive content, heterozigosity ecc…).

  - 2.  Quality of the extracted DNA.

  - 3.  Type of sequence technology.

  - 4.  Reads coverage.

  - 5.  Assembly software.

**The heterozigosity problem
:**

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/Heterozigosity-1.png)

Tipically, more moneys you have, more high-quality data you get and more
likely you will assemble a good quality genome. However moneys are not
all and nowdays big genome projects rely on big constortia with
bioinformaticians extremely experienced in genome assembly and able manually curated them (*e.g* [Darwin Tree of
Life](https://www.darwintreeoflife.org/data/), [Vertebrate Genomes
Project](https://vertebrategenomesproject.org/)).

-----

# How can we evaluate our assembly?

Some, but not all, metrics :

1.  **N50 :** the length of the shortest contig for which longer and
    equal length contigs cover at least 50 % of the assembly. It gives
    an evaluation of the **contiguity** of an assembly.
2.  **[Busco](https://busco.ezlab.org/) :** a set of Benchmarked
    Universal Single-Copy Orthologs (BUSCO) that sould be present along
    a certain taxonomic level. Busco genes are marked as C:complete
    \[D:duplicated\], F:fragmented, M:missing. It gives an evaluation
    about the completness and the redundancy of our assembly.
    Precomputed single copy OG can be found
    [here](https://busco-data.ezlab.org/v5/data/lineages/)
3.  **spectra-cn
    ([KAT](https://kat.readthedocs.io/en/latest/walkthrough.html)) :**
    “the assembly spectra copy number plot checks assembly coherence
    against the content within reads that were used to produce the
    assembly. Basically we represent how many elements of each frequency
    on the read’s spectrum ended up not included in the assembly,
    included once, included twice
etc.”

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/KAT.png)

-----

# Practice

## Contig level assembly :

We will use **[Wtdbg2](https://github.com/ruanjue/wtdbg2)** a very fast
assemblers, but keep in mind that many others exist and each one of them could perform better with different kind of data(*e.g*
[CANU](https://github.com/marbl/canu);
[Falcon](https://github.com/falconry/falcon);[HiFiasm](https://github.com/chhylp123/hifiasm))

#### Assemble reads, costruction of contig layout and edge sequences

    wtdbg2 -x rs -g <EXPECTED_GENOMESIZE> -t <NUMBERS_OF_CORES> -i <FASTQ> -fo <OUT_PREFIX> #Default parameters

#### Produce final consensus in fasta format

    wtpoa-cns -t <NUMBER_OF_CORES> -i <LAYOUT_FILE> -fo <OUT_PREFIX>

See the GitHub page for usefull tips (Remeber to check also the “issue”
page)

## Genome polishing with short and long reads

Softwares : [Minimap2](https://github.com/lh3/minimap2);
[Hypo](https://github.com/kensung-lab/hypo);
[Samtools](http://www.htslib.org/)

**Pre - requisite :** Mapping short and long reads
    data

    minimap2 -ax <PRESET_Options> --MD [other_options] <target.fa> <fastq>
    samtools view -Sb <SAM FILE> > <BAM FILE>
    rm <SAM FILE>
    samtools sort -@<NUM_THREADS> -o <OUT_FILE> <INFILE>
    samtools index <INFILE>
    echo -e "$R1\n$R2" > <READS_PATH>

> OPTIONAL Calulcate mean coverage of short and long reads with
> [Mosdepth](https://github.com/brentp/mosdepth)

**Hypo**

    hypo -d <DRAFT_Contigs> -r <READS_PATH> -s <APPROXIMATE_GENOMESIZE> -c <SHORT_READSCOVERAGE>

## Genome evaluation

#### N50

    assembly-stats <ASSEMBLY>

#### Busco

    busco -m <MODE> -l <LINEAGE> -c <CPU_NUMBER> -o <OUTPUT_NAME> -i <INPUT>

#### KAT

    kat comp -t <NUM_THREADS> -o <OUTPUT_PREFIX> '<FASTQ>...' <ASSEMBLY>

Compare your results with :  
- *[Ruditapes
philippinarum](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/R.philippinarum_KmerSpectra.png)*
- *[Reticulitermes
lucifugus](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/4_GenomeAssembly/Figures/R_lucifugus-main.mx.spectra-cn.png)*

> EXTRA : To remove redundant contigs (*i.e* haplotigs)
> *[Purge\_dups](https://github.com/dfguan/purge_dups)* is one of the
> most used softwares. It uses both long - reads coverage and self-alignments 
> informations to identify and remove poorly collapsed regions
> of the assembly.
