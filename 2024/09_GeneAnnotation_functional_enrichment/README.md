# Annotation

Annotation tools that find similarities between sequences:

* [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [Diamond](https://github.com/bbuchfink/diamond)

Annotation tool using probabilistic models:

* [HMMER](http://hmmer.org/): use hidden Markov models (HMMs) profiles. It is often used together with a profile database, such as Pfam or many of the databases that participate in Interpro.

## Input file

The first thing to functionally explore our data is creating a file that can be annotate. The best one is a list of protein choosen by selecting the longest sequence in each trimmed orthogroup without taking into consideration gaps. This secure the use of the longest protein that is for sure clustered in an orthogroup and implemented in our analyses. To do it we can use the script [longest_protein_OGs](./00_script/longest_protein_OG.sh).

This is the path to the script. Remeber to add the realpath of your Orthogroup_Sequences folder (one of the Orthofinder output) to complete the script.

`/home/PERSONALE/mirko.martini3/00_Lab_CompGeno/2024/08_GeneAnnotation/00_script/longest_protein_OG.sh`

## Databases

* Nr: Non-redundant protein sequences from GenPept, Swissprot, PIR, PDF, PDB, and NCBI RefSeq
* Nt: nucleotide sequences
* Swiss-Prot: manually annotated and reviewd proteins ([UniProt](https://www.uniprot.org/))
* [Pfam](http://pfam.xfam.org/): is a large collection of protein families, each represented by multiple sequence alignments and hidden Markov models
* See also [InterPro consortium](http://www.ebi.ac.uk/interpro/)

## Annotate sequences with Diamond

[Diamond](https://github.com/bbuchfink/diamond) is optimized for large input files of >1 million proteins.

The program may use quite a lot of memory and also temporary disk space. Should the program fail due to running out of either one, you need to set a lower value for the block size parameter -b.

**-b**: Block size in billions of sequence letters to be processed at a time. This is the main parameter for controlling the program’s memory and disk space usage. Bigger numbers will increase the use of memory and temporary disk space, but also improve performance. The program can be expected to use roughly six times this number of memory (in GB).

### Makedb

```bash
diamond makedb --in /var/local/diamond_db/nr.gz --db ./nr_diamond
```

where:

* `--in` file: Path to the input protein reference database ﬁle in FASTA format (may be gzip compressed).
* `--db`: Path to the output DIAMOND database ﬁle.

### Diamond search

The sensitivity can be adjusted using the options `--mid-sensitive`, `--sensitive`, `--more-sensitive`, `--very-sensitive` and `--ultra-sensitive`. Even if more time consuming, the last one is always the recommended one.

(see [Blast manual](https://www.ncbi.nlm.nih.gov/books/NBK279668/#usermanual.BLAST_search_strategies) and [Diamond manual](https://github.com/bbuchfink/diamond))

**Bit-score**: the requires size of a sequence database in which the current match could be found just by chance. The higher the bit-score, the better the sequence similarity. The bit-score gives the same value for hits in databases of different sizes. It is independent of query sequence length and database. The bit-score depends on the raw alignment score. Thus, the higher the bit score, the more highly significant the match is.

**Blast e-value**: number of expected hits of similar quality (score) that could be found just by chance, given the same size of a random database. The E-value (expectation value) is a corrected bit-score adjusted to the sequence database size.

```math
E = (m x n) / (2^bit-score)
```

where:

* m = query sequence length
* n = total database length (sum of all sequences)

**pident**:  % of identical matches

## HMMER Search

* **hmmbuild** build profile from input multiple alignment
* **hmmsearch** search profile against sequence database
* **hmmscan** search sequence against profile database

```bash
hmmscan [-options] /var/local/Pfam/Pfam-A.hmm <seqfile>
```

* **-o** \<f> Direct the main human-readable output to a file `f` instead of the default stdout.
* **--tblout** \<f> Save a simple tabular (space-delimited) file summarizing the per-target output, with one data line per homologous target model found.
* **--domtblout** \<f> Save a simple tabular (space-delimited) file summarizing the per-domain output, with one data line per homologous domain detected in a query sequence for each homologous model.
* **--pfamtblout** \<f> Save an especially succinct tabular (space-delimited) file summarizing the per-target output, with one data line per homologous target model found.
* **-E** \<x> In the per-target output, report target profiles with an Evalue of <= \<x>. The default is 10.0
* **-T** \<x> Instead of thresholding per-profile output on E-value, instead report target profiles with a bit score of >= \<x>.
* **--domE** \<x> In the per-domain output, for target profiles that have already satisfied the per-profile reporting threshold, report individual domains with a conditional E-value of <= \<x>. The default is 10.0.
* **--domT** \<x> Instead of thresholding per-domain output on E-value, instead report domains with a bit score of >= \<x>.

## Gene Ontology

Useful links:

* [geneontology](http://geneontology.org/)
* [QuikGO](https://www.ebi.ac.uk/QuickGO/)
* [amigo](http://amigo.geneontology.org/amigo)

### Get GO terms from protein sequences

GOterms annotation can be performes with a graet variety of programs. [PANNZER](http://ekhidna2.biocenter.helsinki.fi/sanspanz/) and [eggNOG mapper](http://eggnog-mapper.embl.de/) are two common online tool. They are often reduntant, so we will annotate our proteome using the command line program [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/index.html). unfortunately, this cannot be done on the server we are currently using, so as soon as you have your longst_sequence file ready give it to the tutor to perform the annotation.

```bash
/home/PERSONALE/dbs/interproscan-5.65-97.0/interproscan.sh -i <LONGEST_PROTEINS_INPUT> -goterms -pa -b <OUTPUT-FILE-BASE> -cpu <N_CPUS>
```

With:

* -goterms #enables the GO annotation
* -pa #enables the REACTOME pathway annotation
* -b #output file name

### GO enrichment

To see if our genes of interest (for example genes differentially expressed in one condition), show an enrichment in some functional annotation.

We use [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html). topGO requires:

1. **Gene universe** the complete list of genes
2. **Genes of interest**
3. **GO annotation** with the following format ([example](./go2genes.txt)):

    ```text
    gene_ID\tGO_ID1, GO_ID2, GO_ID3, ....
    ```

There are two types of test statistics we can choose: Fisher's exact test which is based on **gene counts**, and a Kolmogorov-Smirnov like test which computes enrichment based on **gene scores**, which consider p-values. There are also many algorithms that can be choosen to perform a ranking of the GO terms obtained. We will use `elim` because it can decrease the redundancy of the annotation result taking into consideration GO term hierarchy.

An example of how to perform an enrichment analysis is present [here](./GO_enrichment.R).

One tools to visualize GO enrichment: [REVIGO](http://revigo.irb.hr/)

## KEGG PATHWAY

KEGG (Kyoto Encyclopedia of Genes and Genomes) is a comprehensive database and bioinformatics resource used to understand high-level functions and utilities of biological systems. One of the greatest strength of KEGG is that almost each known gene or pathways has been associated to a KEGG identifier. K terms work no too dissimilar than GOterms, but gene ones identify the orthogroup of belonging. As GoTerms, also KEGGterm can be analysed with enrichment anaylsis using the R package `clusterProfiler`.

Kegg annotated genes can be used to group together genes tha belong to the same pathway, in order to focus on them if particalarly interesting. Again, annotation program are many, but the one we suggest you to utilise is [KAAS](https://www.genome.jp/kegg/kaas/).

Useful links are:

* [pathways db](https://www.genome.jp/kegg/)
* [KEGG orthology](https://www.genome.jp/kegg/ko.html)

Functional enrichments can be perfomred also on KEGG pathways and KEGG orthologues. If you are intrested I can provide you the script to it. It is not too dissimilar to GO enrichment, but change subtly in a couple of passages.
