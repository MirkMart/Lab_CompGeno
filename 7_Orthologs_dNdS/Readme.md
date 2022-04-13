# Orthologs detection and rate of protein evolution between two species

## Orthologs


> **Orthologs** are genes in different species that evolved from a common ancestral gene by speciation
> **Paralogs** are genes originated by duplication events in one species

![orthologs_paralogs](https://raw.githubusercontent.com/jacopoM28/CompOmics_2022/main/Figures/0.png)

[Orthofinder](https://github.com/davidemms/OrthoFinder)

![Best reciprocal hits](https://raw.githubusercontent.com/jacopoM28/CompOmics_2022/main/Figures/Comparison-of-inparalog-groups-Blast-best-hits-are-used-to-define-the-potential.png)

## Pairwise Sequence Alignment

![nt alignemnt](https://miro.medium.com/max/724/0*r_fPoz7mGydNBmbp.png)

We use **proteins** and not nucleotides to align ortholog sequences, since two orthologs can differ markedly in DNA sequence but still encode the same amino acid sequence.

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_2022/main/Figures/Screenshot_2022-04-07%20Genomics%20and%20Comparative%20Genomics2.png)
--> two organisms which are distantly related may have many differences in their DNA sequence, but have less differences in their amino acid sequence and thus can be compared at this level more easily.

## Retrotranslate aa alignments
To infer rate of protein evolution we need to use **nucleotide** alignments
We use the aa alignement as guide to align nucleotides. For doing this we need 1) the aa alignment obtained before 2) the non aligned nt sequences.

![nt adn aa sequence](https://raw.githubusercontent.com/jacopoM28/CompOmics_2022/main/Figures/Screenshot_2022-04-07%20Genomics%20and%20Comparative%20Genomics.png)

## Rate of protein evolution

dN/dS is the ratio of the number of nonsynonymous substitutions per non-synonymous site (pN) to the number of synonymous substitutions per synonymous site (pS), which can be used as an indicator of selective pressure acting on a protein coding gene. A ratio greater than one implies positive or darwinian selection, less than one implies purifying selection, and a ratio of one indicates neutral selection.
>n= total number of sites
m= number of substitutions
To calculate Ka and ks we need to count the number of synonymous (S) and nonsynonymous counts (N) (S+N=n) and the numbers of synonymous (Sd) and nonsynonymous (Nd) substitutions (Sd +Nd=m). After correcting for multiple substitutions we obtain Ka and ks (the observed n of substitutions underestimates the real number of substitutions happened in time).

Three steps to calculate Ka/Ks:
1)	Count S and N
2)	Count Sd and Nd
3)	Correct for multiple substitutions

> Models for calculating Ka and Ks adopting different substitution models 

## GO enrichment

Apply statistical tests to verify if genes of interest (for example genes with faster evolution) are more often associated to certain biological functions than what would be expected in a random set of genes. 

We use [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html). topGO requires:

1. **Gene universe** the complete list of genes
2. **Genes of interest**
3. **GO annotation** with the following format:

    gene_ID\<TAB>GO:ID1,GO:ID2,GO:ID3, ....

One tools to visualize GO enrichment: [REVIGO](http://revigo.irb.hr/)

---

## Workflow


![workflow](https://raw.githubusercontent.com/jacopoM28/CompOmics_2022/main/Figures/workflow_ortologhi.png)
