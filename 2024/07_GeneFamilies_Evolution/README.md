# GENE FAMILIES EVOLUTIONARY ANALYSES

## INTRO

Evolutionary studies about gene families are becoming very frequent in recent years. For example, in almost all genome paper is present a gene families expansion/contractions analyses. Indeed, has been observed how enormous amount of changes in families size can take place even between closely related organisms. There is much interest in these changes, as even the gain or loss of single genes have been implicated in adaptive divergence between species. In addition, large contractions and expansions of gene families are generally attributed to natural selection.

One of the most popular software is **[CAFE](https://academic.oup.com/bioinformatics/article/22/10/1269/237347)**, even if other apporaches based on gene and species tree exist. Its purpose is to analyze changes in gene family size in a way that accounts for phylogenetic history and provides a statistical foundation for evolutionary inferences. The program uses a birth and death process to model gene gain and loss across a user-specified phylogenetic tree. The distribution of family sizes generated under this model can provide a basis for assessing the significance of the observed family size differences among taxa. In breaf, CAFE try to estimates ancestral states of gene families comparing the tip states (one continouse state for each gene family).

---

## PREPARING THE DATA

Preparing the data for CAFE it's a really tedious process. Fortunatly Orthofinder already provide us the Orthogroup table which is almost in the correct format. Just be sure that species name match between the Orthogroup table and the species tree and that it's present an additional first column with random values. In my case the formatted file is:

```text
NONE  Orthogroup  Aaeg  Aste  Cqui  Llon
NONE  OG0000000 93  0 5 9
NONE  OG0000001 81  0 5 0
NONE  OG0000002 33  21  25  1
NONE  OG0000003 10  5 8 0
NONE  OG0000004 0 0 5 90
NONE  OG0000005 0 0 58  1
NONE  OG0000006 26  0 16  9
NONE  OG0000007 1 0 4 91
NONE  OG0000008 0 0 0 48 
```

Tip to add a NONE column:

```bash
# Adding an extension after -i is usefull to create a backup file automatically.
sed -i.old $'s/^/NONE\t/g' <INFILE> 
```

Tip to remove last column of a file:

```bash
rev <INFILE> | cut -d$'\t' -f 2- | rev > <OG GENE COUNT CAFE READY> 
```

After converting your nexus file into newick we are ready to use CAFE. We can ru the program inferring one or more lambdas. Before any attempt, it is common to estimated a general error model:

```bash
cafe5 -i <OG GENE COUNT> -t <NEWICK TIME TREE> -o Error_model -e
```

With:

- -i #input gene count tsv
- -t #time tree in newick format
- -o #the name of the output folder. It is greatly recommended to do so
- -e #will infer the error model

Then we can run the program specifying the error model just inferred.

```bash
cafe5 -i <OG GENE COUNT> -t <NEWICK TIME TREE> -o CAFE_1Lambda -eBase_error_model.txt
```

with:

- -e<name_error_model> #important not to leave any space. As above, but here you specify the error model just inferred.

To compute mutiple lambdas you must tell CAFE how many different λs there are, and which species or clades share these different λs. The lambdas and their locations are specified in a tree file (that will specified with -y. It is a different tree file). An example tree file:

```text
(Llon:1,(Aste:1,(Cqui:2,Aaeg:2):1):1);
```

To run the 2 lambda analyses taking into consideration the error model:

```bash
cafe5 -i Orthogroups.GeneCount_CAFE.tsv -t TimeTree_CAFE.nwk -o 2L -y TimeTree_CAFE_2l.nwk -eError_model/Base_error_model.txt
```

For a detailed description of CAFE outputs, see the [manual](https://github.com/hahnlab/CAFE5). `*` is used to identify significant changes. Try to use grep to find significantly changed trees or gene family (orthogroups). These one would be great candidates to further functional studies.

Once you have your set of gene of intereset you can perform [enrichment analyses](../09_GeneAnnotation_functional_enrichment/) using TopGO.
