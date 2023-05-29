# GENE FAMILIES EVOLUTIONARY ANALYSES

## INTRO: 

Evolutionary studies about gene families are becoming very frequent in recent years. For example, in almost all genome paper is present a gene families expansion/contractions analyses. Indeed, has been observed how enormous amount of changes in families size can take place even between closely related organisms. There is much interest in these changes, as even the gain or loss of single genes have been implicated in adaptive divergence between species. In addition, large contractions and expansions of gene families are generally attributed to natural selection.

One of the most popular software is **[CAFE](https://academic.oup.com/bioinformatics/article/22/10/1269/237347)**, even if other apporaches based on gene and species tree exist. Its purpose is to analyze changes in gene family size in a way that accounts for phylogenetic history and provides a statistical foundation for evolutionary inferences. The program uses a birth and death process to model gene gain and loss across a user-specified phylogenetic tree. The distribution of family sizes generated under this model can provide a basis for assessing the significance of the observed family size differences among taxa. In breaf, CAFE try to estimates ancestral states of gene families comparing the tip states (one continouse state for each gene family).

---

<br/>

## PREPARING THE DATA

Preparing the data for CAFE it's a really tedious process. Fortunatly Orthofinder already provide us the Orthogroup table which is almost in the correct format. We need just a little bit of bash scripting to change few things:

```
sed -i.old 's/_2//g' <OG GENE COUNT>
sed -i.old 's/_//g'  <OG GENE COUNT>
cat  <OG GENE COUNT> | awk '{print tolower($0)}' > <OG GENE COUNT 2> 
sed -i.old $'s/^/NONE\t/g' <OG GENE COUNT 2> 
rev <OG GENE COUNT 2>  | cut -d$'\t' -f 2- | rev > <OG GENE COUNT CAFE READY> 
```

Now we need to change a bit also our time tree. First of all we need to convert the nexus file into newick. For this operation I usally use figtree. Just open the nexus treefile - ensuring to have divergence time annotated - and save it as newick. Then :

```
cat <TIME TREE> | awk '{print tolower($0)}' > <TIME TREE>
```

Now we are ready to use CAFE

```cafe5 -i <OG GENE COUNT CAFE READY> -t <TIME TREE> -o CAFE_1Lambda -p
```
