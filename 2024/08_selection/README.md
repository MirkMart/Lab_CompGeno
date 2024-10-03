# Evolutionary selection

There are many programs that can dive into the various types of selection. The most straightforward way to study them is using the dN/dS ratio, also known as ω ratio.

From PAML suite documentation:

>The ω ratio is a measure of natural selection acting on the protein. Simplistically, values of ω < 1, = 1, and > 1 means negative purifying selection, neutral evolution, and positive selection. However, the ratio averaged over all sites and all lineages is rarely >1 since positive selection is unlikely to affect all sites over a prolonged time. Thus interest has been focused on detecting positive selection that affects only some lineages or some sites.

The two most common program to infer and study it are CodeML and [Hyphy](http://hyphy.org/#). The first is part of the [PAML suite](http://abacus.gene.ucl.ac.uk/software/paml.html), it is more raw and uses a config file. Once the obstacle of this config file is overcome, it is really easy to run. HyPhy is more user-friendly and presents different programs to test different things. The downside of this collection of models, is that they are perfectly tuned to search only for what they are describing, and it is difficult to adapt them to something that resamble to theri description. Both the programs are very time-consuming, so it is not advisable to run them on each sequence.

## CodeML

CodeML will take as inputs a sequence and a tree. Its result contains the omega ratios of intereset and the likelihood of the testes model. This likelihood then can be compared with a likelihood ratio test with the one of another CodeML run with a different model to identify which run best describe our data.

How CodeML works is extensively described in its [manual](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf).

There are three type of model that can be used:

- **branch models**: allow  the ω ratio to vary among branches in thephylogeny and are useful for detecting positive selection acting on particular lineages.
  
  - free-ratios model (1): assumes an independent ω ratio for each branch. It is parameter rich and discouraged.
  - several ω ratios (2): ave to specify how many ratios and which branches should have which rates in the tree file by using branch labels

- **site models**: allow the ω ratio to vary among sites.
- **branch-site models**: allow ω to vary both among sites in the protein and across branches on the tree and aim to detect positive selection affecting a few sites along particular lineages

### Pre-requisites

We will use the scrip [codeml.sh](./scripts/codeml.sh), which perform both CodeML and LRT analysis—this last one using the second script [LRT.py](./scripts/LRT.py). To use the script, it is important to have a precise folder structure. The script must be saved in the working directoy, which must contain three other folders:

- seq : contains the sequences of the orthogroups of interest (aligned and trimmed).
- tree1 : with as many sequences as those in seq/. These trees must be named after their respecitve sequence, with the extension '.nwk' instead of '.fna'. These trees serve to inform the program about how many omega classes we want to infer and specify how to group branches. This first fod the script will work with these tree identifying the full model.
- tree2: as above (tree1) but with the second model, the reduced one or nested.

In this case, sequences used are not constituted by amino acids, but nucleotides. We must perform what is known as retrotranslation.

Tree file format must follow the guidelines that are reported in the PAML documentation. The most ocmmon notation is using parethesis. Branches that are considered together during the omega inference must be labeled with the same tag, a number from 0 to the number of omega values that we want to compute specified with the symbol '#' (it is not needed to specify 0, it is the default one for not labelled branches). An example is given above (another is [example_terminal.nwk](./example_terminal.nwk)). Branch length can be present and are specified with the symbol ':'. In codeml, these lenghts can be used as starting points for the maximum likelihood inference. Names of the tips of the trees MUST be the same of those in the sequences.

>((Hsa_Human, Hla_gibbon) #1, ((Cgu/Can_colobus, Pne_langur), Mmu_rhesus), (Ssc_squirrelM, Cja_marmoset));

### Run the analysis

```bash
for seq in seq/*.fna; do bash codeml.sh "$seq" codeml.ctl; done
```
