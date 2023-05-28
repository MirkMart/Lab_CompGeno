## Estimate divergence time in a Maximum Likelihood framework

Traditionally, divergence time estimation was performed through **Bayesian analyses** (*e.g* Beast) and this is still the common and prefered practice when analyzing from few to a moderant amount of loci. However nowadays, the large amount avaible data from NGS projects has made the use of software based on maximum likelihood more and more frequent, especially when the divergence time estimation is not the main goal of the research (in that case Beast is still more used). Bayesian analyses is indeed able to naturally take into account and display a lot of different sources of uncertainty, from model selection (actually it is not a model selection but a  **model averanging**), to topology and obviously divergence time. Moreover Bayesian inference has the great advantage (but also disadvantage) to be able to incorporate prior in the form of probability distributions. However, it is also a double edge sword, since the long computational time this type of analysis requires.

As previusly said ML can overcome this issue and, as almost everything, we can carry on this analyses again with IQ-TREE! It implements the least square dating (LSD2) method to build a time tree when you have date information for tips or ancestral nodes ([here](https://academic.oup.com/sysbio/article/65/1/82/2461506) the link for the original paper).

This are the mains functions provided by IQ-TREE:

![](https://github.com/for-giobbe/phy/blob/master/2021/Images/LSD.png)

**Tip dating** is commonly used in two different scenarios:

 * 1. When analyzing viruses
 * 2. When your have a tree builded up including extinct taxon.

**Ancestral dates** are used when you have informations from fossils and you want to constrain, usually between two boundaries, the age of a/multiple specific node/s.

In this tutorial we are going only towards this latest way. Moreover, because it's not a course focused on divergence time estimation and we only need an ultrametric tree for gene family evolutionary analyses, we will use previously estimted divergence time to calibrate our tree. In my case I relied on [this](https://resjournals.onlinelibrary.wiley.com/doi/abs/10.1111/syen.12489#:~:text=Molecular%20divergence%20time%20estimates%20revealed,Jurassic%20(approximately%20197.5%20Mya).) and on [time tree](https://timetree.org/) for the root age. 

The only file that we have to prepare is a date file in which we have to specify the calibration point. It should look something like:

```
taxon1,taxon2 -50
taxon3,taxon4,taxon5 -100
taxon6 -10
```

which, for example, mean that the most recent common ancestor (MRCA) of taxon1 and taxon2 was 50 mya (million year ago) and the MRCA of taxon3, taxon4, taxon5 was 100 mya. Note that **no empty space** should be added to the comma-separated list of taxa, as empty space is used as a separator between taxon list and dates.

Now we arer ready to perform our divergence time estimation

```
iqtree -s ../Aln/concatenated.out --date <CALIBRATION FILE> --date-tip 0 -o <OUTGROUP> -m TESTNEW -nt 6 --prefix Time.Tree --date-options "-u 1" 
```



