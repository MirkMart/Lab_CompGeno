# Computer Environment

### Conda environments :

Most of the software that you are going to use are already installed in different CONDA environments. You can have a look at [THIS](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) cheet - shet for some usefull tips on Conda. 

Environments and related softwares :

  1. **base**  :
  * *bamtools* 
  * *bcftools*
  * *blast*
  * *bowtie2*
  * *bwa*
  * *cd-hit*
  * *diamond*
  * *fastq-dumb*
  * *fastqc*
  * *hmmer*
  * *samtools*
  * *transdecoder*
  * *trimmomatic*
  * *trinity*
  * *wtdbg2*
  * *minimap2*
  * *canu*
  * *assembly_stats*

<br/>
  
  2. **Assembly_tools**
  * *ragtag*
  * *seqkit*
  * *quast*
  * *mummer*
  * *minimap2*
  * *busco*
  
<br/>

  3. **blobtools**
  4. **MAKER**
  5. **Hypo**

<br/>  

  6. **test_env**
  * *Orthofinder*
  * *mafft*
  * *iqtree*
  * *AMAS**

<br/> 
  
  7. **kakscalculator**
 
 <br/> 
 
  8. **Mosdepth**
  * *Mosdepth*
  * *CAFE5*

<br/> 

  9. **GAAS**
  

Some usefull conda commands:

```
conda create --name <ENV_NAME> #Create a target conda environment
conda activate <ENV_NAME> #Activate a target environment
conda deactivate #Deactivate your current environment
conda info --envs #print a list of conda environments
conda list #print a lost of packages installed in your current environment
conda install -c <CONDA_CHANNEL> <PACKAGE_NAME> #install a conda package
```

A good practice is to create AND install necessary packages at the same moment :

```
conda create --name <ENV_NAME> -c <CONDA_CHANNEL> <PACKAGE_NAME>
```

Conda is a package manager based on python. Each conda environment can only have **ONE** specific version of python installed. Now the default version is the latest 3.X but some old software can be run only in python 2.7. To create a python 2.7 environment :

```
conda create --name <ENV_NAME> python=2.7
```

---

### Usefull paths

```
/var/local/ # databases
/home/PERSONALE/jacopo.martelossi2/Data/Reads_Ex #File for fasta/fastq exercises
/home/PERSONALE/jacopo.martelossi2/Data/Anopheles_stephensi #Data for genome assembly and annotation
```

---

### Good practices in bioinformatics :

  1. **Work in a robust and reproducible way**
  2. **Document each step**
  3. **Check everything between computational steps, errors can be silent**
  4. **Code should be readable and organized in a logical way**
  5. **Files, file names and folders organized in a logical way**
  6. **Humans doing rote activities tend to make many mistakes, have your computer do as much of this rote work as possible**
  7. **Internet is your best friend and mentor, google everything that you don't understand!**
  8. **If an error rise first of all try to solve the problem by yourself:** a) read the error message carefully; b) read again the help of the software; c) check for typos, they are everywhere; d) GOOGLE it!

e.g: One easy and fast way to remember your code :

```
echo 'blastp -in myfile.fa -out myfile.blastp -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -max_hsps 1' > Blastp.sh
sh Blastp.sh
```

### Some usefull online resources :

  1. **Stack Overflow**
  2. **BioStars**
  3. **Issue page of GitHub**. Remember to remove ```is:open``` in filters bar.
  
---
  
# Study system

### *Anopheles stephensi* 
Arthropoda; Insecta; Pterygota; Diptera

<br/>

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/Figures/Aste.png)

<br/><br/>

**[LINK](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA629843/) to NCBI genome sequencing related bioproject**

---

# Task list

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/Figures/Tasks.png)



### Disclaimer
This course introduces various software tools for solving common tasks in comparative genomics projects. However, these are just examples, biased by our habits. In real - case scenario and especially at the beginning of your career can be a good practice to explore as better as possible the universe of alternative possibilities. Indeed, biology is **extremely** complex and one software could be the best for a particular dataset but perform extremely worst for another one. In bioinformitcs we are lucky enough to not pay for each analyses, usually. So if you have the possibility try to be curious and test different tools and parameters combination.  

