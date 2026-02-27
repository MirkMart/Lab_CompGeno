# Computer Environment

## Conda environments

Most of the software that you are going to use are already installed in different CONDA environments. You can have a look at [THIS](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) cheat sheet for some usefull tips on Conda.

Environments and related softwares:

## base

- R             4.5.1

## tree

- cafe          5.1.0
- disco         1.4.1
- ete3          3.1.3
- gotree        0.4.5
- HYPHY         2.5.71
- iqtree        3.0.1
- mafft         7.526
- paml          4.10.7
- orthofinder   2.5.5
- raxml-ng      1.2.2
- treeswift     1.1.45

## sequence

Python too high for orthofinder

- agat          1.4.1
- blobtools     1.1.1
- bmge          1.12 (outside all conda environment exist also bmge 2.0)
- BUSCO         6.0.0
- edirect       24.0
- mafft         7.526
- ncbi-datasets 18.3.1
- SRA-tools     3.2.1

[busco datasets on server](/usr/local/share/busco_database):

- arthropoda: /usr/local/share/busco_databases/arthropoda_odb12
- culicidae: /usr/local/share/busco_databases/culicidae_odb12

Use `$BUSCO/name_db` to use a specific db.

## assembly

python 3.13.5

- assembly-stats 1.0.1
- augustus       3.1
- blast          2.16.0
- diamond        2.1.10
- fastqc         0.12.1
- hypo           1.0.3
- maker          3.01.04
- minimap2       2.28
- mosdepth       0.3.10
- multiqc        1.31
- r-base         4.3.3
- samtools       1.21
- spades         4.2.0
- trimmomatic    0.40

## kat

- kat           2.4.2

---

### Usefull paths

```bash
/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads #Reads for genome assembly
/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/01_Anoste_reference #Reference file of Anopheles stephensi from NCBI
```

---

### Good practices in bioinformatics

  1. **Work in a robust and reproducible way**
  2. **Document each step**
  3. **Check everything between computational steps, errors can be silent**
  4. **Code should be readable and organized in a logical way**
  5. **Files, file names and folders organized in a logical way**
  6. **Humans doing rote activities tend to make many mistakes, have your computer do as much of this rote work as possible**
  7. **Internet is your best friend and mentor, google everything that you don't understand!**
  8. **If an error rise first of all try to solve the problem by yourself:** a) read the error message carefully; b) read again the help of the software; c) check for typos, they are everywhere; d) **GOOGLE it**!

e.g: One easy and fast way to remember your code :

```bash
echo 'blastp -in myfile.fa -out myfile.blastp -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -max_hsps 1' > Blastp.sh
sh Blastp.sh
```

### Some usefull online resources

  1. [Stack Overflow](https://stackoverflow.com/)
  2. [BioStars](https://www.biostars.org/)
  3. **Issue page of GitHub**. Remember to remove ```is:open``` in filters bar.

---

## Git and GitHub

One extremely useful resource we now have is Git along with GitHub. More is provided in the dedicated file [00_Git_GitHub](./00_Git_GitHub.md).
