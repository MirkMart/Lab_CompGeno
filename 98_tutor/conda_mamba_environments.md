# Environments in Didattica's Server

This file containes all conda/mamba envirments installed:

## base (last update 02/07/2025)

- R             4.5.1

## tree (last update 04/06/2025)

- gotree        0.4.5
- iqtree        2.3.6
- raxml-ng      1.2.2
- paml          4.10.7
- ete3          3.1.3
- treeswift     1.1.45
- HYPHY         2.5.71
- orthofinder   2.5.5
- mafft         7.526

## sequence (last update 02/07/2025)

Python too high for orthofinder

- blobtools     1.1.1
- bmge          1.12 (outside all conda environment exist also bmge 2.0)
- BUSCO         6.0.0
- edirect       24.0
- mafft         7.526
- ncbi-datasets 18.3.1
- SRA-tools     3.2.1

[where download busco datasets](https://busco-data.ezlab.org/v5/data/lineages/)

[busco datasets on server](/usr/local/share/busco_database):

- arthropoda: /usr/local/share/busco_databases/arthropoda_odb12
- culicidae: /usr/local/share/busco_databases/culicidae_odb12

## assembly (last update 01/09/2025)

python 3.13.5

- augustus      3.5.0
- blast         2.14.1
- diamond       2.1.10
- fastqc        0.12.1
- hypo          1.0.3
- maker         3.01.04
- minimap2      2.28
- mosdepth      0.3.10
- r-base        4.3.3
- ragtag        2.1.0
- samtools      1.21
- trimmomatic   0.36

## kat (last update 29/08/2025)

- kat           2.4.2

```bash
tar -xvzf *.tar.gz
```

## GAAS (last update 25/09/2025)

- agat          1.4.1
- GAAS          1.2.0
