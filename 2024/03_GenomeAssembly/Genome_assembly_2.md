# Genome assembly2

## Contaminants detection

to perform contaminats detection we will use [Blobtools](https://blobtools.readme.io/docs) (this the [original publication](https://f1000research.com/articles/6-1287)).

![Blobtools](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/5_GenomeAssembly2/Figures/Blobtools.png)

Bloobtools is a bioinformatics tool designed to assess and analyse the quality and composition of genome assemblies.It allows for the visualization and classification of contigs based on taxonomic, GC content, and coverage information, making it easier to detect contaminants or misassembled sequences in a genome assembly.

### Preliminary steps

1. Re - mapping short reads data to the assembly

    ```bash
    minimap2 --secondary=no --MD -ax sr -t <NUMBER_CORE> <ASSEMBLY> <FASTQ_R1> <FASTQ_R2> | samtools view -Sb - > <OUT_BAMFILE>.bam
    samtools sort -@10 -o <OUT_SORTED_BAMFILE>.sorted.bam Aste.corrected.renamed-sr.bam
    rm Aste.corrected.renamed-sr.bam
    samtools index Aste.corrected.renamed-sr.sorted.bam
    ```

2. Taxonomic annotation of the contigs

    ```bash
    blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>
    ```

### Construction a BlobDB data structure and assembly visualitation (~6 minutes)

Installation perfomed following [GitHub instructions](https://github.com/DRL/blobtools). **REMEMBER** yo create nodes.db. At the end, the folder was added to the `$PATH` to call the program from wherever. The output will be created appending `.blobDB.json`.

```bash
blobtools create -i <ASSEMBLY> -b <MAPPING_FILE> -t <BLASTN_FILE> -o <OUTPUT_PREFIX>
blobtools view -i <JSON_FILE> -o <OUTPUT_PREFIX>
blobtools plot -i <JSON_FILE> -o <OUTPUT_PREFIX>
```

File in `.png` format can be dowloaded and directly inspectionated. A common way to observe the table file, instead, is:

```bash
grep -v '^##' Anoste_blobDB_table.txt | column -t -s $'\t' | less`
```

Usefull script to remove a list of sequences (works also with multiline fasta):

```bash
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' in.fa | grep -w -v -Ff patterns.txt - | tr "\t" "\n" > out.fa
#The first part ensures the FASTA is in oneline form. Then Filter out everything is contained inside the patterns file, outputting a fasta file oneline with only desired contigs.
```

With:

- -w #ensures to mach the exact name
- -v #invert the recognised pattern
- -Ff #compares directly to file
- - #use as stdin the stdout of the precedent command

To use it we have to create a file where we list all contaminated contigs. To do so, we have to decide based on which parameters we can filter our blobDB_table.

-----

## Reference-based scaffolding

For this purpose we will use [RagTag](https://github.com/malonge/RagTag). This is a collection of software tools for scaffolding and improving modern genome assemblies. Tasks include: Homology-based misassembly correction, Homology-based assembly scaffolding and patching, and Scaffold merging

1. Reference based error correction

    ```bash
    ragtag.py correct -t 20 <REFERENCE_GENOME> <DRAFT_GENOME>
    ```

2. Scaffolding**

    ```bash
    ragtag.py scaffold -C -t 20 -o <OUTPUT_DIR> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>
    ```

    Where : -C = concatenate unplaced contigs and make ‘chr0’

reference genome: `/home/PERSONALE/mirko.martini3/2024/Data/Anopheles_reference/GCF_013141755.1_UCI_ANSTEP_V1.0_chromosomes.fasta`
