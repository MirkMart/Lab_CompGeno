# Genome assembly2

## Contaminant detection

to perform contaminats detection we will use [Blobtools](https://blobtools.readme.io/docs) (this the [original publication](https://f1000research.com/articles/6-1287)).

![Blobtools](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/5_GenomeAssembly2/Figures/Blobtools.png)

Blobtools is a bioinformatics tool that assesses and analyses the quality and composition of genome assemblies. It visualizes and classifies contigs based on taxonomic, GC content, and coverage information, making it easier to detect contaminants or misassembled sequences in a genome assembly.

### Preliminary steps

1. Re-mapping short-read data to the assembly

    ```bash
    minimap2 --secondary=no --MD -ax sr -t <NUMBER_CORE> <ASSEMBLY> <FASTQ_R1> <FASTQ_R2> | samtools view -Sb - > <OUT_BAMFILE>.bam
    samtools sort -@<n_cores> -o <OUT_SORTED_BAMFILE>_sorted.bam Anoste_pol.bam
    rm Anoste_pol.bam
    samtools index Anoste_pol_sorted.bam
    ```

2. Taxonomic annotation of the contigs

    ```bash
    blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>
    ```

### Construction a BlobDB data structure and assembly visualisation (~6 minutes)

Install it by following the [GitHub instructions](https://github.com/DRL/blobtools). **REMEMBER** to create nodes.db. At the end, add the folder to `$PATH` to run the program from anywhere. The output will be created by appending `.blobDB.json`.

```bash
blobtools create -i <ASSEMBLY> -b <MAPPING_FILE> -t <BLASTN_FILE> -o <OUTPUT_PREFIX>
blobtools view -i <JSON_FILE> -o <OUTPUT_PREFIX>
blobtools plot -i <JSON_FILE> -o <OUTPUT_PREFIX>
```

You can download and inspect files in `.png` format directly. A common way to observe the table file, instead, is:

```bash
grep -v '^##' Anoste_blobDB_table.txt | column -t -s $'\t' | less`
```

Useful script to remove a list of sequences (works also with multiline FASTA):

```bash
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' in.fa | grep -w -v -Ff patterns.txt - | tr "\t" "\n" > out.fa
#The first part ensures the FASTA is in one-line form and reformats the file, obtaining">header\tsequence". This makes it easier to process the file with grep. Then it filters out everything in the patterns file, outputting a one-line FASTA file with only the desired contigs.
```

With:

- -w #ensures to match the exact name
- -v #invert the recognised pattern
- -Ff #compares directly to file
- '-' #use stdin as input from the previous command

To use it, we need to create a file listing all contaminated contigs. To do so, we have to decide which parameters we can use to filter our blobDB_table.

-----

## Reference-based scaffolding

For this purpose, we will use [RagTag](https://github.com/malonge/RagTag). This collection of software tools scaffolds and improves modern genome assemblies. Tasks include: Homology-based misassembly correction, Homology-based assembly scaffolding and patching, and Scaffold merging

1. Reference-based error correction

    ```bash
    ragtag.py correct -t 20 <REFERENCE_GENOME> <DRAFT_GENOME>
    ```

2. Scaffolding**

    ```bash
    ragtag.py scaffold -C -t 20 -o <OUTPUT_DIR> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>
    ```

    Where : -C = concatenate unplaced contigs and make ‘chr0’

reference genome: `/home/PERSONALE/mirko.martini3/2024/Data/Anopheles_reference/GCF_013141755.1_UCI_ANSTEP_V1.0_chromosomes.fasta`
