Genome assembly2
================

# Contaminants detection

[Blobtools](https://blobtools.readme.io/docs)

![](https://raw.githubusercontent.com/jacopoM28/CompOmics_Tutorship/main/2023/5_GenomeAssembly2/Figures/Blobtools.png)

### Preliminary steps :

**1. Re - mapping short reads data to the
    assembly**

    minimap2 --secondary=no --MD -ax sr -t <NUMBER_CORE> <ASSEMBLY> <FASTQ_R1> <FASTQ_R2> | samtools view -Sb - > <OUT_BAMFILE>.bam
    samtools sort -@10 -o <OUT_SORTED_BAMFILE>.sorted.bam Aste.corrected.renamed-sr.bam
    rm Aste.corrected.renamed-sr.bam
    samtools index Aste.corrected.renamed-sr.sorted.bam

**2. Taxonomic annotation of the
    contigs**

    blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>

### Construction a BlobDB data structure and assembly visualitation

    ./blobtools create -i <ASSEMBLY> -b <MAPPING_FILE> -t <BLASTN_FILE> -o <OUTPUT_PREFIX>
    ./blobtools view -i <JSON_FILE> -o <OUTPUT_PREFIX>
    ./blobtools plot -i <JSON_FILE> -o <OUTPUT_PREFIX>

<br/> Usefull script to remove a list of sequences (works also with
multiline fasta)
    :

    awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' in.fa | grep -w -v -Ff patterns.txt - | tr "\t" "\n" > out.fa

-----

# Reference - based scaffolding

For this purpose we will use [RagTag](https://github.com/malonge/RagTag)

**1. Reference based error correction**

    ragtag.py correct -t 20 <REFERENCE_GENOME> A<DRAFT_GENOME>

**1.
    Scaffolding**

    ragtag.py scaffold -C -t 20 -o <OUTPUT_DIR> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>

Where : -C = concatenate unplaced contigs and make
‘chr0’
