# Species collection and dataset creation

To start a project in comparative genomics, you need genomes to compare. Once species of intereset have been identified (no more than 5-6 + _Anopheles stephensi_) you need to download their data.

## Download using NCBI datasets

Modern machines can support specific tools. One example is `dataset` which has been designed to work with NCBI datasets, the new NCBI repository for all the genomic, taxonomic, and genetic content of the famous center.

Using NCBI geome datasets, filter species of interest (those already annotated), choose and select the ones you want to use, and list their accession number in a file. Reading attentively the written information, NCBI kindly provides us the command we should run. It is similar to the following one, which has simply been cut to retrieve only the two files we need to continue with our tutorial.

```bash
#|sequence|
datasets download genome accession <ACCESSION_NUMBER> --include gff3,genome
```

Since species we want are at least 4, I recommend you to use a for loop to iterate through many accession number, making the process easier and more direct. Otherwise, to ease and make faster the process, you can use the script [download_dataset](../../00_practice/99_scripts/download_dataset.sh).

You can find the script here `/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/99_scripts/download_dataset.sh`

> Pay attention to the formattation of the input file required by the script.

## Direct download

A second more "direct" approach is using the ftp protocol and the command `wget`.

From the genome dataset page of a spcies, it is possible to reach the FTP site of the species itself, where you can double check if annotation is present (a `*.faa` file should be present, but more in particular you are interested in `.gff`). Manually speaking, from this FTP page you can copy the link of a file you need and use `wget` to retrieve it.

## Longest isoform and translation

In a genome, after annotation, it is not uncommon to have multiple isoforms of the same sequence (in particular if you used a transcriptome to increase the precision of your annotation). To not bias your analyses, you have to remove these isoforms using the program [AGAT](https://github.com/NBISweden/AGAT). AGAT works with several different [pearl script](https://agat.readthedocs.io/en/latest/index.html). To keep longest isoforms you will use [agat_sp_keep_longest_isoform.pl](https://agat.readthedocs.io/en/latest/tools/agat_sp_keep_longest_isoform.html). Then, you will extract, and translate, sequences annotated as genes using [agat_sp_extract_sequences.pl](https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html).

```bash
agat_sp_keep_longest_isoform.pl --gff <GFF_file> -o <output_file>
agat_sp_extract_sequences.pl -g <GFF_longest_file> -f <FASTA_file> -t cds -p --cfs --output <output_file>
#N.B. AGAT use a particular module that wants FASTA file to be wrappend. Here it is important NOT to have single line FASTA. If you already have, try to use the command 'fold'
```

option description:

- --gff/-g #gff file name
- --output/-o #output file name
- -f #fasta file name
- -t #type of sequence to extract. In this case cds (coding sequences)
- -p #enable translation from nucleotide to amino acid
- --cfs #clean_final_stop deletes all * symbols that indicate a STOP codon at the end of the translated sequence (which is expected)

**Remember**: files are many. Use a for cycle (or [a scipt](../../00_practice/99_scripts/AGAT_longest_extract.sh)...)

### Pseudogene deletion

Even though you cannot untrust everything is made by someone else or by a program, it is true that sometimes _trust is good, but control is better_. Annotation files can contains mistakes and sequences labelled with the CDS tag can contain STOP codons. In this case, you are looking at a pseudogene, not a gene, and you must delete it from your proteome/genome since you are interested only in working coding sequences.

I designed the script [pseudogene_find_eliminate](../../00_practice/99_scripts/pseudogene_find_eliminate.sh) to do such work for you, deleting every sequence that contains a `*`. Remember to check the passages that are performed and, moreover, if something is needed before run it.

```bash
bash pseudogene_find_eliminate.sh
```

You can find the script here `/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/99_scripts/pseudogene_find_eliminate.sh`

## Days of Future Past

Since you are already playing with gff file, run also the command

```bash
agat_sp_extract_sequences.pl -g {<GFF_longest_file>} -f {FASTA_file} -t cds --cfs -roo --output {output_file}
```

In this case, you can see that it is absent the option `-p`, so you will not translate nucleotide sequences into amino acid, but it is present `-roo` (--remove_orf_offset). This flag is indispensable because without it all proteins that present their first codon in a frame different from the first (encoded by AGAT with 0) will be extracted erroneously, with heading nucleotides that will change the frame subsequently elaborated. With the `-p` flag, AGAT automatically takes into consideration the protein phase, but this is not true when `-t CDS` is performed alone.
