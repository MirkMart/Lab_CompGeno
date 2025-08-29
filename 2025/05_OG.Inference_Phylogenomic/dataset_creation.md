# Species collection and dataset creation

To start a project in comparative genomics, you need genomes to compare. Once species of intereset have been identified (no more than 5-6 + _Anopheles stephensi_) you need to download their data.

## Download using NCBI datasets

Modern machines can support specific tools. One example is `dataset` which has been desgined to work with NCBI datasets, the new NCBI repository for all the genomic, taxonomic, and genetic content of the famous center.

NCBI kindly provide us the command we should run. It is similar to the following one, which has simply been cut to retrieve only the two files we need to continue with our tutorial.

```bash
#Activate Datasets environment
datasets download genome accession <ACCESSION_NUMBER> --include gff3,genome
```

Since species we want are at least 4, I recommend you to use a for loop to iterate through many accession number, making the process easier and more direct.

## Download using directly FTP

A second more "direct" approach is using the ftp protocol and the command `wget`.

Using NCBI datasets, filter genomes of interest (those already annotated) and list their accession number in a file. By controlling their FTP site, you can double check if annotation is present (`.faa` file should be present, but more in particular you are interested in `.gff`). Manually speaking, from this FTP page you can copy the link of a file you need and yse `wget` to retrieve it.

To ease and make faster the process, you can use the script [download_genome_gff_url_ncbi](./download_genome_gff_url_ncbi.sh) that needs only the abovementioned file listing accession number of interest (optionally you can complete the file with a second column listing the ID you choose for the species and a third with its complete scientific name). This script also needs a file that summarise all assembly present in NCBI that is here `/home/PERSONALE/mirko.martini3/2024/Data/assembly_summary` (originallly, they are located in NCBI FTP site GenBank and RefSeq of folders).

**Remember** not to copy but instead to get used using symlinks (`ln -s <file_to_link> <location_where_to_link>`)

```bash
bash download_genome_gff_url_ncbi.sh <LIST> <ASSEMBLY_SUMMARY_GCF> <ASSEMBLY_SUMMARY_GCA>
```

You can find the script at `/home/PERSONALE/mirko.martini3/00_Lab_CompGeno/2024/05_OG.Inference_Phylogenomic/download_genome_gff_url_ncbi.sh`

## Longest sequence and translation

In a genome, after annotation, it is not uncommon to have multiple isoforms of the same sequence. To not bias your analyses, you can decrease the number of these isoforms using the program [AGAT](https://github.com/NBISweden/AGAT). AGAT works with several different [pearl script](https://agat.readthedocs.io/en/latest/index.html). To keep longest isoforms you will use [agat_sp_keep_longest_isoform.pl](https://agat.readthedocs.io/en/latest/tools/agat_sp_keep_longest_isoform.html). Then, you will extract sequences annotated as genes using [agat_sp_extract_sequences.pl](https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html).

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

**Remember**: files are many. Use a for cycle

### Pseudogene deletion

Even though you cannot untrust everything is made by someone else or by a program, it is ture that sometimes you are true. Annotation files contains mistakes, and sequences labeled with the CDS tag can contain STOP codons inside the sequences. In this case you are looking at a pseudogene, not a gene, and you must delete it from your proteome/genome since you are interested only in working coding sequences.

The first thing to do is to transform your multi line FASTA into a single line FASTA (a common practise since many programs do not recognise multi line FASTA). Specifically in this case it is not demanding since it is the first step performed by the script for psedogene removal.

```bash
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' <multi_line_proteome> > <single_line_proteome>
```

Then run the script [pseudogene_find_eliminate.sh](./pseudogene_find_eliminate.sh), which will do the work for you deleting every sequence that contains a `*`.

```bash
bash pseudogene_find_eliminate.sh
```

You can find the script here `/home/PERSONALE/mirko.martini3/00_Lab_CompGeno/2024/05_OG.Inference_Phylogenomic/pseudogene_find_eliminate.sh`

## Days of Future Past

Since you are already playing with gff file, run also the command

```bash
agat_sp_extract_sequences.pl -g {<GFF_longest_file>} -f {FASTA_file} -t cds --cfs -roo --output {output_file}
```

In this case, you can observe that it is absent the option `-p`, so you will not translate nucleotide sequences into amino acid, but it is present `-roo` (--remove_orf_offset). This flag is indispensable because without it all proteins that present their first codon in a frame different from the first (encoded by AGAT with 0) will be extracted erroneously, with heading nucleotides that will change the frame subsequently elaborated. With the `-p` flag, AGAT automatically takes into consideration the protein phase, but this is not true when `-t CDS` is performed alone.
