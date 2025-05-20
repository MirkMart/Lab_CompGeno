#!/bin/bash

#compute longest isoform and extract from the genome

mkdir 01_gff/00_longest

for gff in 01_gff/*.gff; do 
	name=$(basename -s .gff "$gff")
	agat_sp_keep_longest_isoform.pl --gff "$gff" -o 01_gff/00_longest/"$name"_longest.gff
done

mkdir 02_proteome

for gff in 01_gff/00_longest/*_longest.gff; do
	name=$(basename -s _longest.gff "$gff")
	agat_sp_extract_sequences.pl -g "$gff" -f 00_genome/"$name".fna -t cds -p --cfs --output 02_proteome/"$name".faa
done

mkdir 00_genome/log
mv 00_genome/*.pag 00_genome/log
mv00_genome/*dir 00_genome/log
