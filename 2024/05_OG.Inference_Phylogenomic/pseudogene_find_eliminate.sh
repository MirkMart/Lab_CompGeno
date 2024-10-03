#!/bin/bash/

# list each plausible pseudogene present. 
# raw proteomes should have *.fa extension

mkdir 00_raw_proteomes
mv *.fa 00_raw_proteomes/

for proteome in *.faa; do
	species=$(basename -s .faa "$proteome")
	grep -B1 '*' "$proteome" | grep ">" - >> "$species"_pseudogenes_name.txt
done

# This part wants to eliminate sequences that have been defined as pseudogenes genomes 

for pseudo_file in *_pseudogenes_name.txt; do
	species=$(basename -s _pseudogenes_name.txt "$pseudo_file")
	while IFS=$'\t' read -r header; do
		sed -E -i "/${header}/{N;d;}" "$species".faa # N option loads the next line found after the pattern and put it into pattern space too; d delete the pattern space
	done < "$pseudo_file" 
done 

mkdir 01_pseudogene_name
mv *_pseudogenes_name.txt 01_pseudogene_name/