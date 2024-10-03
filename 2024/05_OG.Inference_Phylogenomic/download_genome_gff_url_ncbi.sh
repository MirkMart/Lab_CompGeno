#!/bin/bash/
set -e
set -u
set -o pipefail

# the script donwloads genomes through built url from assembly_summary.txt
# species.txt = list of complete species names (Bufo bufo)
# list = list of  assembly accession

#species_or_accession="$1"
GCF_summary="$2"
GCA_summary="$3"

#download genomes using accession number and assembly_summary. If provided scientific name change ("$1" ~ SPECIES) in ("$8" ~ SPECIES)

> url.txt

# GCA or GCF genome
while IFS=$'\t' read -r AN ID name; do
	if [[ $AN == GCF* ]]; then
		awk -v SPECIES=^"$AN" 'BEGIN{FS="\t"}{if($1 ~ SPECIES){print $20}}' "$GCF_summary" | awk 'BEGIN{OFS=FS="/"}{print ""$0,$NF"_genomic.fna.gz"}'
	elif [[ $AN == GCA* ]]; then
		awk -v SPECIES=^"$AN" 'BEGIN{FS="\t"}{if($1 ~ SPECIES){print $20}}' "$GCA_summary" | awk 'BEGIN{OFS=FS="/"}{print ""$0,$NF"_genomic.fna.gz"}'
	fi
done < "$1" >> url.txt 

# gff
while IFS=$'\t' read -r AN ID name; do	
	if [[ $AN == GCF* ]]; then
		awk -v list="$AN" 'BEGIN{FS="\t"}{if($1==list) {print $20}}' "$GCF_summary" | sed -r 's/(GC[AF]_[0-9.]*_.*$)/\1\/\1_genomic.gff.gz/g'
	elif [[ $AN == GCA* ]]; then
		awk -v list="$AN" 'BEGIN{FS="\t"}{if($1==list) {print $20}}' "$GCA_summary" | sed -r 's/(GC[AF]_[0-9.]*_.*$)/\1\/\1_genomic.gff.gz/g'
	fi
done < "$1" >> url.txt

for url in $(cat url.txt); do	
	wget "$url"
done

rm url.txt

