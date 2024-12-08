#!/bin/bash

# this script looks in every OGs looking for the longest protein in order to create a fasta file ready to be used to infer GO and KEGG terms

> longest_protein_OGs.txt

for orthogroup in *_trimmed.fa; do
	orthoname=$(basename -s _trimmed.fa "$orthogroup")
	maxlen=0
	maxname=""
	#create an orthogroup file without any gap and while read it load two different variables that correspond to consecutive lines
	while IFS= read -r name && IFS= read -r sequence; do
		#get the length of the sequence
		lenprotein=${#sequence}
		#if the measured length is greater to the one already stored, substitute the max length and sequence name variable
		if (( lenprotein > maxlen )); then
			maxname=${name#>}
			maxlen=${lenprotein}
		fi
	done < <(sed 's/-//g' "$orthogroup")
	#retrive che complete sequence from where original orthgroups are stored
	sequence=$(grep -A1 "$maxname" <original_orthogroups_folder>${orthoname}.fa | tail -n1)
	printf ">%s@%s\n%s\n" "$orthoname" "$maxname" "$sequence" >> longest_protein_OGs.txt
done
