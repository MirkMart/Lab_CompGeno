#!/bin/bash

#Extract the longest protein for each gene. NB: if more isoforms are of the same length it just take the first occurence.

feature_file=$1
protein_file=$2

####Creation of a summary table with only the longest isoform for gene

cat "$1" | awk 'BEGIN{FS="\t";OFS="\t"}($1=="CDS"){print $11,$13,$16,$19}' | sort -k3,3 -k4,4nr | awk 'BEGIN{FS="\t";OFS="\t";l=0;g=0}{if($3==g) {if ($4>=l) {print $0; g=$3; l=$4}} else {print $0; g=$3; l=$4} }' | awk -F"\t" '!seen[$3]++' > Longest.Isoforms.txt;

####Extraction of the longest isoform from .protein.faa

awk -F"\t" '{print$1}' Longest.Isoforms.txt > Longest_Isoform.tmp
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$2" | grep -Ff Longest_Isoform.tmp - | tr "\t" "\n" > "$2"_NoIsoforms
