#!/bin/bash

#This script aims to extract a bunch of fasta sequences from a bogger bigger file using headers

#$1=list of loci to grep in $2 
#$2=sequences

#name the variable
list="$1"
#create the name of the output file
name=$(basename -s .txt "$list")

#for each header listed in list, grep the header and the following line from the bigger fasta file and store it in a new file
for locus in $(cat list); do
       grep -w -A1 "$locus" "$2"
done > "$name".fasta
