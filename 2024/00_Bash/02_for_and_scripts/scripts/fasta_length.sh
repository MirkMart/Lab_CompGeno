#!/bin/bash

#This script aims to count the legth of FASTA sequences

#$1=fasta file

#The first part check each line for the presence of ">". If it is present, awk prints the entire line. Otherwise, awk counts the characters
#The second part performs a different task whether the line number is odd or even. If odd prints again the entire line with a trailing tab (no more a new line); if it is even simply prints the line next to the first

name=$(basename -s .fa "$1")
awk '{if ($1~">") print $0; else print length}' "$1" | awk 'NR%2==1 {printf $0"\t"} NR%2==0 { print $0}' > "$name"_length.tsv
