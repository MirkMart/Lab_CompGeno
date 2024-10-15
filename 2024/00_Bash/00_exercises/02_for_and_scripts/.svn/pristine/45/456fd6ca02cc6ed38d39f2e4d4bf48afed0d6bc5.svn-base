#!/bin/bash
#$1=fasta file

awk '{if ($1~">") print $0; else print length}' "$1" | awk 'NR%2==1 {printf $0"\t"} NR%2==0 { print $0}' > "$1"_length
