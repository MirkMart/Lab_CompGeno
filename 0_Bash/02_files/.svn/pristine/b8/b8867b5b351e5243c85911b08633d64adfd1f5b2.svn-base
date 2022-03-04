#!/bin/bash
#$1=list of loci to grep in $2 
#$2=sequences


file="$1" ; name=$(cat $file) ; for locus in $name ; do grep -w -A1 $locus "$2" ; done > "$1".fasta
