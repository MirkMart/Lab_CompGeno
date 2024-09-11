#!/bin/bash
#$1=list to grep in $2 
#$2=where to search


file="$1" ; name=$(cat $file) ; for locus in $name ; do grep -w $locus "$2" ; done > "$1"_out
