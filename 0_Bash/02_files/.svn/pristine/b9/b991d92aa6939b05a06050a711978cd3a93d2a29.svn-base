#!/bin/bash
#$1=fasta to unwrapp

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $1 > $1".tmp" && mv $1".tmp" $1
