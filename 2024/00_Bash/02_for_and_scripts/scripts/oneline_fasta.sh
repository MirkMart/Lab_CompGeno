#!/bin/bash

#This scripts unwrapp multiline fasta creating oneline files.

#$1=fasta to unwrapp

#The first part works on lines that are not headers---do not start with > (!= ^>)--- and concatenates them. The second part whenever finds an header creates goes to a new line. At the end of the entire script additional new line are added to properly format the file. A the end, the file multiline is replaced with the oneline one having the same name
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $1 > $1".tmp" && mv $1".tmp" $1
