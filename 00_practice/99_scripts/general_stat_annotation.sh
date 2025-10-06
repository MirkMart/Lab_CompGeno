#!/bin/bash

#This script summerises general information about the new inferred gff

gff=$1

agat_sp_statistics.pl --gff "$gff"  -o ${gff%.all.gff}_stat.txt
agat_sq_repeats_analyzer.pl -i "$gff" -o ${gff%.all.gff}_repeat_stat.txt
