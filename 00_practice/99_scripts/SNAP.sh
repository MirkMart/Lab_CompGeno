#!/bin/bash

#this scripts wants to execute all passages needed to compile a SNAP prediction

datastore=$(realpath "$1")
species=$(basename "$datastore" | grep -o -E "[A-Z][a-z]{5}")
number=$(basename "$datastore" | grep -o -E "[0-9]+")

mkdir SNAP
cd SNAP

maker2zff -c 0 -e 0 -l 80 -x 0.1 -d "$datastore"            #To extract gene models based on mutiple filter criterion
                                                            #genome.dna = sequence
                                                            #genome.ann = gene models
fathom genome.ann genome.dna -gene-stats > fathom_stats.txt #Print some summary statistics of the selected gene models
fathom genome.ann genome.dna -validate                      #Validate gene models and print summary statistics
fathom genome.ann genome.dna -categorize 1000               #Extract gene modeles together with 1000 bp at both ends (flanking sequences) for training
fathom uni.ann uni.dna -export 1000 -plus                   #Export and convert uni genes to strand plus

mkdir forge
cd forge

forge ../export.ann ../export.dna                           #Call the parameter estimation program, better in another directory   

cd ..

hmm-assembler.pl "${species}_prediction${number}" forge > "${species}_prediction${number}.hmm"
