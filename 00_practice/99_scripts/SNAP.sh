#!/bin/bash

#this scripts wants to execute all passages needed to compile a SNAP prediction

datastore=$(realpath "$1")
species=$(basename "$datastore" | grep -o -E "[A-Z][a-z]{5}")
number=$(basename "$datastore" | grep -o -E "[0-9]+")

mkdir SNAP
cd SNAP

maker2zff -c 0 -e 0 -l 80 -x 0.1 -d "$datastore"
fathom genome.ann genome.dna -gene-stats > fathom_stats.txt
fathom genome.ann genome.dna -validate 
fathom genome.ann genome.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus

mkdir forge
cd forge

forge ../export.ann ../export.dna

cd ..

hmm-assembler.pl "${species}_prediction${number}" forge > "${species}_prediction${number}.hmm"
