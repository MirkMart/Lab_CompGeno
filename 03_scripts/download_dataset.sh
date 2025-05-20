#!/bin/bash

#This script downlads genomes and relative gff from NCBI datasets creating ready to use folders

AN2name=$1

mkdir 00_genome
mkdir 01_gff

while IFS=$'\t' read -r AN sname ID; do
	echo $AN
	#download specifying the name of the folder
	datasets download genome accession "$AN" --filename "$ID".zip --include genome,gff3
	#unzip specifying the name of the folder
	unzip "$ID".zip -d "$ID"
	#rename the two file of interest
	mv "$ID"/ncbi_dataset/data/"$AN"/*.fna 00_genome/"$ID".fna
	mv "$ID"/ncbi_dataset/data/"$AN"/*.gff 01_gff/"$ID".gff
	#delete the folder
	rm -r "$ID"/
done < "$AN2name"
