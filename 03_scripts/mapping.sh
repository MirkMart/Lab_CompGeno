#!/bin/bash

#This script aims to automatically map short reads on an assembly

#variables
assembly=$1
fastq1=$2
fastq2=$3
longread=$4
cores=$5

#name creation
outname=$(basename -s .fasta "$assembly")

#mapping short reads on assembly
#minimap2 -ax sr --MD -t "$cores" "$assembly" "$fastq1" "$fastq2" > "$outname"_sr.sam
#minimap2 -ax map-pb --MD -t "$cores" "$assembly" "$longread"  > "$outname"_pb.sam

#converting sam into bam
samtools view -Sb "$outname"_sr.sam > "$outname"_sr.bam
samtools view -Sb "$outname"_pb.sam > "$outname"_pb.bam

rm "$outname"_sr.sam
rm "$outname"_pb.sam

#sorting and indexing bam file
samtools sort -@"$cores" -o "$outname"_sr_sorted.bam "$outname"_sr.bam
samtools sort -@"$cores" -o "$outname"_pb_sorted.bam "$outname"_pb.bam

samtools index "$outname"_sr_sorted.bam
samtools index "$outname"_pb_sorted.bam

rm "$outname"_sr.bam
rm "$outname"_pb.bam

