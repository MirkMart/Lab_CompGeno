#!/bin/bash

#This script wants to collect gff and fasta files after a genome annotation run with MAKER

maker_datastore=$1

fasta_merge -d "$maker_datastore"
gff3_merge -d "$maker_datastore"
