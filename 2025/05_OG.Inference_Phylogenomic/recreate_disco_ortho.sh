#!/bin/bash

#create output folder
if [ ! -d "../01_disco_single_OG" ]; then
  mkdir -p "../01_disco_single_OG"
fi

OG_folder=$1

# recreate orthogroups using original ones
for tree in *.nwk; do
	OG=$(basename -s _disco.nwk "$tree")
	while IFS= read -r sequence; do
		grep -A1 "$sequence" "$OG_folder"/"$OG".fa >> ../01_disco_single_OG/"$OG".fa
	done < <(grep -o -E "[A-Z][a-z]{5}.[^:]+" "$tree")
done

