#!/bin/bash

if [ ! -d "../disco_single_OG" ]; then
  mkdir -p "../disco_single_OG"
fi

# recreate orthogroups using original ones
for tree in *.tre; do
	name=$(basename -s .tre "$tree")
	OG=$(basename -s .tre "$tree")
	while IFS= read -r sequence; do
		grep -A1 "$sequence" ORIGINAL_FOLDER/"$OG".fa >> ../disco_single_OG/"$name".fa
	done < <(grep -o -E "[A-Z][a-z]{4}.[^:]+" "$tree")
done

