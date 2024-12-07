#!/bin/bash/

# list each plausible pseudogene present. 

mkdir raw_proteomes
mv *.faa raw_proteomes/

#make proteome one line
cd raw_proteomes
for proteome in *.faa; do
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$proteome" > ../${proteome/.faa}".faa"
done
cd ..

#Extract all pseudogene names
mkdir 00_pseudogene_name

for proteome in *.faa; do
	species=$(basename -s .faa "$proteome")
	grep -B1 '*' "$proteome" | grep ">" >> 00_pseudogene_name/"$species"_pseudogenes_name.txt
done

#removes sequences identified as pseudogenes

for pseudo_file in 00_pseudogene_name/*_pseudogenes_name.txt; do
	species=$(basename -s _pseudogenes_name.txt "$pseudo_file")
	while IFS=$'\t' read -r header; do
		sed -E -i "/${header}/{N;d;}" "$species".faa # N option loads the next line found after the pattern and put it into pattern space too; d delete the pattern space
	done < "$pseudo_file" 
done

mv 00_pseudogene_name ../00_genome
