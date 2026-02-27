#!/bin/bash

# this script aims to automatically execute many codeml analyses. It needs the aligned and trimmed sequence of the orthogroup and the two trees that indicate how many omega values to compute.

sequence=$1
rawconfig=$2

#create variable, modify name and path.
ortho=$(basename -s .fna "$sequence")
#the second part of this sed is needed to use it in the following sed command, in order to escale each / which, in turn, could create incomprehansion with the command
sequence=$(realpath "$sequence" | sed 's/\//\\\//g')
tree1=$(realpath tree1/"$ortho".nwk | sed 's/\//\\\//g')
tree2=$(realpath tree2/"$ortho".nwk | sed 's/\//\\\//g')

#create control file for each model we want to test. Create intermediate folders if they are not present
if [ ! -d "ctl" ]; then
  mkdir -p "ctl"
fi

if [ ! -d "codeml_results" ]; then
  mkdir -p "codeml_results"
fi

echo "Compiling ${ortho} configuration files" 
sed -E "s/SEQUENCE/$sequence/; s/TREEFILE/$tree1/; s/OUTFILE/codeml_results\/${ortho}_codeml1.txt/" "$rawconfig" > ctl/"$ortho"_1.ctl
sed -E "s/SEQUENCE/$sequence/; s/TREEFILE/$tree2/; s/OUTFILE/codeml_results\/${ortho}_codeml2.txt/" "$rawconfig" > ctl/"$ortho"_2.ctl

source activate IQ_tree

#run the codeml analysis
if [ ! -d "log" ]; then
  mkdir -p "log"
fi

echo "Running codeml model 1 on ${ortho}"
codeml ctl/"$ortho"_1.ctl > log/"$ortho"_1.log
echo "Running codeml model 2 on ${ortho}"
codeml ctl/"$ortho"_2.ctl > log/"$ortho"_2.log

source activate Python

#run the likelihood ratio test in order to understand which model was the best
echo "Running LRT on ${ortho}"
python3 LRT.py codeml_results/${ortho}_codeml1.txt codeml_results/${ortho}_codeml2.txt

echo "${ortho} analyses COMPLETED!"