# For tutors

## Conda init

Remeber to `conda init` before starting the course, maybe parallele to GitHub. After `init` you need to restart the shell. You need to delete every screen that was created before, because it will be without conda.

## blastn for bloobtool

When blastn for Bloobtool use Sidius.

```bash
blastn -query Anoste_polass.fa -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms title' -max_target_seqs 25 -max_hsps 1 -num_threads 8 -evalue 1e-25 -out Anoste_blast.tsv
```

Available were 8 CPUs. started around 13.20 7/10/2024. End: 21.25 9/10/2025

BLASTDB variable set to `export BLASTDB=/DATABIG/dbs/NCBI/nt`. Only /DATABIG/dbs/NCBI was not sufficient for calling nt. Database last update 22/12/2023

### Liliana

Downloaded database on Liliana `home/PERSONALE/dbs/blastdb`. changed ownership with:

```bash
sudo chown -R 'PERSONALE\mirko.martini3':'PERSONALE\domain^users' blastdb/`
#'' used to take literally what was inside
#PERSONALE\domain^users is the group
```

All databases of intereset must be put in the same folder. This folder does not have to contain sufolders with each database, but all files must be stored together. This is the only way to use directly the name of the database in the command line `nt` or `nr`. Taxdb was automaticcaly dowloaded with nt. The variable `$BLASTDB` is:

```bash
export BLASTDB=/home/PERSONALE/dbs/blastdb
```

## Maker RepeatMasker

During MAKER run it seems to preceed even though there is an import error for the module `h5py` when running `/usr/local/anaconda3/envs/MAKER/share/RepeatMasker/famdb.py`. Indeed even if the conda environment is python2.7, the shebang recalls python3. For this reason I sudo changed it from `#!/usr/bin/env python3` to `#!/usr/local/anaconda3/envs/MAKER/bin/python` (path obtained with `which python` inside MAKER environment).

This did not resolve the problem. There was a new one about syntax. For this reason shebang returned to be the original one, but I deactivated conda and install h5py for python3 of the entire system.

## Install NCBI-Datasets

Then follow the command line suggested by the genome you want to download

## Change headers of Orthogroups without species (ORTHOFINDER REQUIRE SPECIES NAME. It is automatically appended only in tree)

```bash
for fa in *.fa; do for header in $(grep ">" "$fa"); do species=$(grep -oP ".{7}(?=${header/\>/})" ../../00_Results_Dec03/Gene_Trees/${fa/_aligned_output.fa/_tree.txt}); sed -i "s/$header/>$species${header/\>/}/" $fa; done; done
```
