# For tutors

## Conda init

Remeber to `conda init` before starting the course, maybe parallel to GitHub. After `init` you need to restart the shell. You need to delete every screen that was created before, because it will be without conda.

## blastn for bloobtool

When blastn for Bloobtool use Sidius (137.204.142.152).

```bash
export BLASTDB=/DATABIG/dbs/NCBI/nt
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

True in the precedent server:

> During MAKER run it seems to preceed even though there is an import error for the module `h5py` when running `/usr/local/anaconda3/envs/MAKER/share/RepeatMasker/famdb.py`. Indeed even if the conda environment is python2.7, the shebang recalls python3. For this reason I sudo changed it from `#!/usr/bin/env python3` to `#!/usr/local/anaconda3/envs/MAKER/bin/python` (path obtained with `which python` inside MAKER environment).
> This did not resolve the problem. There was a new one about syntax. For this reason shebang returned to be the original one, but I deactivated conda and install h5py for python3 of the entire system.

In any case, installed h5py in the system even this time (without any conda and in the sequence conda). This has been done before any error was printed in the stdout.

```bash
sudo apt install python3-h5py
mamba install h5py #h5py already present in the environment
```

## Install NCBI-Datasets

Then follow the command line suggested by the genome you want to download

## Change headers of Orthogroups without species (ORTHOFINDER REQUIRE SPECIES NAME. It is automatically appended only in tree)

```bash
for fa in *.fa; do for header in $(grep ">" "$fa"); do species=$(grep -oP ".{7}(?=${header/\>/})" ../../00_Results_Dec03/Gene_Trees/${fa/_aligned_output.fa/_tree.txt}); sed -i "s/$header/>$species${header/\>/}/" $fa; done; done
```

## Augustus

Command runs with sudo permision because it needs to write here `/opt/miniforge3/envs/sequence/config/species`. Probably this folder must be moved inside the environment with maker (NOT SEQUENCE)

```bash
busco -i ../../03_GenomeAssembly/03_scaffolding/Anoste_chr.fasta -c 30 -l /usr/local/share/busco_databases/culicidae_odb12 --augustus --long -m genome --out Anoste_BuscoTraining_cu --augustus_parameters='--progress=true'
```

Options:

- `--long`: helps Augustus improving its performance even if adds time to the computation. Useful for non-model species.
- `--augustus`: activates Augustus annotation.
- `--progress=true`: show a progress meter.

Once BUSCO has finished, inside the folders of outputs there will be a folder named 'retraining_parameters' (Anoste_BuscoTraining_cu/run_culicidae_odb12/augustus_output/retraining_parameters). Inside this forder there is 'BUSCO_Anoste_BuscoTraining_cu' with are all the files neeed to create a new reference species in augustus config folder. Thus I copied them in '/opt/miniforge3/envs/assembly/config/species/Anoste_cu' creating this new model species.

Lastly, change the name of the old files in the configuration file.

```bash
sed -i 's/BUSCO_Anoste_BuscoTraining_cu/Anoste_cu/' Anoste_cu_parameters.cfg
```

## Transcriptome

Good [resource](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004393).

SRR_accession_list is obtained from [SRA run selector](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP585377&o=acc_s%3Aa&s=SRR33576268,SRR33576267,SRR33576266,SRR33576264,SRR33576263,SRR33576265).

```bash
while read SRR; do prefetch "$SRR"; fastq-dump --split-3 "$SRR"; done < SRR_accession_list.txt #split is not default differently than as said
# Quality control reads
for i in *.fastq; do fastqc $i & done
# Sum into a single report
multiqc .
# Trim
for i in *_1.fastq; do name=$(sed 's/_1.fastq//' <<< "$i"); trimmomatic PE -threads 80 -phred33 "$name"_1.fastq "$name"_2.fastq "$name"_1_paired.fastq "$name"_1_unpaired.fastq "$name"_2_paired.fastq "$name"_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:9 2> "$name"_trimmomatic.stats; done
```
