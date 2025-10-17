# Complete code of the course

> word between |...| in comments defineS the environment where the command was launched

```bash
## Download SRA anda data
#|sequence|
prefetch SRR11672506
fasterq-dump SRR11672506
gzip SRR11672506.fastq
prefetch SRR11672503
fasterq-dump SRR11672503
gzip SRR11672503_1.fastq
gzip SRR11672503_2.fastq

datasets download genome accession GCF_013141755.1 --include gff3,rna,cds,protein,genome,seq-report
unzip ncbi.zip

## Fastqc: view reads quality
#|assembly|
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
multiqc .

## Trimmomatic: remove reads or parts of reads based on the quality score
#|assembly|
trimmomatic PE -threads 20 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log

#which adapter? [adapters illumina](https://github.com/usadellab/Trimmomatic/issues/71)

#Check after trimming
fastqc SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq

## KAT hist: compute k-mer frequency
#|kat|
kat hist -o Anoste_kmer27 -t 4 -m 27 SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq
#Remove comments from Anoste_kmer27 and upload to GenomeScope2 to take the lenght parameter (len) for the following assembly

## Assemble reads, construction of contig layout and edge sequences
#|assembly|
wtdbg2 -x rs -g 227054799 -t 6 -i SRR11672506.fastq.gz -o Anoste #Contig Assembly
wtpoa-cns -t 6 -i Anoste.ctg.lay.gz -o Anoste_raw.fasta #Consensus Sequences

## BUSCO
#|sequence|
export NUMEXPR_MAX_THREADS=80
busco -i ../../03_GenomeAssembly/Anoste_raw.fasta -m geno -l /usr/local/share/busco_databases/diptera_odb12 --cpu 20 -o ./Anoste_raw
busco -i ../../03_GenomeAssembly/Anoste_raw.fasta -m geno -l /usr/local/share/busco_databases/culicidae_odb12 --cpu 20 -o ./Anoste_raw
#|assembly|
assembly-stats Anoste_raw.fasta > Anoste_raw.stats

## KAT
#|kat|
kat comp -o Anoste_raw -t 40 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_raw.fasta


## Genome polishing
#|assembly|
#map reads on assembly
#mapping launched using the script mapping.sh
bash ../../../03_scripts/mapping.sh Anoste_raw.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq SRR11672506.fastq.gz 80

## Hypo
#|assembly|
mosdepth -n --fast-mode -t 80 --by 500 Anoste_raw_sr Anoste_raw_sr_sorted.bam
mosdepth -n --fast-mode -t 80 --by 500 Anoste_raw_pb Anoste_raw_pb_sorted.bam
zcat Anoste_raw_sr.regions.bed.gz | awk '{sum += $4;count++} END {print sum / count}' > coverge_sr.txt

realpath ../../../Fastqc/SRR11672503_1_paired.fastq.gz > Sr.path
realpath ../../../Fastqc/SRR11672503_2_paired.fastq.gz >> Sr.path

R1 = $(realpath ../../../Fastqc/SRR11672503_1_paired.fastq.gz)
R2 = $(realpath ../../../Fastqc/SRR11672503_2_paired.fastq.gz)
echo -e "$R1\n$R2" > Sr.path

hypo -d Anoste_raw.fasta -r @Sr.path -s 227m -c 136 -B Anoste_raw_pb_sorted.bam -b Anoste_raw_sr_sorted.bam -t 80

## Polished Genome statistics
#|sequence|
#genome must be folded because one line fasta interfere with bbtools, blocking BUSCO
fold -w 100 Anoste_pol.fasta > Anoste_fold.fasta #then changed name from fold to pol removing the original one
busco -i ../../03_GenomeAssembly/01_polishing/Anoste_pol.fasta -m geno -l /usr/local/share/busco_databases/culicidae_odb12 --cpu 40 -o ./Anoste_pol
#|kat|
kat comp -o Anoste_pol -t 40 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_pol.fasta
#|assembly|
assembly-stats Anoste_pol.fasta > Anoste_pol.stats

## Second mapping
minimap2 --secondary=no --MD -ax sr -t 6 Anoste_pol.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq | samtools view -Sb - > Anoste_pol.bam
samtools sort -@10 -o Anoste_pol_sorted.bam Anoste_pol.bam
rm Anoste_pol.bam
samtools index Anoste_pol_sorted.bam

## BLAST
#|run in Sidius|
blastn -query Anoste_pol.fasta -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms title' -max_target_seqs 25 -max_hsps 1 -num_threads 10 -evalue 1e-25 -out Anoste_blast.tsv
#outfmt: blastn tabular output format 6
#max_target_seqs: 25 target sequences but only 1 hsps so we force the one hsps with minimum e-value


## Contaminats detection - BlobTools Workflow A

#|sequence|
# Follow [passages](https://github.com/DRL/blobtools) to create node and names dbs
blobtools create -i Anoste_pol.fasta -b Anoste_pol_sorted.bam -t Anoste_blast.tsv -o Anoste
blobtools view -i Anoste.blobDB.json -r phylum -o Anoste_phylum #view the contaminants contigs to remove
# Also family and genus
blobtools plot -i Anoste.blobDB.json -r phylum -o Anoste_phylum #view contigs in the plot and search for contaminants

## Removing contaminants filtering by GC content
grep -v "#" Anoste_phylum.Anoste.blobDB.table.txt | awk '$3 > 0.53' | grep -v "Arthropoda" | cut -f1 > contaminants_contigs.txt
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Anoste_pol.fasta | grep -w -v -Ff contaminants_contigs.txt - | tr "\t" "\n" > Anoste_noContaminants.fasta

#|assembly|
assembly-stats Anoste_noContaminants.fasta > Anoste_noContaminants.stats

## Error correction and scaffolding
#|assembly|
ragtag.py correct -t 50 GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna Anoste_noContaminants.fasta
ragtag.py scaffold -C -t 50 -o ragtag_output/ GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna ragtag_output/ragtag.correct.fasta
# -C creates a 0-chromosome that comprises all contigs that did not map on the reference assembly, this 0-chromosome will have to be removed
# Reference Genome: /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/01_Anoste_reference/GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna

# Create genome with only main chromosomes
grep -A1 "NC_050201.1_RagTag" ragtag.scaffold.fasta >> ../Anoste_chr.fasta
grep -A1 "NC_050202.1_RagTag" ragtag.scaffold.fasta >> ../Anoste_chr.fasta
grep -A1 "NC_050203.1_RagTag" ragtag.scaffold.fasta >> ../Anoste_chr.fasta
fold -w 80 Anoste_chr.fasta > Anoste_chr_fold.fasta #then renamed

## Genome Annotation
#|assembly|
maker -CTL
#creates three control files (.ctl) 
    # maker_bopts: parameters for the softwares used by MAKER, we mantain default parameters
    # maker_exe: paths to all softwares
    # maker_opts: MAKER run parameters, here we need to modify:
        # genome=path/to/genome (create with realpath)
            # est=path/to/assembled/RNASeq (don't use it)
        # protein=/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/02_Annotation/*
        # model_org=*leave blank
        # rmlib=/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/02_Annotation/Anoste_RepeatModeler_library.fa
        # protein2genome=1
            # unmask=1 is used to annotate and compare the unmasked genome to masked genome in order to search for genes that were included in the TE library (don't use it)
        #phred_stats=1 (AED - Annotation Edit Distance: general measure of how well the predicted gene is supported by external evidence (uses Jakkard distance: calculated by overlapping level)
        #min_protein=50
            #alt_splice=1 on,y with RNA-Seq (don't use it)
            #split_hit: if the same protein produces two alignments with a distance greater then split_hit parameter, the two alignments are considered separated hits
maker -base Anoste_rnd1 Anoste_rnd1.ctl

fasta_merge -d Anoste_rnd1_master_datastore_index.log 
gff3_merge -d Anoste_rnd1_master_datastore_index.log

awk '$2 == "protein2genome"' Anoste_rnd1.all.gff > protein2genome.gff
awk '$2 == "repeatmasker"' Anoste_rnd1.all.gff > RepeatMasker.gff

#|GAAS|
agat_sp_statistics.pl --gff Anoste_rnd1.all.gff -o Anoste_rnd1.statistics.txt #Summary statistics of gene models
agat_sq_repeats_analyzer.pl --gff Anoste_rnd1.all.gff -o Anoste_rnd1_repeats.txt #Summary statistics of repeats

#Launch BUSCO on protein-mode on the predicted proteome
busco -i Anoste_rnd1.all.maker.proteins.fasta -m prot -l $BUSCO/culicidae_odb12 --cpu 6 -o Anoste_rnd1

#|assembly|
bash ~/Lab_CompGeno/03_scripts/SNAP.sh Anoste_rnd1_master_datastore_index.log

maker -CTL
# maker_opts
    #change protein_gff as realpath to protein2gff
    #change rm_gff as realpath to RepeatMask
    #change snaphmm as realpath to Anoste_snap-hmm
    #augustus_species=Anoste_cu
    #est2genome=0
    #protein2genome=0
    #AED_threshold=0.5 #if you want to filter models. You can also do it later
maker -base Anoste_rnd2 Anoste_rnd2.ctl

#|GAAS|
gaas_maker_merge_outputs_from_datastore.pl Anoste_rnd2.maker.output/
# In this case there will be multiple fasta files
# At the end of this process we will use the maker.protein and maker.transcript files of proteins predicted by maker

## Final statistics
#|GAAS|
agat_sq_repeats_analyzer.pl --gff maker_mix.gff -o Anoste_rnd2_repeats.txt #Summary statistics of repeats
AED_cdf_generator.pl -b 0.025 maker_mix.gff > AED_maker_mix.stats #plot in R with ggplot2
# Cumulative distribution of AED: what percentage of our dataset has a value equal or lower then an x value (AED)
# AED is better when lower (indicates the quality of annotation)

## Launch BUSCO on protein-mode on the predicted proteome
#|sequence|
busco -i maker_annotation.proteins.fasta -m prot -l $BUSCO/culicidae_odb12 --cpu 6 -o Anoste_rnd2/

## Expanding our Protein Dataset from NCBI
# Create dataset "AN    sname   ID"
#|sequence|
bash ../99_scripts/download_dataset.sh dataset.txt 
#|GAAS|
# extract and translate longest isoform
bash ../../99_scripts/AGAT_longest_extract.sh
# remove speudogenes if present
bash ../../99_scripts/pseudogene_find_eliminate.sh
# extract nt CDS if ω will be computed
for i in *; do agat_sp_extract_sequences.pl -g $i -f ../../00_genome/${i/_longest.gff/.fna} -t cds --cfs -roo --output ../../02_proteome/01_nt/${i/_longest.gff/.fna} & done

#Re-format the downloaded proteoms headers keeping the id and the species id_name (es. Anoste|XP_007). For our annotation proteine with an incremental numeration will do it.
for prote in *.faa; do sed -i -E "s/(>.[^ ]+) gene\=gene-(.[^ ]+) (.+$)/>"${prote/.faa/}"\|\2/" "$prote"; done
for geno in *.fna; do sed -i -E "s/(>.[^ ]+) gene\=gene-(.[^ ]+) (.+$)/>"${geno/.fna/}"\|\2/" "$geno"; done

#Re-format the Anoste_proteome by keeping gene id and adding Anoste
awk '/^>/ {printf ">Anoste|protein%d\n", ++count; next} {print}' Anoste.faa > Anoste_changed
mv Anoste_changed Anoste.faa
awk '/^>/ {printf ">Anoste|protein%d\n", ++count; next} {print}' ../../../../04_GenomeAnnotation/01_round2/maker_output_processed_Anoste_rnd2/maker_annotation.transcripts.fasta > Anoste.fna

## Orthology inference
#|orthofinder|
orthofinder -t 60 -a 60 -f ../00_data/03_dataset/02_proteome/
#dir must be a directory containing all proteoms with an extension of either .fasta, .fa, .faa. etc (see orthofinder [github](https://github.com/davidemms/OrthoFinder))

cd Orthogroups/
#Orthogroups.tsv : For every orthogroup there are proteins associated
#Orthogroups.GeneCount.tsv : Percentage or number of proteins per species

#incrementing single copy orthogroups pruning those multicopy. Then delete those that did not pass our filter and are empty
while IFS=' ' read -r OG tree; do python3 /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/99_scripts/disco.py -i <(echo "$tree") -o ../../../01_DISCO/${OG/:/}.nwk -d "|" -m 4 --remove_in_paralogs --keep-labels --verbose > disco.log; done < <(sed -E 's/[A-Z][a-z]{5}_//g; s/n[0-9]*+//g' Resolved_Gene_Trees.txt)
find . -type f -empty -delete
#delete those that are already known to be single copy complete
for fa in 00_orthofinder/Results_Nov30/Single_Copy_Orthologue_Sequences/*.fa; do name=$(basename $fa); rm ${name/.fa/_disco.nwk}; done
#recreate new orthogroups
bash /home/PERSONALE/mirko.martini3/00_Lab_CompGeno/2024/05_OG.Inference_Phylogenomic/recreate_disco_ortho.sh ../../00_orthofinder/Results_Nov30/Orthogroup_Sequences

#Align single copy orthogroups
for fa in *; do mafft --auto --anysymbol "$fa" > 02_aligned/00_single_complete/${fa/.fa/_aligned.fa}

for fa in *; do bmge -i "$fa" -t AA -m BLOSUM30 -e 0.5 -g 0.4 -of 03_trimmed/00_fasta/00_single_complete/${fa/_aligned.fa/_trimmed.fa} -oh 03_trimmed/01_html/${fa/_aligned.fa/_trimmed.html}; done
#Trimming reduces alignment dimension, reducing computational time

# NCBI alignment viewer: used to see alignment coverage (view variable regions and core regions)
    #View one trimmed and one non-trimmed sequence: After trimming all major drops in coverage are removed

# Concatenating protein fasta
for fa in *; do sed -i -E 's/\|.+$//g' "$fa"; done
python3 AMAS.py concat -y nexus -i *.fa -f fasta -d aa -c 10 -t cosncatenated.nexus
#Creates a concatenated file and a partitioning file in nexus format(.txt)

## Phylogenetic inference (~40 hours with 25 CPUs)
iqtree -m TESTNEW -bb 1000 -s concatenated.nexus --prefix species_tree -nt AUTO    
    # -p: optional if you want to infere based on partition (one per each gene)
    # -bb 1000: one thousand of bootstrap, although this is ultrafast bootstrap so you need at least 1000 replicates (100 with std bootstrap)

# species_tree.treefile (newick format)
    # itol online (https://itol.embl.de/upload.cgi) to view treefile file

## Divergence Time Estimation
#Find secondary calibration information on timetree (https://timetree.org/) using adjusted_time
#Create time estimation file
vi calibration.txt

iqtree -s ../../02_Orthogroups/03_trimmed/00_fasta/concatenated.nexus --date calibration.txt --date-tip 0 -o Dromel -m TESTNEW -nt 6 --prefix time_tree --date-options "-u 1"
    #-u : branches shorter than 1 collapse in a politomy

## CAFE
mkdir CAFE
cp Orthofinder/Orthogroups/Orthogroups.GeneCount.tsv
sed 's/^/NONE\t/g' Orthogroups.GeneCount.tsv | cut -f1,2,3,4,5,6,7,8 > GeneCount.CAFE.tsv
#remove 'protein suffix' in nano text editor

#Export nexus time-tree in newick format from ITOL
cp TimeTree.nwk

conda activate CAFE
/home/PERSONALE/jacopo.martelossi2/Software/CAFE5/bin/cafe5 --infile GeneCounts.CAFE.tsv --tree TimeTree.nwk -e -o CAFE_Error
#CAFE_Error/Base_error_model.txt: 3rd line, value on the right is a percentage number (0.01) that represents the error model stimated
/home/PERSONALE/jacopo.martelossi2/Software/CAFE5/bin/cafe5 --infile GeneCounts.CAFE.tsv --tree TimeTree.nwk -eCAFE_Error/Base_error_model.txt -o CAFE_Error_Corrected
#for summarise all results for multiple CAFE analysis in order to confirm convergence
for folder in */; do lnL=$(grep "lnL" ${folder}/Base_results.txt | grep -oE "[0-9]*\.[0-9]*"); L=$(grep "Lambda" ${folder}/Base_results.txt | grep -oE "[0-9]*\.[0-9]*"); E=$(grep "Epsilon" ${folder}/Base_results.txt | grep -oE "[0-9]*\.[0-9]*"); echo -e "$lnL\t$L\t$E" >> 1L0G_results.tsv; done
#Base_clade_results: for each node, how many gene famlies increase and how many decrease
#Base_asr.tre -> ITOL -> Advanced: Node IDs: Display: to see which number corresponds to which node
#to extract only tree with significant changes. This file can be visualised in FigTree
echo $'#nexus\nbegin trees;' > Significant_trees.tre
grep "*" Gamma_asr.tre >> Significant_trees.tre
echo "end;">>Significant_trees.tre
#Base_family_results.txt: family id - p-value - is there a significant change in the branch?
cut -f3 Base_family_results.txt | sort | unique -c #to see how many families have significant changes

grep "Anoste<5>\*" Base_asr.tre | wc -l #to find all significant modifications in Anoste branch
varAnoste=$(grep "<5>\*" Base_asr.tre | cut -d ":" -f8 | cut -d"_" -f2)
varAaeg=$(grep "<1>\*" Base_asr.tre | cut -d ":" -f1 | cut -d"_" -f2)
varAsub=$(grep "<2>\*" Base_asr.tre | cut -d ":" -f2 | cut -d"_" -f2)
varCqui=$(grep "<4>\*" Base_asr.tre | cut -d ":" -f6 | cut -d"_" -f2)
varWsmi=$(grep "<3>\*" Base_asr.tre | cut -d ":" -f4 | cut -d"_" -f2)

var8=$(grep "<5>\*" Base_asr.tre | cut -d":" -f9 | cut -d"_" -f2) #previous node (8) value for significant families
varOG=$(grep "<5>\*" Base_asr.tre | cut -d" " -f4)

pAnoste <(printf %s "$varOG") <(printf %s "$varAnoste") <(printf %s "$var8") > Anoste_Significant.txt
#upload to R and create another column that subtracts the first and the second column: based on if it's negative or positive it is a contraction or expansion
awk '{print $0, $2 - $3}' Anoste_Significant.txt > Anoste_Significant.txt.diff
sort -k4,4n Aaeg_Significant.txt.diff | grep -v "-" | awk '{print $1}' > Anoste_Significant.GO
while read line; do grep "Anoste" Orthofinder/orthogroup_Sequences/"$line"*.fa; done < Anoste_Significant.GO > Anoste_Significant_Genes.txt

#GO annotation with interproscan. Performed on 134.142.204.241
interproscan.sh -i longest_protein_OGs.fa -goterms -pa -cpu 40 -appl PANTHER

#Panzer2: reasearch by homology of GO-Terms (http://ekhidna2.biocenter.helsinki.fi/sanspanz/)
    #upload all Anoste proteins
    #GO prediction: Arthropoda

http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//VvjXh3TeoFe/index.html
wget GO.out

#Add ":" before each GO id
awk '{print $1, "GO:"$3}' GO.out > Gene_GO.tsv

#Sort GO id for each Gene
awk '{ a[$1]=a[$1]","$2; } END { for (i in a) {sub(/,/,"",a[i]);printf "%s %s\n",i,a[i] } }' Gene_GO.tsv > Gene_GO_annotation.tsv

## Retrotranslation

Rename orthogroups using precedent headers

```bash
for fa in *.fa; do for header in $(grep ">" "$fa"); do new=$(grep "$header" ../../../00_orthofinder/Results_Nov30/Single_Copy_Orthologue_Sequences/${fa/_trimmed/}); sed -i.old "s/${header}/${new}/" "$fa"; done; done
```

nucleotide orthogroups (both from single copy and DISCO)

```bash
for trimmed in *.fa; do bash ../../../create_nucleo_orthogroups.sh "$trimmed" ../../../04_cds_orthogroups /home/PERSONALE/mirko.martini3/01_2024/00_Data/05_CDS/03_nucleo_cds; done
```

## R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install('topGO', force = TRUE)

library(topGO)
library(tidyr)
library(plyr)

geneID2GO <- readMappings(file = "Aaeg_GO_annotation.tsv")
geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("Aaeg_significant_genes.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5)

Expanded_BP_classic_fisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")

Expanded_BP_allRes <- GenTable(myGOdata, 
                               Classic_Fisher = Expanded_BP_classic_fisher,
                               orderBy = "Classic_Fisher", topNodes=1000)

Expanded_BP_allRes$Classic_Fisher <- as.numeric(Expanded_BP_allRes$Classic_Fisher)
Expanded_BP_allRes <- subset(Expanded_BP_allRes, Classic_Fisher < 0.05)

#Basic Bar Plot: top 20
top10_BP_allRes <- head(Expanded_BP_allRes, 20)

allRes_df <- as.data.frame(top10_BP_allRes)

ggplot(allRes_df, aes(x = reorder(Term, -log10(Classic_Fisher)), y = -log10(Classic_Fisher))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Enriched GO Terms of Aedes aegypti",
       x = "GO Term",
       y = "-log10(p-value)") +
  theme_minimal()