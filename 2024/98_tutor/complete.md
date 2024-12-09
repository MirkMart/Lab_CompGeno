# Complete code of the course

## Fastqc: view reads quality

```bash
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
```

## Trimmomatic: remove reads or parts of reads based on the quality score

```bash
trimmomatic PE -threads 20 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/usr/local/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic
```

## KAT hist: compute k-mer frequency
kat hist -o Anoste_kmer27 -t 4 -m 27 Trimmomatic/SRR11672503_1_paired.fastq Trimmomatic/SRR11672503_2_paired.fastq
# Remove comments from Anoste_kmer27 and upload to GenomeScope2 to take the lenght parameter (len) for the following assembly
    #http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=0HnIKI37lVX1nZliWksR

## Assemble reads, construction of contig layout and edge sequences
wtdbg2 -x rs -g GenomeScope_lenght -t 6 -i SRR11672506.fastq.gz -o Anoste #Contig Assembly
wtpoa-cns -t 6 -i .ctg.lay.gz -o Anoste #Consensus Sequences

## BUSCO
conda activate Assembly_tools
busco -i Anoste.raw.fa -m geno -l ../../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Anoste_busco

## Genome polishing
#short-reads(sr)
minimap2 -ax sr --MD -t 6 ../Consensus/Anoste.raw.fa ../../../Fastqc/SRR11672503_1.fastq.gz ../../../Fastqc/SRR11672503_2.fastq.gz > Anoste.raw-sr.sam
#long-reads PacBio (pb)
minimap2 -ax map-pb --MD -t 6 ../Consensus/Anoste.raw.fa ../SRR11672506.fastq.gz > Anoste.raw-pb.sam

##KAT
kat comp -o Anoste.kat -t 6 '../../../Fastqc/SRR11672503_1_paired.fastq.gz ../../../Fastqc/SRR11672503_2_paired.fastq.gz' ../Consensus/Anoste.raw.fa

## Samtools script
SAM=$1
BAM=$2
samtools view -Sb "$SAM" > "$BAM"
rm "$SAM"
samtools sort -@6 -o "${BAM/.bam/.sorted.bam}" $BAM
rm "$BAM"
samtools index "${BAM/.bam/.sorted.bam}"

## Hypo
realpath ../../../Fastqc/SRR11672503_1_paired.fastq.gz > Sr.path
realpath ../../../Fastqc/SRR11672503_2_paired.fastq.gz >> Sr.path

conda activate Mosdepth
mosdepth -n --fast-mode --by 500 Anoste.raw-sr.sorted.wgs Anoste.raw-sr.sorted.bam
zcat Anoste.raw-sr.sorted.wgs.regions.bed.gz | awk '{sum += $4;count++} END {print sum / count}'

conda activate Hypo
hypo -d ../Consensus/Anoste.raw.fa -r @Sr.path -s 224m -c 136 -B ../Polishing/Anoste.raw-pb.sorted.bam -b ../Polishing/Anoste.raw-sr.sorted.bam

# Final Genome statistics
conda activate Assembly_tools
busco -i Hypo.fa -m geno -l ../../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Anoste_busco
kat comp -o Anoste.kat -t 6 '../../../Fastqc/SRR11672503_1_paired.fastq.gz ../../../Fastqc/SRR11672503_2_paired.fastq.gz' Hypo.fa


## Secondary mapping
minimap2 --secondary=no --MD -ax sr -t 6 ../Primary_assembly/Hypo/hypo_Anoste.raw.fasta ../../Fastqc/SRR11672503_1_paired.fastq.gz ../../Fastqc/SRR11672503_2_paired.fastq.gz | samtools view -Sb - > Anoste.corrected.renamed-sr.bam
samtools sort -@10 -o Anoste.corrected.renamed-sr.sorted.bam Anoste.corrected.renamed-sr.bam
rm Anoste.corrected.renamed-sr.bam
samtools index Anoste.corrected.renamed-sr.sorted.bam

## BLAST
blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>
#outfmt: blastn tabular output format 6
#max_target_seqs: 25 target sequences but only 1 hsps so we force the one hsps with minimum e-value


## Contaminats detection - BlobTools Workflow A
conda activate blobtools

../../blobtools/blobtools create -i ../Assembly/Primary_assembly/Hypo/hypo_Anoste.raw.fasta -b ../Assembly/Primary_assembly/Polishing/Anoste.raw-sr.sorted.bam -t hypo_Anoste.raw.blastn -o Anoste.corrected1
../../blobtools/blobtools view -i Anoste.corrected1.blobDB.json -o Anoste-phylum #view the contaminants contigs to remove
../../blobtools/blobtools plot -i Anoste.corrected1.blobDB.json -o Anoste-phylum #view contigs in the plot and search for contaminants

## Removing contaminants filtering by GC content
grep -v "#" Anoste-phylum.Anoste.corrected.blobDB.table.txt | awk '$3 > 0.53' | grep -v "Arthropoda" | cut -f1 > contaminants_contigs.txt
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' ../Assembly/Primary_assembly/Hypo/hypo_Anoste.raw.fasta | grep -w -v -Ff contaminants_contigs.txt - | tr "\t" "\n" > Anoste.corrected-noCont.fa

## Error correction and scaffolding
conda activate Assembly_tools
ragtag.py correct -t 20 <REFERENCE_GENOME> Anoste.corrected-noCont.fa
ragtag.py scaffold -C -t 20 -o ragtag_output/ <REFERENCE_GENOME> ragtag_output/Anoste.correct.fa
    # -C creates a 0-chromosome that comprises all contigs that did not map on the reference assembly, this 0-chromosome will have to be removed
# Reference Genome: /home/PERSONALE/jacopo.martelossi2/Data/Anopheles_stephensi/NCBI_GCF/GCF_013141755.1_UCI_ANSTEP_V1.0_chromosomes.fasta

## Genome Annotation
conda activate MAKER
maker -CTL
#creates three control files (.ctl) 
    # maker_bopts: parameters for the softwares used by MAKER, we mantain default parameters
    # maker_exe: paths to all softwares
    # maker_opts: MAKER run parameters, here we need to modify:
        # genome=path/to/genome (create with realpath)
            # est=path/to/assembled/RNASeq (don't use it)
        # protein=/home/PERSONALE/jacopo.martelossi2/Data/Anopheles_stephensi/Proteomes/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa
        # model_org=*leave blank
        # rmlib=/home/PERSONALE/jacopo.martelossi2/Old_Analyses/Anopheles_stephensi/Genome_annotation/Anoste-families.fa
        # protein2genome=1
            # unmask=1 is used to annotate and compare the unmasked genome to masked genome in order to search for genes that were included in the TE library (don't use it)
        #phred_stats=1 (AED - Annotation Edit Distance: general measure of how well the predicted gene is supported by external evidence (uses Jakkard distance: calculated by overlapping level)
        #min_protein=50
            #alt_splice=1 on,y with RNA-Seq (don't use it)
            #split_hit: if the same protein produces two alignments with a distance greater then split_hit parameter, the two alignments are considered separated hits
maker -base Anoste_rnd1

fasta_merge -d Anoste_rnd1_mAnoster_datastore_index.log
gff3_merge -d Anoste_rnd1_mAnoster_datastore_index.log

awk '$2 == "protein2genome"' Anoste_rnd1.all.gff > protein2genome.gff
awk '$2 == "repeatmasker"' Anoste_rnd1.all.gff > RepeatMasker.gff

conda activate GAAS
agat_sp_statistics.pl --gff Anoste_rnd1.all.gff -o Anoste_rnd1.statistics #Summary statistics of gene models
agat_sq_repeats_analyzer.pl -i -o Anoste_rnd1.repeats.analyzer #Summary statistics of repeats                    ?????

#Launch BUSCO on protein-mode on the predicted proteome
conda activate Assembly_tools
busco -i Anoste_rnd1.all.maker.protein.fasta -m prot -l ../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Anoste.proteins.busco

conda activate MAKER
maker2zff -c 0 -e 0  -l 80 -x 0.1 -d Anoste_rn1.maker.output/Anoste_rnd1_mAnoster_datastore_index.log    #To extract gene models based on mutiple filter criterion
    #genome.dna = sequence
    #genome.ann = gene models
cut -d" " -f5 genome.ann | sort -u | grep -v ">" | wc -l    #count the number of genes, if too much we can put more filters in maker2zff flags(-l, -x, -o)

fathom genome.ann genome.dna -gene-stats > gene.stats   #Print some summary statistics of the selected gene models
fathom genome.ann genome.dna -validate > gene.validate  #Validate gene models and print summary statistics
fathom genome.ann genome.dna -categorize 1000           #Extract gene modeles together with 1000 bp at both ends (flanking sequences) for training
fathom uni.ann uni.dna -export 1000 -plus               #Export and convert uni genes to strand plus

mkdir Params
forge ../export.ann ../export.dna                       #Call the parameter estimation program, better in another directory
cd ..
hmm-assembler.pl Anoste_snap Params/ > Anoste_snap.hmm

mkdir rnd2
ln -s Anoste_rnd1.maker.output/protein2genome.gff
ln -s Anoste_rnd1.maker.output/RepeatMasker.gff
ln -s Anoste_rnd1.maker.output/Anoste.ragtag_scaffolds.chr.fa

maker -CTL
# maker_opts
    #change genome in maker_opts
    #change protein_gff as realpath to protein2gff
    #change rm_gff as realpath to RepeatMask
    #change snaphmm as realpath to Anoste_snap-hmm
    #augustus_species=Anoste
    #model_org=*leave blank
    #change augustus_species=Anoste
    #est2genome=0
    #protein2genome=0
    #min_protein=50
    #AED_threshold=0.5
maker -base Anoste_rnd2

fasta_merge -d Anoste_rnd2_mAnoster_datastore_index.log  # In this case there will be multiple fasta files
gff3_merge -d Anoste_rnd2_mAnoster_datastore_index.log
# At the end of this process we will use the maker.protein and maker.transcript files of proteins predicted by maker

## Final statistics
conda activate GAAS
agat_sp_statistics.pl --gff Anoste_rnd2.all.gff -o Anoste_rnd2.statistics #Summary statistics of gene models
agat_sq_repeats_analyzer.pl  #Summary statistics of repeats

## Launch BUSCO on protein-mode on the predicted proteome
conda activate Assembly_tools
busco -i Anoste_rnd2.all.maker.protein.fasta -m prot -l ../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Anoste.proteins.busco

conda activate MAKER
/home/PERSONALE/jacopo.martelossi2/scripts/AED_cdf_generator.pl -b 0.025 Anoste_rnd2.all.maker.gff > AED.stats #plot in R with ggplot2
    # -b :intervals every 0.025
# Cumulative distribution of AED: what percentage of our dataset has a value equal or lower then an x value (AED)
# AED is better when lower (indicates the quality of annotation)

## Expanding our Protein Dataset from NCBI
wget link_to_NCBI file.protein.faa.gz
wget link_to_NCBI feature_table.txt.gz
gunzip file.protein.faa.gz
gunzip feature_table.txt.gz

## Isoforms will be considered as gene copies, so they have to be removed

```bash
bash /home/PERSONALE/jacopo.martelossi2/scripts/Longest_Isoform.bash feature_table.txt file.protein.faa
#Check feature_count on coding_genes, should be equal to grep -c ">" of protein.NoIsoform.faa
```

Re-format the downloaded proteoms headers keeping the id and the species id_name (es. Anoste|XP_007). For our annotation proteine with an incremental numeration will do it.

```bash
for prote in *.faa; do sed -i -E "s/>(.[^-]+)-(.+) (.+$)/>${prote/.faa/}\|\2/" "$geno"
for geno in *.fna; do sed -i -E "s/>(.[^-]+)-(.+) (.+$)/>${geno/_cds.fna/}\|\2/" "$geno"
```

#Re-format the Anoste_proteom by keeping gene id and adding Anoste
awk '/^>/ {printf ">protein%d\n", ++count; next} {print}' Anoste_anno2.all.maker.proteins.fasta > Anoste.faa


## Orthology inference
conda activate test_env
orthofinder -f proteomes_folder/
#dir1 must be a directory containing all proteoms with an extension of either .fasta or .prot

cd Orthogroups/
#Orthogroups.tsv : For every orthogroup there are proteins associated
#Orthogroups.GeneCount.tsv : Percentage or number of proteins per species

#incrementing single copy orthogroups pruning those multicopy. Then delete those that did not pass our filter and are empty
for tree in *; do echo "$tree"; python3 01_DISCO/disco.py -i "$tree" -o 01_DISCO/${tree/_tree.txt/_disco.nwk} -d "|" -m 3 --remove_in_paralogs --single_tree --keep-labels --verbose 2>/dev/null; done > 01_DISCO/disco.log
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
/home/PERSONALE/dbs/interproscan-5.65-97.0/interproscan.sh -i longest_protein_OGs.fa -goterms -pa -b longest_compgeno/ -cpu 40

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