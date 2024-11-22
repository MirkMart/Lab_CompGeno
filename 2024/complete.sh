#!/bin/bash

## Fastqc: view reads quality
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz

## Trimmomatic: remove reads or parts of reads based on the quality score
trimmomatic PE -threads 20 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/usr/local/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic

## KAT hist: compute k-mer frequency
kat hist -o Aste_kmer27 -t 4 -m 27 Trimmomatic/SRR11672503_1_paired.fastq Trimmomatic/SRR11672503_2_paired.fastq
# Remove comments from Aste_kmer27 and upload to GenomeScope2 to take the lenght parameter (len) for the following assembly
    #http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=0HnIKI37lVX1nZliWksR

## Assemble reads, construction of contig layout and edge sequences
wtdbg2 -x rs -g GenomeScope_lenght -t 6 -i SRR11672506.fastq.gz -o Aste #Contig Assembly
wtpoa-cns -t 6 -i .ctg.lay.gz -o Aste #Consensus Sequences

## BUSCO
conda activate Assembly_tools
busco -i Aste.raw.fa -m geno -l ../../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Aste_busco

## Genome polishing
#short-reads(sr)
minimap2 -ax sr --MD -t 6 ../Consensus/Aste.raw.fa ../../../Fastqc/SRR11672503_1.fastq.gz ../../../Fastqc/SRR11672503_2.fastq.gz > Aste.raw-sr.sam
#long-reads PacBio (pb)
minimap2 -ax map-pb --MD -t 6 ../Consensus/Aste.raw.fa ../SRR11672506.fastq.gz > Aste.raw-pb.sam

##KAT
kat comp -o Aste.kat -t 6 '../../../Fastqc/SRR11672503_1_paired.fastq.gz ../../../Fastqc/SRR11672503_2_paired.fastq.gz' ../Consensus/Aste.raw.fa

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
mosdepth -n --fast-mode --by 500 Aste.raw-sr.sorted.wgs Aste.raw-sr.sorted.bam
zcat Aste.raw-sr.sorted.wgs.regions.bed.gz | awk '{sum += $4;count++} END {print sum / count}'

conda activate Hypo
hypo -d ../Consensus/Aste.raw.fa -r @Sr.path -s 224m -c 136 -B ../Polishing/Aste.raw-pb.sorted.bam -b ../Polishing/Aste.raw-sr.sorted.bam

# Final Genome statistics
conda activate Assembly_tools
busco -i Hypo.fa -m geno -l ../../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Aste_busco
kat comp -o Aste.kat -t 6 '../../../Fastqc/SRR11672503_1_paired.fastq.gz ../../../Fastqc/SRR11672503_2_paired.fastq.gz' Hypo.fa


## Secondary mapping
minimap2 --secondary=no --MD -ax sr -t 6 ../Primary_assembly/Hypo/hypo_Aste.raw.fasta ../../Fastqc/SRR11672503_1_paired.fastq.gz ../../Fastqc/SRR11672503_2_paired.fastq.gz | samtools view -Sb - > Aste.corrected.renamed-sr.bam
samtools sort -@10 -o Aste.corrected.renamed-sr.sorted.bam Aste.corrected.renamed-sr.bam
rm Aste.corrected.renamed-sr.bam
samtools index Aste.corrected.renamed-sr.sorted.bam

## BLAST
blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>
#outfmt: blastn tabular output format 6
#max_target_seqs: 25 target sequences but only 1 hsps so we force the one hsps with minimum e-value


## Contaminats detection - BlobTools Workflow A
conda activate blobtools

../../blobtools/blobtools create -i ../Assembly/Primary_assembly/Hypo/hypo_Aste.raw.fasta -b ../Assembly/Primary_assembly/Polishing/Aste.raw-sr.sorted.bam -t hypo_Aste.raw.blastn -o Aste.corrected1
../../blobtools/blobtools view -i Aste.corrected1.blobDB.json -o Aste-phylum #view the contaminants contigs to remove
../../blobtools/blobtools plot -i Aste.corrected1.blobDB.json -o Aste-phylum #view contigs in the plot and search for contaminants

## Removing contaminants filtering by GC content
grep -v "#" Aste-phylum.Aste.corrected.blobDB.table.txt | awk '$3 > 0.53' | grep -v "Arthropoda" | cut -f1 > contaminants_contigs.txt
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' ../Assembly/Primary_assembly/Hypo/hypo_Aste.raw.fasta | grep -w -v -Ff contaminants_contigs.txt - | tr "\t" "\n" > Aste.corrected-noCont.fa

## Error correction and scaffolding
conda activate Assembly_tools
ragtag.py correct -t 20 <REFERENCE_GENOME> Aste.corrected-noCont.fa
ragtag.py scaffold -C -t 20 -o ragtag_output/ <REFERENCE_GENOME> ragtag_output/Aste.correct.fa
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
        # rmlib=/home/PERSONALE/jacopo.martelossi2/Old_Analyses/Anopheles_stephensi/Genome_annotation/Aste-families.fa
        # protein2genome=1
            # unmask=1 is used to annotate and compare the unmasked genome to masked genome in order to search for genes that were included in the TE library (don't use it)
        #phred_stats=1 (AED - Annotation Edit Distance: general measure of how well the predicted gene is supported by external evidence (uses Jakkard distance: calculated by overlapping level)
        #min_protein=50
            #alt_splice=1 on,y with RNA-Seq (don't use it)
            #split_hit: if the same protein produces two alignments with a distance greater then split_hit parameter, the two alignments are considered separated hits
maker -base Aste_rnd1

fasta_merge -d Aste_rnd1_master_datastore_index.log
gff3_merge -d Aste_rnd1_master_datastore_index.log

awk '$2 == "protein2genome"' Aste_rnd1.all.gff > protein2genome.gff
awk '$2 == "repeatmasker"' Aste_rnd1.all.gff > RepeatMasker.gff

conda activate GAAS
agat_sp_statistics.pl --gff Aste_rnd1.all.gff -o Aste_rnd1.statistics #Summary statistics of gene models
agat_sq_repeats_analyzer.pl -i -o Aste_rnd1.repeats.analyzer #Summary statistics of repeats                    ?????

#Launch BUSCO on protein-mode on the predicted proteome
conda activate Assembly_tools
busco -i Aste_rnd1.all.maker.protein.fasta -m prot -l ../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Aste.proteins.busco

conda activate MAKER
maker2zff -c 0 -e 0  -l 80 -x 0.1 -d Aste_rn1.maker.output/Aste_rnd1_master_datastore_index.log    #To extract gene models based on mutiple filter criterion
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
hmm-assembler.pl Aste_snap Params/ > Aste_snap.hmm

mkdir rnd2
ln -s Aste_rnd1.maker.output/protein2genome.gff
ln -s Aste_rnd1.maker.output/RepeatMasker.gff
ln -s Aste_rnd1.maker.output/Aste.ragtag_scaffolds.chr.fa

maker -CTL
# maker_opts
    #change genome in maker_opts
    #change protein_gff as realpath to protein2gff
    #change rm_gff as realpath to RepeatMask
    #change snaphmm as realpath to Aste_snap-hmm
    #augustus_species=Aste
    #model_org=*leave blank
    #change augustus_species=Aste
    #est2genome=0
    #protein2genome=0
    #min_protein=50
    #AED_threshold=0.5
maker -base Aste_rnd2

fasta_merge -d Aste_rnd2_master_datastore_index.log  # In this case there will be multiple fasta files
gff3_merge -d Aste_rnd2_master_datastore_index.log
# At the end of this process we will use the maker.protein and maker.transcript files of proteins predicted by maker

#Final statistics
conda activate GAAS
agat_sp_statistics.pl --gff Aste_rnd2.all.gff -o Aste_rnd2.statistics #Summary statistics of gene models
agat_sq_repeats_analyzer.pl  #Summary statistics of repeats

#Launch BUSCO on protein-mode on the predicted proteome
conda activate Assembly_tools
busco -i Aste_rnd2.all.maker.protein.fasta -m prot -l ../../../../../PERSONALE/jacopo.martelossi2/dbs/diptera_odb10 --cpu 6 -o Aste.proteins.busco

conda activate MAKER
/home/PERSONALE/jacopo.martelossi2/scripts/AED_cdf_generator.pl -b 0.025 Aste_rnd2.all.maker.gff > AED.stats #plot in R with ggplot2
    # -b :intervals every 0.025
# Cumulative distribution of AED: what percentage of our dataset has a value equal or lower then an x value (AED)
# AED is better when lower (indicates the quality of annotation)

## Expanding our Protein Dataset from NCBI
wget link_to_NCBI file.protein.faa.gz
wget link_to_NCBI feature_table.txt.gz
gunzip file.protein.faa.gz
gunzip feature_table.txt.gz

#Isoforms will be considered as gene copies, so they have to be removed
bash /home/PERSONALE/jacopo.martelossi2/scripts/Longest_Isoform.bash feature_table.txt file.protein.faa
#Check feature_count on coding_genes, should be equal to grep -c ">" of protein.NoIsoform.faa

#Re-format the downloaded proteoms headers keeping the id and the species id_name (es. XP_007 | Aalb)
sed -E 's/^>([^ ]+).*/>\1|Species/' your_input.fasta > modified_output.fasta

#Re-format the Aste_proteom by keeping gene id and adding Aste
sed 's/>augustus_masked-\(.*\)-gene-\(.*\)-mRNA-1 protein/>gene-\2|Aste/' your_input.fasta > your_output.fasta


## Orthology inference
conda activate test_env
orthofinder -f proteomes_folder/
#dir1 must be a directory containing all proteoms with an extension of either .fasta or .prot

cd Orthogroups/
#Orthogroups.tsv : For every orthogroup there are proteins associated
#Orthogroups.GeneCount.tsv : Percentage or number of proteins per species

#Align single copy orthogroups
cd Analysis/Tree/Aln/
ln -s Single_Copy_Orthologue_Sequences/
for i in OG000*; do mafft --auto --thread 10 "$i" > "${i/.fa/.mafft}"; done

for i in *.mafft; do trimal -in "$i" -gappyout -out "${i/.mafft/.trimmed.mafft}"; done
#Trimming reduces alignment dimension, reducing computational time

# NCBI alignment viewer: used to see alignment coverage (view variable regions and core regions)
    #View one trimmed and one non-trimmed sequence: After trimming all major drops in coverage are removed

# Concatenating protein fasta
for i in *trimmed.mafft; do sed 's/>.*|/>/' "$i" > "${i/.mafft/.renamed.mafft}"; done
AMAS.py concat -y nexus -i *renamed.mafft -f fasta -d aa
#Creates a concatenated file and a partitioning file in nexus format(.txt)

## Phylogenetic inference
conda activate test_env
cd Tree/iqtree
iqtree -s concatenated.out -spp partitions.txt -bb 1000 -nt AUTO
    # -spp: proportional branches based on the partition file, this allows to estimate longer or shorter branches for each partition
    # -bb 1000: one thousand of bootstrap, although this is ultrafast bootstrap so you need at least 1000 replicates (100 with std bootstrap)

# Partitions.txt.treefile (newick format)
    # itol online (https://itol.embl.de/upload.cgi) to view treefile file

## Divergence Time Estimation
cd Tree/Time_tree
#Find secondary calibration information on timetree (https://timetree.org/) using adjusted_time
#Create time estimation file
nano Dates.txt
    #taxon1,taxon2  -calibration_time (scale is not important as long as it is consistent)
    #Aste,Llon  -241

iqtree -s concatenated.out --date Dates.txt --date-tip 0 -o "<taxon>" -spp partitions.txt.best_scheme.nex -nt 6 --prefix Time.Tree --date-options "-u 1 -k" -te partitions.txt.treefile
    #-u : branches shorter than 1 collapse in a politomy
# Time.Tree.timetree.nex (nexus format) to view on ITOL

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
#Base_clade_results: for each node, how many gene famlies increase and how many decrease
#Base_asr.tre -> ITOL -> Advanced: Node IDs: Display: to see which number corresponds to which node
#Base_family_results.txt: family id - p-value - is there a significant change in the branch?
cut -f3 Base_family_results.txt | sort | unique -c #to see how many families have significant changes

grep "Aste<5>\*" Base_asr.tre | wc -l #to find all significant modifications in Aste branch
varAste=$(grep "<5>\*" Base_asr.tre | cut -d ":" -f8 | cut -d"_" -f2)
varAaeg=$(grep "<1>\*" Base_asr.tre | cut -d ":" -f1 | cut -d"_" -f2)
varAsub=$(grep "<2>\*" Base_asr.tre | cut -d ":" -f2 | cut -d"_" -f2)
varCqui=$(grep "<4>\*" Base_asr.tre | cut -d ":" -f6 | cut -d"_" -f2)
varWsmi=$(grep "<3>\*" Base_asr.tre | cut -d ":" -f4 | cut -d"_" -f2)

var8=$(grep "<5>\*" Base_asr.tre | cut -d":" -f9 | cut -d"_" -f2) #previous node (8) value for significant families
varOG=$(grep "<5>\*" Base_asr.tre | cut -d" " -f4)

paste <(printf %s "$varOG") <(printf %s "$varAste") <(printf %s "$var8") > Aste_Significant.txt
#upload to R and create another column that subtracts the first and the second column: based on if it's negative or positive it is a contraction or expansion
awk '{print $0, $2 - $3}' Aste_Significant.txt > Aste_Significant.txt.diff
sort -k4,4n Aaeg_Significant.txt.diff | grep -v "-" | awk '{print $1}' > Aste_Significant.GO
while read line; do grep "Aste" Orthofinder/orthogroup_Sequences/"$line"*.fa; done < Aste_Significant.GO > Aste_Significant_Genes.txt


#Panzer2: reasearch by homology of GO-Terms (http://ekhidna2.biocenter.helsinki.fi/sanspanz/)
    #upload all Aste proteins
    #GO prediction: Arthropoda

http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//VvjXh3TeoFe/index.html
wget GO.out

#Add ":" before each GO id
awk '{print $1, "GO:"$3}' GO.out > Gene_GO.tsv

#Sort GO id for each Gene
awk '{ a[$1]=a[$1]","$2; } END { for (i in a) {sub(/,/,"",a[i]);printf "%s %s\n",i,a[i] } }' Gene_GO.tsv > Gene_GO_annotation.tsv

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