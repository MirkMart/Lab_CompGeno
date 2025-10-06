# Genome annotation

Once you have a de-novo assembled genome, what you usually want to know is what genomic features are on that genome. For this purpose we need an annotation pipeline. But what can we annotate on our genome?

1. Repetitive elements (Transposons, low complexity regions, Microsatellistes)
2. Protein - coding Genes
3. Pseudo - genes
4. Annotation of tRNA, ncRNA, ecc…
5. UTR regions

Eukariotic genomes are very complex (results of evolution by molecular tinkering) and gene structures can greatly vary between lineages and even between closely relatetd species (exon length, intron length, gene length, splicing site sequence patterns, ecc…). In model species exist highly curated models able to describe each feature (with the form of HMM profiles). This models can be directly used for *ab - initio* gene prediction. However for non model species most of the time these patterns are unknown. To partially overcome this issue, we can use **external evidences** to directly annotates some features on the genome (*evidence - based* gene prediction) and/or as a source for *ab - initio* gene predictors training.

-----

## Genome annotation (easy) workflow

1. Annotation of repetitive elements in the genome. Repetitive regions and transposons must be excluded from gene annotation because they have completly different structures.
2. Alignment of high quality external evidences (Transcriptomes and proteins).
3. Extraction of an initial set of highly confident gene models.
4. Training of *ab-initio* gene predictors softwares.
5. Re - annotation of the genome through trained models. **Steps 3, 4 and 5 are usually repeated multiple times to improve gene prediction at each iteration**
6. Collection of all gene models and creation of a consensus set based on support given by external evidences (usually gene models low or not supported by external evidences are discared).
7. Evaluation of our final gene set.

Luckly, multiple genome annotation pipeline exist to facilitate this process. One of the most used is **[MAKER](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018)**. Another [great step-wise tutorial](./maker_pipeline.pdf) is loaded in this folder.

![MAKER](../99_Figures/MAKER_Apollo_view.jpeg)

-----

## 1\. Transposons annotation

For a quick and dirty transposons annotation we usually used **[RepeatModeler](https://www.repeatmasker.org/RepeatModeler/)** and **[RepeatMasker](http://www.repeatmasker.org/)**. Basically, RepeatModeler will collect a set of consensus sequences rapresentative of the repeats present in our genome. This library can than be used to annotate the genome. However, take in mind that this approach is usually good only for a raw transposons annotation but is not enough if the goal of our study is an indeep analyses of TEs.

Due to computationally limits I have alredy run RepeatModeler for you, the path of the consensus library is:

> /home/PERSONALE/mirko.martini3/01_2024/00_Data/02_annotation/Anoste_RepeatModeler_library.fa

RepeatMasker will be run inside Maker and we will summarize results at the end.

## 2.1 MAKER rnd 1 (~24.5 hours with 40 CPUs)

First of all you need to collect yout external evidences, for example on NCBI or Uniprot. Unfortunately we don’t have a transcriptome but in real scenarios is almost mandatory. Then you need to create and compile the maker configuration file necessary for the run.

```bash
maker -CTL
```

Three files are generated:

- maker_bopts: parameters for the softwares used by MAKER, we mantain default parameters. Here, a great source to understand each [input parameter](https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained) in MAKER.
- maker_exe: paths to all softwares
- maker_opts: MAKER run parameters, here we need to modify

Now we need to set the paths of the input files in the `opts.ctl` file.

```ctl
    #-----Genome (these are always required)
    genome= <GENOME> #genome sequence (fasta file or fasta embeded in GFF3 file)
    organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic
    
    #-----Re-annotation Using MAKER Derived GFF3
    maker_gff= #MAKER derived GFF3 file
    est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
    altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
    protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
    rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
    model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
    pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
    other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no
    
    #-----EST Evidence (for best results provide a file for at least one)
    est= #set of ESTs or assembled mRNA-seq in fasta format
    altest= #EST/cDNA sequence file in fasta format from an alternate organism
    est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
    altest_gff= #aligned ESTs from a closly relate species in GFF3 format
    
    #-----Protein Homology Evidence (for best results provide a file for at least one)
    protein= <PROTEOME> #protein sequence file in fasta format (i.e. from mutiple oransisms)
    protein_gff=  #aligned protein homology evidence from an external GFF3 file
    
    #-----Repeat Masking (leave values blank to skip repeat masking)
    model_org=<EMPTY> #select a model organism for RepBase masking in RepeatMasker
    rmlib= <RepeatModeler_library> #provide an organism specific repeat library in fasta format for RepeatMasker
    repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
    rm_gff= #pre-identified repeat elements from an external GFF3 file
    prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
    softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)
    
    #-----Gene Prediction
    snaphmm= #SNAP HMM file
    gmhmm= #GeneMark HMM file
    augustus_species= #Augustus gene prediction species model
    fgenesh_par_file= #FGENESH parameter file
    pred_gff= #ab-initio predictions from an external GFF3 file
    model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
    est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
    protein2genome=<0 OR 1> #infer predictions from protein homology, 1 = yes, 0 = no
    trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
    snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
    unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
    
    #-----Other Annotation Feature Types (features MAKER doesn't recognize)
    other_gff= #extra features to pass-through to final MAKER generated GFF3 file
    
    #-----External Application Behavior Options
    alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
    cpus=<CPUS> #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)
    
    #-----MAKER Behavior Options
    max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
    min_contig=1 #skip genome contigs below this length (under 10kb are often useless)
    
    pred_flank=200 #flank for extending evidence clusters sent to gene predictors
    pred_stats=<0 or 1> #report AED and QI statistics for all predictions as well as models
    AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
    min_protein=50 #require at least this many amino acids in predicted proteins
    alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
    always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
    map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
    keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
    
    split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
    single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
    single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
    correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
    
    tries=2 #number of times to try a contig if there is a failure for some reason
    clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
    clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
    TMP= #specify a directory other than the system default temporary directory for temporary files
```

Changes to do:

- add genome.
- add proteome (external evidence. Multiple ones can be listed using commas).
- remove model_org.
- provide RepeatModeler library.
- protein2genome= from 0 to 1 #This change enable maker to use protein homology to infer prediction (needed since we are feeding the program with a known proteome)
- decide number of CPUs to use
- phred_stats from 0 to 1 #This enables the program to calculate phred score for each annotation.
- add minimum length of protein allowed to be added to the annotation (min_protein=50)
- alternative_splice from 0 to 1
- split_hit empty

External evidence (proteomes):
`/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/02_Annotation`

**To be known**: AED is a measure of annotation quality, indicating how well the predicted gene matches external evidence. It ranges from 0 (perfect match) to 1 (no match).

Then just type:

```bash
maker -base <OUTPUT PREFIX>
```

To create a [gff3](http://www.ensembl.org/info/website/upload/gff3.html) and a transcript/protein file we can use:

```bash
fasta_merge -d <DATASTORE INDEX FILE>
gff3_merge -d <DATASTORE INDEX FILE>
```

Finally we must extract all evidences from the gff file using bash (awk on "protein2genome" and "repeatmasker" strings). Try to do it by your self.

We can also summarize results of Evidence-based gene annotation and RepeatMasker using the [AGAT](https://github.com/NBISweden/AGAT) package, a very usefull set of perl scripts to manage gff3 files, print the help and run the script :

```bash
agat_sp_statistics.pl --gff file.gff -o <output_file>        #Summary statistics of gene models
agat_sq_repeats_analyzer.pl -i <input_file> -o <output_file> #Summary statistics of repeats
```

Now try to summarize the results (in R) to see the genome proportion covered by each TE class and compare results with other species and assemblies.

We can also assess the completeness of our first proteome using BUSCO. In this way, we can study how proteome quality changes along our pipeline.

## 2.2 SNAP and Augustus training

**NB:** take in mind that without transcripts evidences gene predictors training won’t work very well, but the process would be the same.

### SNAP

SNAP is a gene prediction tool designed to identify protein-coding genes in genomic DNA sequences, particularly in newly sequenced genomes. SNAP uses a probabilistic approach based on HMMs to model various features of genes and can be trained using a set of known gene models to improve its accuracy for a specific organism. SNAP is often used alongside other gene prediction tools and results from multiple tools can be combined to create a consensus annotation.

Fathom is part of te SNAP package and it is used to prepare, train, and elaborate our gene prediction.

```bash
#I suggest you to create a small script with all these commands because are really fast
maker2zff -c 0 -e 0  -l 80 -x 0.1 -d <datastore_index> #To extract gene models based on mutiple filter criterion. It transforms maker output into .zff files, the format for SNAP
fathom <.ANN FILE> <.DNA FILE> -gene-stats #Print some summary statistics of the selected gene models
fathom <.ANN FILE> <.DNA FILE> -validate #Validate gene models and print summary statistics. output is in the stdout and cannot be redirected
fathom <.ANN FILE> <.DNA FILE> -categorize 1000 #Extract gene modeles together with 1000 bp at both ends for training
fathom <UNI.ANN FILE> <UNI.DNA FILE> -export 1000 -plus #Export and convert unique genes to strand plus. Tipically, only the unique genes are sued for training and testing. The value limits the intergenic sequence at the ends.
forge <export.ann> <export.dna> #call the parameter estimation program, better in another directory
hmm-assembler.pl <NAME> <FORGE DIRECTORY> > <OUTPUT HMM FILE>
```

With:

- -c #The fraction of splice sites confirmed by an EST alignment, default 0.5 (50% of the splice sites should be supported by EST data).
- -e #The fraction of exons that overlap an EST alignment, default 0.5 (50% of the exons in a predicted gene must overlap with an EST for the gene to be kept).
- -l #The min length of the protein sequence produced by the mRNA. It filters out very short, potentially false predictions (e.g., fragments or incomplete genes).
- -x #Max AED to allow 0.5 is default

### Augustus

Augustus training is more complex and computationally intensive, but is one of the most well - performing predictors. For a full augusut training some good options are:

1. Select only well supported gene models.
2. Remove incomplete and short gene models.
3. Remove gene models spaced by lass than 1 - 2 kb.
4. Remove gene models with similar (80% similarity) protein products.

Then you will usually separate your gene models in a training and a testing set. You will use the training set to train augustus and the testing set to test the quality of the models. Finally, the trained models can be used. Again AGAT can help us in this process.

However, in this course we will use a more straigthforward way, training Augustus inside Busco.

```bash
busco -i ../../03_GenomeAssembly/03_scaffolding/Anoste_chr.fasta -c 30 -l /usr/local/share/busco_databases/culicidae_odb12 --augustus --long -m genome --out Anoste_cu --augustus_parameters='--progress=true'
```

Since also this process is quite computanionally intensive and the trained models must be present in a specific path were Augustus will search for them, I have already trained for you augusutus. The name of the models is **Aste** (just use Aste when it will be required).

## 2.3 MAKER rnd 2 (~12h with 40 CPUs)

Now you must change the maker config file specifing (copy of the old one):

1. Add the path of repeats and proteins alignment files (protein_gff and rm_gff).
2. Remove the external protein inputs.
3. Add the path of the SNAP models.
4. The name of the Augustus models -> **Anoste_cu**.
5. Deactive `protein2genome` and `est2genome`.
6. In the case you want to perform a third round, remember to maintain pred_stats from 0 to 1.

Finally, run again MAKER and collect the gff and protein/transcript files.

-----

## Evaluation of a gene set

To evaluate our genome annotation we have mutiple options :

1. Compare gene statistics with bibliography knowledge (*e.g* Mean gene length, mean exon length, mean number of exon and introns per gene, ecc…)
2. Busco on predicted proteins.
3. Summarize AED values.
4. Align the de-novo assembled transctiptome to the genome-based predicted one.

To print sumamry statistics:

```bash
agat_sp_statistics.pl #script for summary statistic, it will take some minutes to run
```

For Busco you should already know how to run it.

`/home/PERSONALE/mirko.martini3/01_2024/03_scripts/AED_cdf_generator.pl` : script for summarize cumulative distribution of AED values.

`/home/PERSONALE/mirko.martini3/01_2024/03_scripts/quality_filter.perl` : script for filter gene models based on AED values.

-----

**NOTES :**

1. Remeber that if you want, you can re - train SNAP based on the new gene models, just redo part 2.2. Then use newly generated HMM models for a new MAKER run. Finally use busco and calculate summary gene stats to see if you get better results.

2. Our annotation pipeline was quite simple. In real case scenario you will use more sophisticated pipelines. Beside performing a complete Augustus training, you should use transcripts evidences coming from the same specie, you can use more proteins (more evidences are always better), you can use more *ab - inito* gene predictors (*e.g* [GeneMark](http://exon.gatech.edu/GeneMark/)) and you can play with MAKER Behavior Options (*e.g* min protein length; split hit length).

3. More annotation pipeline exists, like **[BRAKER](https://github.com/Gaius-Augustus/BRAKER)** as well as softwares to create a consensus set of gene models based on mutiple gene predictors and transcript/protein alignments (*e.g* [EvidenceModeler](https://evidencemodeler.github.io/)).
