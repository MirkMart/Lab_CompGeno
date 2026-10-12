# Laboratorio di Genomica Comparata (Comparative Genomics Laboratory) AA 2026-2027

This repository collects all the information and materials we developed to teach and support the "Laboratorio di Genomica Comparata" course at the University of Bologna for the 2026-2027 academic year. A brief excursus about what is present in this folder is:

- [00_Bash](./00_Bash/): Here are listed ALL the commands that we will use, their explanations, and some examples of their use. This is meant to be a compendium of how to talk to our computer, and nobody expects you to know it by heart.
- [01_ComputerEnv_ProjectWorkflow](./01_ComputerEnv_ProjectWorkflow/): This folder includes how the environments are named and organised on our machine and some information about our study system, _Anopheles stephensii_.
- [02_KmerBased_GenomeSurvey](./02_KmerBased_GenomeSurvey/): _a good start is half of the battle_. Good-quality input is essential to get good-quality results. Here we present the FASTQ format and how to assess sequencer output quality.
- [03_GenomeAssembly](./03_GenomeAssembly/): This folder includes the descriptions of how to assemble our genome. We have contig-level assembly, genome polishing, and genome evaluation to know how well we performed the first three passages. Finally, we perform contaminant detection and reference-based scaffolding.
- [04_GenomeAnnotation](./04_GenomeAnnotation/): Here we perform the genome annotation of our brand new _Anopheles stephensii_ genome, discovering its features.
- [05_OG.Inference_Phylogenomic](./05_OG.Inference_Phylogenomic/): Here we infer orthologs and then use the resulting orthgroups to infer our species tree.
- [06_DivergenceTime_Estimation](./06_DivergenceTime_Estimation/): Using all the information we gathered so far, we estimate the divergence time between our species of interest.
- [07_GeneFamilies_Evolution](./07_GeneFamilies_Evolution/): genes evolve in several ways. Gene duplication and deletion—expanding or contracting gene families—is only one way. Using CAFE, we can explore this.
- [08_GeneAnnotation](./08_GeneAnnotation/): This folder includes explanations and tips on how to annotate genes in our genome (personally my favorite part). This functional part may be broad, but it is meant to be so to let everyone's creativity thrive. We will use GO terms, codeml to search for signs of selection in our genes of interest, and the KEGG database to understand how everything connects.

- [99_Figures](./99_Figures/): The folder contains some figures that are used in the README files across this course.
