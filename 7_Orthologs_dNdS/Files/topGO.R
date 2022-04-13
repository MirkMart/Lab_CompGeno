
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")

#older R versions
#source("http://bioconductor.org/biocLite.R")
#biocLite()

#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")

library(topGO)

geneID2GO=readMappings(file="gene_universe.txt")
#The object returned by readMappings is a named list of character vectors. The list names give the genes identifiers. Each element of the list is a character vector and contains the GO identifiers annotated to the specific gene.

str(head(geneID2GO))
geneNames <- names(geneID2GO)
head(geneNames)
int_genes= read.table("geneinterest.txt", sep = "\t")
gene_int_list <- as.vector(int_genes$V1)
str(gene_int_list)
geneList <- factor(as.integer(geneNames %in% gene_int_list))
str(geneList)
names(geneList) <- geneNames
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
#ontology: character string specifying the ontology of interest (BP, MF or CC)
#annFUN.gene2GO this function is used when the annotations are provided as a gene-to-GOs mapping
#nodeSize: an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated genes


resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#runTest is used to apply the specified test statistic and method to the data
#classic algorithm = each GO category is tested independently. 
#elim algorithm = removes the genes mapped to significant GO terms from more general (higher level) GO terms.
#weigth algorithm = genes annotated to a GO term receive weights based on the scores of neighboring GO terms.

resultFis
allRes <- GenTable(GOdata, classicFisher = resultFis,ranksOf = "classicFisher", topNodes = 100)

write.table(allRes,file="topGO_BP.txt", quote=FALSE, row.names=FALSE, sep = "\t")