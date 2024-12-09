#go2genes = go term annotation for all genese
#list_expanded_genes.txt = list of genes of interest

library(tidyverse)
library(topGO)

# the following commented passages are part of the original script from where I developed the function `GOenrichment`

# geneID2GO <- readMappings(file = "go2genes.txt")
# geneUniverse <- names(geneID2GO)

# genesOfInterest <- read.table("list_expanded_genes.txt",header=FALSE)
# genesOfInterest <- as.character(genesOfInterest$V1)
# geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
# names(geneList) <- geneUniverse
# myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5)

# Expanded_BP_classic_fisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")

# Expanded_BP_allRes <- GenTable(myGOdata, Classic_Fisher = Expanded_BP_classic_fisher, orderBy = "Classic_Fisher", topNodes=1000)

# Expanded_BP_allRes$Classic_Fisher <- as.numeric(Expanded_BP_allRes$Classic_Fisher)
# Expanded_BP_allRes <- subset(Expanded_BP_allRes, Classic_Fisher < 0.05)


GOenrichment <- function(trait, trait_name) {
  genesOfInterest <- as.character(trait$V1) #as vector not character 
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  print(trait_name)
  
  ontology_values = c("BP", "MF", "CC")
  
  GOdata_list <- lapply(ontology_values, function(ontology_value) {
    GOdata_name <- paste("GOdata_", ontology_value, sep = "")
    # annot = annFUN.gene2GO this imparts the program which annotation it should use. In this case, it is specified that it will be in gene2GO format and provided by the user.
    # gene2GO = gene_universe is the argument used to tell where is the annotation
    assign(GOdata_name, new("topGOdata", ontology=ontology_value, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene_universe))
  })
  
  elim_list <- lapply(seq_along(ontology_values), function(i) {
    elim_name <- paste("elim_", ontology_values[i], sep="")
    assign(elim_name, runTest(GOdata_list[[i]], algorithm="elim", statistic="fisher"))
  })
  
  results_elim <- function(GO_data, elim_data) {
    resulte <- GenTable(GO_data, Classic_Fisher = elim_data, orderBy = "Classic_Fisher", topNodes=1000, numChar=1000)
    resulte$Classic_Fisher <- as.numeric(resulte$Classic_Fisher)
    resulte <- subset(resulte, Classic_Fisher < 0.05)
    return(resulte)
  }
  
  results_elim_list <- lapply(seq_along(ontology_values), function(i) {
    resulte_name <- paste("resulte_", ontology_values[i], sep="")
    assign(resulte_name, envir = .GlobalEnv, results_elim(GOdata_list[[i]], elim_list[[i]]))
  })
  
  write_elim_results <- function(result, ontology_value, raw_trait_name) {
    trait <- gsub("_", "", gsub("^s", "", raw_trait_name))
    table_name <- paste("02_enrichment/topGOe_", trait, "_", ontology_value, ".txt", sep="")
    write.table(result, file=table_name, quote=F, sep = "\t", row.names = F)
  }
  
  lapply(seq_along(ontology_values), function(i) {
    write_elim_results(results_elim_list[[i]], ontology_values[i], trait_name)
  })
}

#this particular syntax has been necessary since it was impossible to give the function the trait name it was computing.
GO_enrichment <- function(list) {
  lapply(seq_along(list), function(i) {
  GOenrichment(list[[i]], names(list)[i])
  })
}

#Final function to run
GO_enrichment()