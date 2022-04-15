if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NOISeq")

library(NOISeq)

matrix_dg <- read.table("pernoiseq.txt", sep = " ",header = TRUE, row.names = 1)
View(matrix_dg)
factors <- data.frame(dg_conditions = c("Aste","Aste","Aste","Hill","Hill","Hill"))

dg_noiseq <- readData(data=matrix_dg, factors = factors)
dg_noiseq
head(assayData(dg_noiseq)$exprs)
#dg_noiseq <- addData(dg_noiseq, length = mylength)

# look at data quality

mysaturation = dat(dg_noiseq, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)
mycountsbio = dat(dg_noiseq, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")

#filter loci with low counts

myfilt10 = filtered.data(matrix_dg, factor = factors$dg_conditions, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")
dg_noiseq <- readData(data=myfilt10, factors = factors)
mycountsbio = dat(dg_noiseq, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
dg_TMM10 = tmm(assayData(dg_noiseq)$exprs, long = 1000, lc = 0)
#lc = 0 and long = 1000 -->no length correction is applied; if we use FPKM --> lc = 1
write.table(dg_TMM10, file="dg_tmm10", append = FALSE, eol="\n", quote = FALSE)


#Differential expression between two conditions

dg_noiseq_aste_hill <- readData(data=dg_TMM10, factors = factors)
mynoiseqbio_aste_hill=noiseqbio(dg_noiseq_aste_hill, k=0.1, norm="n", filter=0, factor="dg_conditions")
#k is used to replace the zero values in the expression matrix with other non-zero value in order to avoid indetermination in some calculations such as fold-change
mynoiseqbio_aste_hill_deg = degenes(mynoiseqbio_aste_hill, q = 0.95, M = NULL)
#q=1-FDR
#M if = "up" --> show up-regulated in condition 1; if = "down" --> show down-regulated in condition 1, if = NULL --> show all differentially expressed features
write.table(mynoiseqbio_aste_hill_deg, file="aste_hill_deg", append = FALSE, eol="\n", quote = FALSE)


#plot the average expression value and highlight the feature differentially expressed
DE.plot(mynoiseqbio_aste_hill, q = 0.95, graphic = "expr", log.scale = TRUE)
dev.copy2pdf(file= "Aste_Hill_expr.pdf")

#plot the log-FC , M=log2FC, D= |exprCond1 - exprCond2|
DE.plot(mynoiseqbio_aste_hill, q = 0.95, graphic = "MD")
dev.copy2pdf(file= "Aste_Hill_DM.pdf")