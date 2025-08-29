library("seqinr")

#Create an empty vector for sequence names and read the kept.msk file.

sequence <- character() #this is the character vector containing all sequences names
kept.msk <- scan(file=snakemake@input$mask,what="character",quiet=TRUE)

sequence <- snakemake@wildcards$orthogroup 
# #add the name of the elaborated sequence

#Read conserved sites from kept.msk
current.first.kept <- 2 #always third position in kept mask since it is only one mask
current.last.kept <- length(kept.msk)
current.kept <- as.numeric(kept.msk[current.first.kept:current.last.kept]) #vector of positions on which the script works

#Back-translate aminoacids into nucleotides.

current.aminoacids <- read.fasta(file=snakemake@input$amino,seqtype="AA",forceDNAtolower=FALSE,strip.desc=TRUE)
current.nucleotides <- read.fasta(file=snakemake@input$nucle,seqtype="DNA",forceDNAtolower=FALSE,strip.desc=TRUE)

if (!all(getName(current.aminoacids) == getName(current.nucleotides))) {
	message(paste(sequence,": warning! Different names were detected between aminoacid and nucleotide file!",sep=""))
}

num.current.taxa <- length(current.aminoacids)
for (t in 1:num.current.taxa) {
	current.back.translation <- character()
	current.masked.back.translation <- character()
	gaps.aminoacids <- 0
	for (s in 1:length(current.aminoacids[[t]])) {
		if (current.aminoacids[[t]][s] == "-") {
			gaps.aminoacids <- gaps.aminoacids+1
			current.back.translation <- c(current.back.translation,"-","-","-")
			}
	  else if (any(is.na(current.nucleotides[[t]][(s-gaps.aminoacids)*3-2]),is.na(current.nucleotides[[t]][(s-gaps.aminoacids)*3-1]),is.na(current.nucleotides[[t]][(s-gaps.aminoacids)*3]))) {
			message(paste(sequence,": warning! NA to be inserted!",sep=""))
			current.back.translation <- c(current.back.translation,current.nucleotides[[t]][(s-gaps.aminoacids)*3-2],current.nucleotides[[t]][(s-gaps.aminoacids)*3-1],current.nucleotides[[t]][(s-gaps.aminoacids)*3])
			}
		else current.back.translation <- c(current.back.translation,current.nucleotides[[t]][(s-gaps.aminoacids)*3-2],current.nucleotides[[t]][(s-gaps.aminoacids)*3-1],current.nucleotides[[t]][(s-gaps.aminoacids)*3])

#Write a codon into masked back-translated files if it is a conserved site

		if (s %in% current.kept) {
			if (current.aminoacids[[t]][s] == "-") {
				current.masked.back.translation <- c(current.masked.back.translation,"-","-","-")
				}
			else if (any(is.na(current.nucleotides[[t]][(s-gaps.aminoacids)*3-2]),is.na(current.nucleotides[[t]][(s-gaps.aminoacids)*3-1]),is.na(current.nucleotides[[t]][(s-gaps.aminoacids)*3]))) {
				message(paste(sequence,": warning! NA to be inserted!",sep=""))
				current.masked.back.translation <- c(current.masked.back.translation,current.nucleotides[[t]][(s-gaps.aminoacids)*3-2],current.nucleotides[[t]][(s-gaps.aminoacids)*3-1],current.nucleotides[[t]][(s-gaps.aminoacids)*3])
				}
			else current.masked.back.translation <- c(current.masked.back.translation,current.nucleotides[[t]][(s-gaps.aminoacids)*3-2],current.nucleotides[[t]][(s-gaps.aminoacids)*3-1],current.nucleotides[[t]][(s-gaps.aminoacids)*3])
			}
		}
	write.fasta(current.back.translation,getName(current.aminoacids[t]),snakemake@output$retro,open="a",nbchar=length(current.back.translation))
	write.fasta(current.masked.back.translation,getName(current.aminoacids[t]),snakemake@output$retro_mask,open="a",nbchar=length(current.masked.back.translation))
	}
message(paste("Orthogroup",sequence,"back-translated!",sep=" "))
	
