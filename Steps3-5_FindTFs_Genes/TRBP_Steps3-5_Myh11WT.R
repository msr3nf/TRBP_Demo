setwd("~/Desktop/") #change to directory of interest

#load necessary packages
library(Biostrings)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
library(biomaRt)
library(org.Mm.eg.db)
library(jsonlite)
library(tidyverse)
library(httr)
library(stringr)
library(data.table)
library(purrr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)
library(dplyr)

#####STEP 3:#####
#read in PPRs for genes in the genome of interest (i.e. mm10)
##to see how this file was generated and to create your own PPR library, refer to TRBP_Steps1-2_Myh11WT.R on our GitHub Repo
WholeGenome_List <- read.csv("WholeGenome_List_mm10_20kPPRseq_03162022.csv") #this is a large file, so we recommend subsetting for a smaller list of genes or using a high-performance computing system (HPCS) 
WholeGenome_List <- WholeGenome_List[,-1]

##select the PPR for your gene of interest...for our demo, we chose Myh11
Myh11_PPR <- subset(WholeGenome_List, subset = mgi_symbol== "Myh11")
write.csv(Myh11_PPR, "Myh11_TF_list_20kBPregion.csv")

#notice that there are 3 exon_id's associated with Myh11 according to GenomicRange
##if there are multiple exon_ids for the same gene, then enter the gene name (i.e. Myh11) in the UCSC genome browser and select the information for the exon_id available on UCSC genome
###because Myh11's PPR is on the negative DNA strand, we must select Chr16:14291365 as the TSS (verified on UCSC genome)

PPR_start <- 14291365-10000
PPR_end <- 14291365+10000
#REST API from UCSC Genome for retrieving TFs assocaited with PPR
responseTF = fromJSON("https://api.genome.ucsc.edu/getData/track?genome=mm10;track=hub_186875_JASPAR2022_TFBS_mm10;chrom=chr16;start=14281365;end=14301365")
Myh11_TF_list <- responseTF$hub_186875_JASPAR2022_TFBS_mm10
colnames(Myh11_TF_list)[4] <- "MatrixID"
colnames(Myh11_TF_list)[7] <- "TF_name"
write.csv(Myh11_TF_list, "Myh11_TF_list_20kBPregion.csv")
#find DNA binding sequences for Myh11-associated TFs
Myh11_TF_list <- read.csv("Myh11_TF_list_20kBPregion.csv") 
Myh11_TF_list <- Myh11_TF_list[,-1]

for (i in 1:nrow(Myh11_TF_list)) {
  chrom <- Myh11_TF_list$chrom[[i]]
  chromStart <- Myh11_TF_list$chromStart[[i]]
  chromEnd <- Myh11_TF_list$chromEnd[[i]]
  TF_name <- Myh11_TF_list$TF_name[[i]]
  
  my.dnastring <- as.character(Biostrings::getSeq(Mmusculus, chrom, chromStart, chromEnd))
  Myh11_TF_list$TF_sequence[[i]] <- my.dnastring
  
  print(paste0(i, "_", TF_name, "_done"))
}

write.csv(Myh11_TF_list, "Myh11_TF_list_20kBPregion_03162022.csv")
Myh11_TF_list <- read.csv("Myh11_TF_list_20kBPregion_03162022.csv")


Myh11_TF_list <- Myh11_TF_list[,-1]
Myh11_TF_list['Mouse_Genes_with_TF_sequence'] <- NA
Myh11_TF_list['NEWMouse_Genes_with_TF_sequence'] <- NA


Myh11_TF_list <- subset(Myh11_TF_list, subset = chromStart > 14281365) #this is to further select for TFs that bind completely within the PPR specified
Myh11_TF_list <- subset(Myh11_TF_list, subset = chromEnd < 14301365 )

#####STEP 4:#####
#Identify all genes with PPRs that can bind Myh11 TFs --> this code does not require parallel processing (i.e. run locally using small gene lists/PPR sizes)
##We do not reccomend this method for whole genome queries
for (j in 1:nrow(Myh11_TF_list)) {
  result <- WholeGenome_List[WholeGenome_List$PPR_sequence %like% Myh11_TF_list$TF_sequence[[j]], ]
  result_mouse_gene <- paste(as.character(unique(result[[18]])), sep="' '", collapse=", ") #col #16 has all gene names
  Myh11_TF_list$NEWMouse_Genes_with_TF_sequence[[j]] <- result_mouse_gene
  
  
  print(paste0(j, "_", Myh11_TF_list$TF_name[[j]], "_done"))
}
write.csv(Myh11_TF_list, "Myh11_TF_list_wMouseGenesONLY_031720222.csv")

##The recommended, faster method that is optimal for HPCS is shown below:

unique_Myh11_TFseq <- as.data.frame(unique(Myh11_TF_list$TF_sequence))
colnames(unique_Myh11_TFseq)[1] <- "TF_sequence"
unique_Myh11_TFseq$Mouse_Genes_with_TF_sequence <- NA

seqMatch <- function(x){
  result <- WholeGenome_List[WholeGenome_List$PPR_sequence %like% x, ]
  result_mouse_gene <- paste(as.character(unique(result[[17]])), sep="' '", collapse=", ") #col 17 has all gene names --> check this prior to running
  return(result_mouse_gene)
}

numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1

options(mc.cores=numCores)

#this only runs for the first 10000 unique TF sequences
##for faster results when working with larger data queries (i.e. 30968 unique TF sequences), we recommend "splitting up" the number of unique TF sequences, as shown below
###to run this, create two new R files containing lines 87-99 in addition to each "tier" --> now you can run all 3 in parallel (x2)
unique_Myh11_TFseq$Mouse_Genes_with_TF_sequence[1:10000] <- unlist(mclapply(unique_Myh11_TFseq$TF_sequence[1:30968], seqMatch))
write.csv(unique_Myh11_TFseq, "unique_Myh11_TFseq_10000_03172022.csv")

#tier 2
unique_Myh11_TFseq$Mouse_Genes_with_TF_sequence[10001:20000] <- unlist(mclapply(unique_Myh11_TFseq$TF_sequence[10001:20000], seqMatch))
write.csv(unique_Myh11_TFseq, "unique_Myh11_TFseq_20000_03172022.csv")

#tier 3
unique_Myh11_TFseq$Mouse_Genes_with_TF_sequence[20001:30968] <- unlist(mclapply(unique_Myh11_TFseq$TF_sequence[20001:30968], seqMatch))
write.csv(unique_Myh11_TFseq, "unique_Myh11_TFseq_30000_03172022.csv")

#import all 3 files which took ~9.5 hrs total to run
tier1 <- read.csv("unique_Myh11_TFseq_10000_03172022.csv")
tier2 <- read.csv("unique_Myh11_TFseq_20000_03172022.csv")
tier3 <- read.csv("unique_Myh11_TFseq_30000_03172022.csv")

allTiers <- rbind(tier1, tier2, tier3) 
allTiers <- allTiers[,-1]

#check that the split TF sequences align --> sum should be 30968
length(intersect(allTiers$TF_sequence, Myh11_TF_list$TF_sequence)) #checks out --> should be 30968 with all files

#match gene lists by Myh11-specific TF seq
i=1
resultsTable <- Myh11_TF_list
resultsTable <- resultsTable[-(1:75532),] #creating empty data results table

for (i in 1:nrow(allTiers)) {
  data <- subset(Myh11_TF_list, subset = TF_sequence == allTiers$TF_sequence[[i]])
  data$Mouse_Genes_with_TF_sequence <- allTiers$Mouse_Genes_with_TF_sequence[[i]]
  resultsTable <- rbind(resultsTable, data)
  print(paste0(i, "_done"))
  
}

length(intersect(resultsTable$TF_sequence, allTiers$TF_sequence)) #all check out!
write.csv(resultsTable, "Myh11_TF_list_wMouseGenesONLY_031820222_FINAL.csv")

#now you have a preliminary results table showing TFs with all potentially co-regulated genes whose PPRs can bind TFs!
