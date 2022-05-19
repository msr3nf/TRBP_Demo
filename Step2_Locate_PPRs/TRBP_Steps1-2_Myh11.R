#This code shows users how we generated PPRs for our input gene and gene list of interest
##Our input gene is the murine Myh11
###Our gene list is the entire murine genome (version mm10)

setwd("~/Desktop/") #change to directory of interest

#load necessary packages
library(jsonlite)
library(tidyverse)
library(httr)
library(stringr)
library(data.table)
library(purrr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)
library(dplyr)

#

#follow vignette from GenomicRanges for deriving genes' 5' UTR

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all_txdb_ids <- keys(TxDb.Mmusculus.UCSC.mm10.knownGene) #THIS IS THE FINAL LIST THAT I CAN USE
all_txdb_ids_table <- as.data.frame(all_txdb_ids)
colnames(all_txdb_ids_table)[1] <- "list"
check_list <- intersect(all_txdb_ids_table$list,mouseGenome$entrezgene_id) #THIS IS THE FINAL LIST of 220225 genes I will use
check_list <- as.character(check_list)
txbygene <- transcriptsBy(txdb, "gene")[check_list] #we have 

map <- relist(unlist(txbygene, use.names=FALSE)$tx_id, txbygene)
map
map1 <- as.data.frame(map)
colnames(map1)[2] <- "entrezgene_id"


fiveUTR <- fiveUTRsByTranscript(txdb)
txid <- unlist(map, use.names=FALSE)
fiveUTR <- fiveUTR[names(fiveUTR) %in% txid]

length(fiveUTR)
fiveUTR
UTR_List <- as.data.frame(fiveUTR) #must break up files to get 1000 lines at a time for UCSC genome
colnames(UTR_List)[2] <- "value"

UTR_List1 <- merge(UTR_List, map1, by = "value")
UTR_List2 <- merge(UTR_List1, mouseGenome, by = "entrezgene_id")

All_genes_5UTR <- UTR_List2[UTR_List2[,"exon_rank"] == 1,]
write.csv(All_genes_5UTR, "WholeGenome_mm10_5UTR_Coord.csv")

WholeGenome_List <- read.csv("WholeGenome_mm10_5UTR_Coord.csv")
WholeGenome_List <- WholeGenome_List[,-1]
WholeGenome_List['PPR_sequence'] <- NA
length(unique(WholeGenome_List$entrezgene_id)) #19325 total number genes
