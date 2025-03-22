#This code shows users how we generated PPRs for our input gene and gene list of interest
##Our input gene is the murine Myh11
###Our gene list is the entire murine genome (version mm10)

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
library(purrr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)
library(dplyr)

#####STEP 1:#####
##We have already defined our gene of interest as Myh11, so we will now use the Ensembl database to gather information for our gene list of interest, which consists of all genes in the murine genome (including our input gene, Myh11)
{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("mmusculus_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  anno <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "entrezgene_id","mgi_symbol", "description", "external_gene_name"), mart=ensembl)
}

view(anno)
mouseGenome <- as.data.frame(anno)
mouseGenome <- subset(mouseGenome, subset = entrezgene_id != "NA") #we need valid entrezgene ids to proceed with GenomicRanges to charactreize PPRs, so select for genes that have documented entrezgene id's 


##follow vignette from GenomicRanges for deriving genes' 5' UTR: http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all_txdb_ids <- keys(TxDb.Mmusculus.UCSC.mm10.knownGene) 
all_txdb_ids_table <- as.data.frame(all_txdb_ids)
colnames(all_txdb_ids_table)[1] <- "list"
check_list <- intersect(all_txdb_ids_table$list,mouseGenome$entrezgene_id) #This is the final list of 220225 genes I will use for the mouse genome
check_list <- as.character(check_list)
txbygene <- transcriptsBy(txdb, "gene")[check_list] 

map <- relist(unlist(txbygene, use.names=FALSE)$tx_id, txbygene)
map
map1 <- as.data.frame(map)
colnames(map1)[2] <- "entrezgene_id"


fiveUTR <- fiveUTRsByTranscript(txdb)
txid <- unlist(map, use.names=FALSE)
fiveUTR <- fiveUTR[names(fiveUTR) %in% txid]

length(fiveUTR)
fiveUTR
UTR_List <- as.data.frame(fiveUTR) 
colnames(UTR_List)[2] <- "value"

UTR_List1 <- merge(UTR_List, map1, by = "value")
UTR_List2 <- merge(UTR_List1, mouseGenome, by = "entrezgene_id")

All_genes_5UTR <- UTR_List2[UTR_List2[,"exon_rank"] == 1,]
write.csv(All_genes_5UTR, "WholeGenome_mm10_5UTR_Coord.csv")

WholeGenome_List <- read.csv("WholeGenome_mm10_5UTR_Coord.csv")
WholeGenome_List <- WholeGenome_List[,-1]
WholeGenome_List['PPR_sequence'] <- NA
length(unique(WholeGenome_List$entrezgene_id)) #19325 total number genes; note there are repeats for different gene versions

#####STEP 2:#####
#Note: If determining PPR information for large list of genes, such as an organism's entire genome, we do NOT recommend running the following loop locally
##Therefore, we strongly suggest running this code on a high-performance computing system 
###When determining the nucleotide strings for PPRs, this loop accounts for whether the TSS is on the negative or positive strand of DNA

for (i in 1:nrow(WholeGenome_List)) {
  if (WholeGenome_List$strand[[i]] == "+") {
    chrom <- WholeGenome_List$seqnames[[i]]
    chromStart <- WholeGenome_List$start[[i]]-10000 #this is where users can modify how many nucleotides upstream...
    chromEnd <- WholeGenome_List$start[[i]]+10000 #... or downstream from the TSS comprise the PPR 
    Gene_name <- WholeGenome_List$external_gene_name[[i]]
    my.dnastring <- as.character(Biostrings::getSeq(Mmusculus, chrom, chromStart, chromEnd))
    
  }
  
  if (WholeGenome_List$strand[[i]] == "-") {
    chrom <- WholeGenome_List$seqnames[[i]]
    chromStart <- WholeGenome_List$end[[i]]-10000
    chromEnd <- WholeGenome_List$end[[i]]+10000
    Gene_name <- WholeGenome_List$external_gene_name[[i]]
    my.dnastring <- as.character(Biostrings::getSeq(Mmusculus, chrom, chromStart, chromEnd))
    
  }
  print(paste0(i, "_", Gene_name, "_done"))
  
  WholeGenome_List$PPR_sequence[[i]] <- my.dnastring
  
}
length(unique(WholeGenome_List$external_gene_name)) #19321 unique gene IDs

WholeGenome_List1 <- apply(WholeGenome_List,2,as.character)
write.csv(WholeGenome_List1, "WholeGenome_List_mm10_20kPPRseq_03162022.csv")

##to conserve file size, users can remove extraneous information using the following sample code:
WholeGenome_List <- WholeGenome_List[!duplicated(WholeGenome_List$PPR_sequence), ]
WholeGenome_List <- WholeGenome_List[,-1]
WholeGenome_List <- WholeGenome_List[,-c(2,3,7,10)]
WholeGenome_List <- WholeGenome_List[,-c(8)]
write.csv(WholeGenome_List1, "WholeGenome_List_mm10_20kPPRseq_03222025.trimmed.csv")


