setwd("~/Desktop/") #may have to delete this?
library(jsonlite)
library(tidyverse)
library(httr)
library(stringr)
library(data.table)
library(purrr)
library(GenomicFeatures)
library(dplyr)

Myh11_TFs_Mut_CArG3 <- read.csv("Myh11_TFs_Mut_CArG3.csv")
WholeGenome_List <- read.csv("WholeGenome_List_mm10_20kPPRseq_03162022.csv")

unique_Myh11_TFs_Mut_CArG3_TFseq <- as.data.frame(unique(Myh11_TFs_Mut_CArG3$TF_sequence))
colnames(unique_Myh11_TFs_Mut_CArG3_TFseq)[1] <- "TF_sequence"
unique_Myh11_TFs_Mut_CArG3_TFseq$Mouse_Genes_with_TF_sequence <- NA

seqMatch <- function(x){
  result <- WholeGenome_List[WholeGenome_List$PPR_sequence %like% x, ]
  result_mouse_gene <- paste(as.character(unique(result[[18]])), sep="' '", collapse=", ") #col #16 has all gene names
  return(result_mouse_gene)
}

options(mc.cores=5)

unique_Myh11_TFs_Mut_CArG3_TFseq$Mouse_Genes_with_TF_sequence[1:30928] <- unlist(mclapply(unique_Myh11_TFs_Mut_CArG3_TFseq$TF_sequence[1:30928], seqMatch))
write.csv(unique_Myh11_TFs_Mut_CArG3_TFseq, "unique_Myh11_TFs_Mut_CArG3_TFseq_05032022.csv")


