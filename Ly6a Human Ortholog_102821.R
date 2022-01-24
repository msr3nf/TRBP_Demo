##read in files 
setwd("~/Desktop/R Pdgfrb:Dual")
wholeGenomePPR <- read.csv("wholeMM10genome_PPRseq_no_repPPRseq.csv")
##according to table, Ly6a PPR coord= Chr15:74998348-74996848 --> neg strand, so switch start/end chrom
#get all known TF names + exact binding sequences for this Ly6a PPR region:
library(jsonlite)
library(tidyverse)
library(httr)
library(stringr)
library(data.table)
library(purrr)


responseTF = fromJSON("https://api.genome.ucsc.edu/getData/track?genome=mm10;track=hub_186875_JASPAR2022_TFBS_mm10;chrom=chr15;start=74996848;end=74998348")
Ly6a_TF_list <- responseTF$hub_186875_JASPAR2022_TFBS_mm10
colnames(Ly6a_TF_list)[4] <- "TF_name"
write.csv(Ly6a_TF_list, "Ly6a_TF_list.csv")

##script for Rivanna
#convert chrom. TF locations to binding sequences
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

Ly6a_TF_list['TF_sequence'] <- NA
#if Ly6a PPR is ONLY on - strand, then only consider TFs on - strand?
Ly6a_TF_list <- subset(Ly6a_TF_list, subset=strand=="-")

for (i in 1:nrow(Ly6a_TF_list)) {
    chrom <- Ly6a_TF_list$chrom[[i]]
    chromStart <- Ly6a_TF_list$chromStart[[i]]
    chromEnd <- Ly6a_TF_list$chromEnd[[i]]
    TF_name <- Ly6a_TF_list$TF_name[[i]]
    
    my.dnastring <- as.character(Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, chrom, chromStart, chromEnd))
    Ly6a_TF_list$TF_sequence[[i]] <- my.dnastring
  
    print(paste0(i, "_", TF_name, "_done"))
}

write.csv(Ly6a_TF_list, "Ly6a_TF_list_negStrand.csv")

for (j in 1:nrow(Ly6a_TF_list)) {
  result <- wholeGenomePPR[wholeGenomePPR$PPR.Sequence %like% Ly6a_TF_list$TF_sequence[[j]], ]
  result_mouse_gene <- paste(as.character(result[[16]]), sep="' '", collapse=", ") #col #16 has all gene names
  Ly6a_TF_list$Mouse_Genes_with_TF_sequence[[j]] <- result_mouse_gene
  #eliminate repeats in this results
  humanGenes <- as.data.frame(convertMouseGeneList(result[[16]])) #format result as table entry
  Ly6a_TF_list$Human_Genes_with_TF_sequence[[j]] <-paste(as.character(humanGenes[[1]]), sep="' '", collapse=", ") #this also eliminates repeats! for human set
  #convert "result" to human ortholog  
  print(paste0(j, "_", Ly6a_TF_list$TF_name[[j]], "_done"))
}

#read in Rivanna pt. 1 resutls
ly6a_Human_Mouse_result <- read.csv("Ly6a_TF_list_human_mouse_wholeGenome_geneList.csv")

#read in ATACseq control data as nt sequences
setwd("~/Desktop/R Pdgfrb:Dual") #GSM1876023_020805.2_ATAC_PDGFBB_peaks.bed
controlHCASMC_peakSeq <- read.csv("controlHCASMC_ATACseq.csv")
PDGF_HCASMC <- read.table("GSM1876027_CA1508_ATAC_PDGFBB_peaks.txt")

controlHCASMC_peakSeq['peakGene'] <- NA
#convert peaks to genes based on chr coordinates

for (i in 1:nrow(controlHCASMC_peakSeq)) {
  chrom <- controlHCASMC_peakSeq$chrNum[[i]]
  chromStart <- controlHCASMC_peakSeq$chrStart[[i]]
  chromEnd <- controlHCASMC_peakSeq$chrEnd[[i]]
  peakNum <- controlHCASMC_peakSeq$peakNum[[i]]
  
  responseGene = fromJSON(paste("https://api.genome.ucsc.edu/getData/track?genome=hg19;track=ncbiRefSeq;chrom=",chrom,";start=",chromStart,";end=",chromEnd,sep = ""))
  my.gene <- responseGene$ncbiRefSeq
  my.geneName <- tail(names(sort(table(my.gene$name2))), 1)
  
  if (is.null(my.geneName) == FALSE) { #fun fact i do NOT need the != part
    controlHCASMC_peakSeq$peakGene[[i]] <- tail(names(sort(table(my.gene$name2))), 1)
    print(paste0(i, "_", peakNum, "_done"))
  }

}

##
#determine which of these control peaks contain binding sequences for Ly6a-specific TFs
ly6a_Human_Mouse_result['ATAC_peaks'] <- NA

j=2
for (j in 1:nrow(ly6a_Human_Mouse_result)) {
  result <- controlHCASMC_peakSeq[controlHCASMC_peakSeq$peakSequence %like% ly6a_Human_Mouse_result$TF_sequence[[j]], ]
  result_peak <- paste(as.character(result[[5]]), sep="' '", collapse=", ")
  
  ly6a_Human_Mouse_result$ATAC_peaks[[j]] <- result_peak
  
  print(paste0(j, "_", Lgals3_TF_list$TF_name[[j]], "_done"))
}
setwd("~/Desktop/R Pdgfrb:Dual")
test <- read.csv("controlHCASMC_peakSeq_first25000.csv")




#attempting to parallelize code: DON'T NEED THIS IF YOU JUST USE ONE DATABASE
library(parallel)
library(foreach)
install.packages("doParallel",dependencies=TRUE)
library(doParallel)




#make table output for initial 215 results:
library(stringr)
rivannaOutput <- read.csv("Ly6a_TF_list_WITH_HUMAN_PREDICTED_GENESETS_completed.csv")
#clean up table:
rivannaOutput <- subset(rivannaOutput, select= -c(1)) #x2
rivannaOutput["Num_of_HumanGenes"] <- NA
rivannaOutput["COPY_HumanGeneSets"] <- NA
rivannaOutput$COPY_HumanGeneSets <- rivannaOutput$Human_Genes_with_TF_sequence #save copy of human gene set to reference after merge

for (j in 1:nrow(rivannaOutput)) {
  rivannaOutput$Num_of_HumanGenes[[j]] <- str_count(rivannaOutput$Human_Genes_with_TF_sequence[[j]], "\\w+") #"\\w+" allows you to count the WORDS with str_count
  print(paste0(j, "_", rivannaOutput$TF_name[[j]], "_done"))
}

write.csv(rivannaOutput, "Ly6a_TF_list_WITH_HUMAN_PREDICTED_GENESETS_completed.csv")




topGeneSets <- as.data.frame(table(rivannaOutput$Human_Genes_with_TF_sequence))
colnames(topGeneSets)[1] <- "GeneSet"
colnames(topGeneSets)[2] <- "Frequency"

export <- merge(x=topGeneSets, y=rivannaOutput, by.x='GeneSet', by.y='Human_Genes_with_TF_sequence')
write.csv(export, "Ly6a_TF_list_HUMAN_PREDICTED_GENESETS_num_112321.csv")

export <- read.csv("Ly6a_TF_list_HUMAN_PREDICTED_GENESETS_num_112321.csv")
export <- subset(export, select=-c(1))

final <- export[!is.na(export$COPY_HumanGeneSets), ]

#determining which bulk/lesion cues upregulate gene sets (need to repeat for downregulated)
library(tidyverse)
library(dplyr)


geneList <- as.data.frame(final$Mouse_Genes_with_TF_sequence[[282]]) #need to keep as mouse genes for transcriptome match
colnames(geneList)[1] <- "gene"
geneList <- geneList %>%
     tidyr::separate_rows("gene", sep = ", ") #this separates comma-separated data in a row into own columns :)
geneList <- geneList[!duplicated(geneList), ] #removes duplicated mouse genes
geneList1 <- as.data.frame(geneList$gene)
colnames(geneList1)[1] <- "gene"


setwd("~/Desktop/R Pdgfrb:Dual/2021-08-07_Merge_GFA_E080_E081_E082_WT.YFP.POS.MEDIA.cells_Gene 2/bulkRNAseq_P-0.05_Fold-1")
L_bulk_up = list.files(pattern="UP_Gene") #make sure previous versions of code/files not saved (otherwise out of bounds error)
df_bulk_up <- lapply(setNames(L_bulk_up, tools::file_path_sans_ext(basename(L_bulk_up))), read.csv) 
my_gene <- c("gene")
genes_bulk_up = lapply(df_bulk_up, "[", , my_gene)

pairs_up <- expand.grid(val1 = genes_bulk_up, val2 = geneList1)

results_up <- mapply(
  function(x, y) length(intersect(x, y)), pairs_up$val1,pairs_up$val2
)


NAMES_results_up <- mapply(
  function(x, y) intersect(x, y), pairs_up$val1,pairs_up$val2
)




