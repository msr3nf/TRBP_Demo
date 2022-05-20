setwd("~/Desktop")
#the first part of this script is similar to the 04292022 R script, but the second half shows how we computed TF changes between WT and mutated Myh11 PPRs
##because this is the original script written, there is some repetitive code that was generated while troubleshooting but may be helpful for users attempting to understand how our algorithms were developed
library(jsonlite)
library(tidyverse)
library(httr)
library(stringr)
library(data.table)
library(purrr)
library(GenomicFeatures)
library(dplyr)

Myh11_TFs_Mut_ALL <- read.csv("Myh11_TFs_Mut_ALL.csv")
WholeGenome_List <- read.csv("WholeGenome_List_mm10_20kPPRseq_03162022.csv")

length(unique(WholeGenome_List$))
WholeGenome_List <- WholeGenome_List[!duplicated(WholeGenome_List$PPR_sequence), ]
WholeGenome_List <- WholeGenome_List[,-1]
WholeGenome_List <- WholeGenome_List[,-c(2,3,7,10)]
WholeGenome_List <- WholeGenome_List[,-c(8)]

dir.create("~/Desktop/R Pdgfrb:Dual/20kb_PPR_mm10")
setwd("~/Desktop/R Pdgfrb:Dual/20kb_PPR_mm10")

write.csv(WholeGenome_List[1:5000,], "WholeGenome_List_mm10_20kPPRseq_pt1.csv")

unique_Mut_ALL_TFseq <- as.data.frame(unique(Myh11_TFs_Mut_ALL$TF_sequence))
colnames(unique_Mut_ALL_TFseq)[1] <- "TF_sequence"
unique_Mut_ALL_TFseq$Mouse_Genes_with_TF_sequence <- NA

seqMatch <- function(x){
  result <- WholeGenome_List[WholeGenome_List$PPR_sequence %like% x, ]
  result_mouse_gene <- paste(as.character(unique(result[[18]])), sep="' '", collapse=", ") #col #16 has all gene names
  return(result_mouse_gene)
}

numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1
options(mc.cores=numCores)

unique_Mut_ALL_TFseq$Mouse_Genes_with_TF_sequence[1:10] <- unlist(mclapply(unique_Mut_ALL_TFseq$TF_sequence[1:10], seqMatch))
write.csv(unique_Mut_ALL_TFseq, "unique_Mut_ALL_TFseq_TEST_05032022.csv")

#control shift C comments out (and in) chunks of code :)

##alternate method:

library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

unique_Mut_ALL_TFseq <- read.csv("Myh11_TFs_Mut_ALL.csv")
Myh11_TF_list <- read.csv("Myh11_TF_list_wMouseGenesONLY_031820222_FINAL.csv")
Mut_All_genes <- subset(Myh11_TF_list, subset = TF_sequence %in% unique_Mut_ALL_TFseq$TF_sequence)


{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("mmusculus_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoMouse <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "entrezgene_id", 
                                    "mgi_symbol", "description", "external_gene_name"), mart=ensembl)
}


mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
homologues <- getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"), filters = "external_gene_name", values = annoMouse$mgi_symbol, mart = mouse)

percentTable <- homologues
length(unique(percentTable$mgi_symbol))
colnames(percentTable)[1] <- "mgi_symbol" #note that these are NOT exlusively protein-coding genes
percentTable["Frequency"] <- NA
percentTable["Frequency"] <- 0
percentTable["Associated_Myh11_TFs"] <- NA

resultsTable <- Mut_All_genes

dataSet <- resultsTable
geneList <- (unlist(strsplit(dataSet$Mouse_Genes_with_TF_sequence, " "))) #add TF_name
geneList <- gsub(",", "", gsub("([a-zA-Z]),", "\\1 ", geneList))
geneList <- as.data.table(geneList)
colnames(geneList)[1] <- "gene"
geneList$gene <- str_trim(geneList$gene, "right")

frequencyTable <- table(geneList$gene)
frequencyTable <- as.data.table(frequencyTable) #this is my frequency table!
colnames(frequencyTable)[1] <- "mgi_symbol"
colnames(frequencyTable)[2] <- "Frequency"

library(dplyr)
finalTable <- percentTable %>% inner_join(frequencyTable,by="mgi_symbol")
finalTable <- finalTable[,-3]
colnames(finalTable)[4] <- "Frequency"
finalTable['% of Myh11 TFs that bind Mouse gene’s PPR'] <- NA


for (i in 1:nrow(finalTable)) {
  finalTable$`% of Myh11 TFs that bind Mouse gene’s PPR`[[i]] <- finalTable$Frequency[[i]]/(length(Mut_All_genes$X))*100
}

write.csv(finalTable, "Mut_All_geneFreq_05042022.csv")
finalTable <- read.csv("Mut_All_geneFreq_05042022.csv")

#####rerun for just Myh11 PPR of interest --> contains all identified reg. regions#####

#step1: identify where regulatory region is located
Myh11_TFs_WT <- read.csv("Myh11_TF_list_wMouseGenesONLY_031820222_FINAL.csv")

#upstream SRF region= chr16:14289628-14289644
#downstream SRF region= chr16:chr16:14292512-14292528

# Myh11_TFs_upperBound <- subset(Myh11_TFs_WT, subset = chromStart > 14289628 & chromEnd < 14292528)
# rm(Myh11_TFs_upperBound) # rm() clears specific items from environment 
# rm(Myh11_TFs_narrow)
Myh11_TFs_WTnarrow <- subset(Myh11_TFs_WT, subset = chromStart >= 14289628 & chromEnd <= 14292528) #9700 TFs in 2900 bp region

{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("mmusculus_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoMouse <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "entrezgene_id", 
                                    "mgi_symbol", "description", "external_gene_name"), mart=ensembl)
}


mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
homologues <- getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"), filters = "external_gene_name", values = annoMouse$mgi_symbol, mart = mouse)

percentTable <- homologues
length(unique(percentTable$mgi_symbol))
colnames(percentTable)[1] <- "mgi_symbol" #note that these are NOT exlusively protein-coding genes
percentTable["Frequency"] <- NA
percentTable["Frequency"] <- 0
percentTable["Associated_Myh11_TFs"] <- NA

resultsTable <- Myh11_TFs_WTnarrow
resultsTable <- 

dataSet <- resultsTable
geneList <- (unlist(strsplit(dataSet$Mouse_Genes_with_TF_sequence, " "))) #add TF_name
geneList <- gsub(",", "", gsub("([a-zA-Z]),", "\\1 ", geneList))
geneList <- as.data.table(geneList)
colnames(geneList)[1] <- "gene"
geneList$gene <- str_trim(geneList$gene, "right")

frequencyTable <- table(geneList$gene)
frequencyTable <- as.data.table(frequencyTable) #this is my frequency table!
colnames(frequencyTable)[1] <- "mgi_symbol"
colnames(frequencyTable)[2] <- "Frequency"

library(dplyr)
finalTable <- percentTable %>% inner_join(frequencyTable,by="mgi_symbol")
finalTable <- finalTable[,-3]
colnames(finalTable)[4] <- "Frequency"
finalTable['% of Myh11 TFs that bind Mouse gene’s PPR'] <- NA


for (i in 1:nrow(finalTable)) {
  finalTable$`% of Myh11 TFs that bind Mouse gene’s PPR`[[i]] <- finalTable$Frequency[[i]]/(length(Myh11_TFs_WTnarrow$X))*100
}

write.csv(finalTable, "narrowPPR_Myh11_WT_05062022.csv")
narrowPPR_Myh11_WT <- read.csv("narrowPPR_Myh11_WT_05062022.csv")

####determine for Mut_ALL PPR
Myh11_TFs_Mut_ALL <- read.csv("Myh11_TFs_Mut_ALL.csv")
Myh11_TFs_Mut_ALLnarrow <- subset(Myh11_TFs_Mut_ALL, subset = chromStart >= 14289628 & chromEnd <= 14292528) #9700 TFs in 2900 bp region

{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("mmusculus_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoMouse <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "entrezgene_id", 
                                    "mgi_symbol", "description", "external_gene_name"), mart=ensembl)
}


mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
homologues <- getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"), filters = "external_gene_name", values = annoMouse$mgi_symbol, mart = mouse)

percentTable <- homologues
length(unique(percentTable$mgi_symbol))
colnames(percentTable)[1] <- "mgi_symbol" #note that these are NOT exlusively protein-coding genes
percentTable["Frequency"] <- NA
percentTable["Frequency"] <- 0
percentTable["Associated_Myh11_TFs"] <- NA

resultsTable <- subset(Myh11_TFs_WTnarrow, subset = TF_sequence %in% Myh11_TFs_Mut_ALLnarrow$TF_sequence)

length(unique(resultsTable$TF_name)) #687 unique TFs
length(unique(Myh11_TFs_WTnarrow$TF_name)) #690 unique TFs


dataSet <- resultsTable
geneList <- (unlist(strsplit(dataSet$Mouse_Genes_with_TF_sequence, " "))) #add TF_name
geneList <- gsub(",", "", gsub("([a-zA-Z]),", "\\1 ", geneList))
geneList <- as.data.table(geneList)
colnames(geneList)[1] <- "gene"
geneList$gene <- str_trim(geneList$gene, "right")

frequencyTable <- table(geneList$gene)
frequencyTable <- as.data.table(frequencyTable) #this is my frequency table!
colnames(frequencyTable)[1] <- "mgi_symbol"
colnames(frequencyTable)[2] <- "Frequency"

library(dplyr)
finalTable <- percentTable %>% inner_join(frequencyTable,by="mgi_symbol")
finalTable <- finalTable[,-3]
colnames(finalTable)[4] <- "Frequency"
finalTable['% of Myh11 TFs that bind Mouse gene’s PPR'] <- NA


for (i in 1:nrow(finalTable)) {
  finalTable$`% of Myh11 TFs that bind Mouse gene’s PPR`[[i]] <- finalTable$Frequency[[i]]/(length(Myh11_TFs_WTnarrow$X))*100
}



write.csv(finalTable, "narrowPPR_Myh11_mutAll_05072022.csv")
finalTable_MutAll <- read.csv("narrowPPR_Myh11_mutAll_05072022.csv")
finalTable_MutAll <- finalTable_MutAll[order(-finalTable_MutAll$X..of.Myh11.TFs.that.bind.Mouse.gene.s.PPR),] 
write.csv(finalTable_MutAll, "narrowPPR_Myh11_mutAll_05072022.csv")

finalTable_WT <- read.csv("narrowPPR_Myh11_WT_05062022.csv")
finalTable_WT <- finalTable_MutWT[order(-finalTable_MutWT$X..of.Myh11.TFs.that.bind.Mouse.gene.s.PPR),] 
write.csv(finalTable_WT, "narrowPPR_Myh11_WT_05062022.csv")

#####verify accessibility with ATAC-seq#####
PDGFBB_genes <- read.csv("output_PDGFBB_ATAC_wSeq.csv")
HC_SMC <- read.csv("carotid_VSMC_cluster_markers copy.csv") 
finalTable_MutAll['HumanCarotid_scRNA-seq_cluster'] <- NA
finalTable_MutAll <- finalTable_MutAll[,-c(1:2,5)]
finalTable_WT['HumanCarotid_scRNA-seq_cluster'] <- NA
finalTable_WT <- finalTable_WT[,-c(1:2,5)]


for (j in 1:nrow(finalTable_MutAll)) {
  for (i in 1:ncol(HC_SMC)) {
    HC_SMC_common <- length(intersect(finalTable_MutAll$hsapiens_homolog_associated_gene_name[[j]], HC_SMC[[i]]))
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      finalTable_MutAll$`HumanCarotid_scRNA-seq_cluster`[[j]] <- paste0(finalTable_MutAll$`HumanCarotid_scRNA-seq_cluster`[[j]], clusterNum, sep=" ", collapse=", ")
    }
  }
  print(paste0(finalTable_MutAll$mgi_symbol[[j]],"_",j,"_done"))
}


for (j in 1:nrow(finalTable_WT)) {
  for (i in 1:ncol(HC_SMC)) {
    HC_SMC_common <- length(intersect(finalTable_WT$hsapiens_homolog_associated_gene_name[[j]], HC_SMC[[i]]))
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      finalTable_WT$`HumanCarotid_scRNA-seq_cluster`[[j]] <- paste0(finalTable_WT$`HumanCarotid_scRNA-seq_cluster`[[j]], clusterNum, sep=" ", collapse=", ")
    }
  }
  print(paste0(finalTable_WT$mgi_symbol[[j]],"_",j,"_done"))
}

humanPred_WT <- subset(finalTable_WT, subset = `HumanCarotid_scRNA-seq_cluster` != " " )
humanPred_mutAll <- subset(finalTable_MutAll, subset = `HumanCarotid_scRNA-seq_cluster` != " " )

length(intersect(PDGFBB_genes$gene, humanPred_WT$hsapiens_homolog_associated_gene_name))
ATACval_WT <- as.data.frame(intersect(PDGFBB_genes$gene, humanPred_WT$hsapiens_homolog_associated_gene_name))
colnames(ATACval_WT)[1] <- "gene"

humanPred_scATAC_WT <- subset(humanPred_WT, subset = `hsapiens_homolog_associated_gene_name` %in% ATACval_WT$gene )
length(unique(humanPred_scATAC_WT$hsapiens_homolog_associated_gene_name)) #138 unique genes
write.csv(humanPred_scATAC_WT, "Myh11_humanPred_scATAC_WT.csv")
write.csv(humanPred_WT, "Myh11_scRNAC_WT.csv")
write.csv(humanPred_mutAll, "Myh11_scRNAC_mutAll.csv")



###repeat for MutAll
length(intersect(PDGFBB_genes$gene, humanPred_mutAll$hsapiens_homolog_associated_gene_name))
ATACval_mutAll <- as.data.frame(intersect(PDGFBB_genes$gene, humanPred_mutAll$hsapiens_homolog_associated_gene_name))
colnames(ATACval_mutAll)[1] <- "gene"

humanPred_scATAC_mutAll <- subset(humanPred_mutAll, subset = `hsapiens_homolog_associated_gene_name` %in% ATACval_mutAll$gene )
length(unique(humanPred_scATAC_mutAll$hsapiens_homolog_associated_gene_name)) #138 unique genes
write.csv(humanPred_scATAC_mutAll, "Myh11_humanPred_scATAC_mutAll.csv")

setdiff(humanPred_scATAC_mutAll[2:21,]$mgi_symbol, humanPred_scATAC_WT[2:21,]$mgi_symbol) #"loss" of Adgrf5

##NEED TO REPEAT FOR TFS
TOMTOM <- read.csv("TOMTOM_TFnames_PDGFBB_01242022.csv")
Ly6a_TFs <- read.csv("Ly6a_TF_list_20kbPPRseq_03282022.csv")
humanTFs <- as.data.frame(intersect(TOMTOM$TF_name, Ly6a_TFs$TF_name))
colnames(humanTFs)[1] <- "TF_name"

####need WT SMCs +MURINE scRNA-seq data

percentDiff <- finalTable_WT
colnames(percentDiff)[7] <- "WT_percent_Myh11_TFs" 
percentDiff <- percentDiff[,-c(1:2,4:5)]

colnames(finalTable_MutAll)[7] <- "MutAll_percent_Myh11_TFs" 
finalTable_MutAll <- finalTable_MutAll[,-c(1:2,4:5)]

percentDiff <- inner_join(percentDiff, finalTable_MutAll, by="mgi_symbol")
percentDiff <- percentDiff[,-c(2,6)]

percentDiff <- percentDiff[!duplicated(percentDiff$mgi_symbol),]
percentDiff['Δfrequency'] <- NA
percentDiff['Δpercent'] <- NA

i=1
for (i in 1:nrow(percentDiff)){
  percentDiff$Δfrequency[i] <- percentDiff$Frequency.x[i]-percentDiff$Frequency.y[i]
  percentDiff$Δpercent[i] <- percentDiff$WT_percent_Myh11_TFs[i]-percentDiff$MutAll_percent_Myh11_TFs[i]
  
}
write.csv(percentDiff, "Myh11_WT_mutAll_percentFreqDiff_noRNAseq_05082022.csv")
###compare to mouse RNA-seq
library(Seurat)
library(tidyverse)
library(MAST)
library(data.table)
library(biomaRt)
library(ReactomePA)
library(org.Mm.eg.db)
library(formattable)
library(dplyr)
library(DESeq2)
library(stringr)
library(recount)
library(gage)
library(pathview)
library(gageData)
library(pheatmap)
library(clusterProfiler)
all_plaque <- readRDS(file = "healthy.disease.OwensLab.only.making.final.graphs15dim.new.idents.new.idents.v2.rds")

datasetUSE = subset(all_plaque, subset = origin == "SMC_Klf4_WT_eYFP_Positive")
{
  print("Loading annotation tables")
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("mmusculus_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  anno <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "entrezgene_id","mgi_symbol", "description", "external_gene_name"), mart=ensembl)
}

firstgroup=1
for (firstgroup in 1:7) {  
      #DE analysis on different groups
      groupMarkers = FindMarkers(datasetUSE, ident.1 = firstgroup, min.pct = 0.25)
      groupMarkers = setDT(groupMarkers, keep.rownames = TRUE)[]
      colnames(groupMarkers)[1] <- "gene"
      print("groupMarkers")
      
      #Annotation of markers
      anno_groupMarkers <-inner_join(groupMarkers, anno, by=c("gene"= "external_gene_name")) 
      print("anno_groupMarkers")
      
      #Subsetting for significant genes
      anno_up_groups = subset(anno_groupMarkers, subset = p_val_adj <0.05 & avg_log2FC>0)
      anno_down_groups = subset(anno_groupMarkers, subset = p_val_adj <0.05 & avg_log2FC<0)
      print("subset")
  
      Final_UP_Group = anno_up_groups
      Final_DOWN_Group = anno_down_groups
        
      write.csv(Final_UP_Group, paste0("scRNA_up_groups.",firstgroup,".csv"))
      write.csv(Final_DOWN_Group, paste0("scRNA_down_groups.",firstgroup,".csv"))
      
      print(paste("Completed",firstgroup))
}


#####determine which genes experience most TF reg changes#####
setwd("~/Desktop/R Pdgfrb:Dual/scRNA-seq_murineMarkers")
my_col = c("gene")
L_sc_up = list.files(pattern=paste0("scRNA_up")) #make sure previous versions of code/files not saved (otherwise out of bounds error)
df_sc_up <- lapply(setNames(L_sc_up, tools::file_path_sans_ext(basename(L_sc_up))), read.csv) 
desc_sc_up = lapply(df_sc_up, "[", , my_col)

L_sc_down = list.files(pattern=paste0("scRNA_down")) #make sure previous versions of code/files not saved (otherwise out of bounds error)
df_sc_down <- lapply(setNames(L_sc_down, tools::file_path_sans_ext(basename(L_sc_down))), read.csv) 
desc_sc_down = lapply(df_sc_down, "[", , my_col)
setwd('..')

percentDiff['Murine_RNA-seq_cluster_UP'] <- NA
percentDiff['Murine_RNA-seq_cluster_DOWN'] <- NA


for (j in 1:nrow(percentDiff)) {
  for (i in 1:length(desc_sc_up)){
    if (percentDiff$mgi_symbol[[j]] %in% desc_sc_up[[i]] == TRUE)
      
      percentDiff$`Murine_RNA-seq_cluster_UP`[[j]] <- paste0(percentDiff$`Murine_RNA-seq_cluster_UP`[[j]],"cluster_",i,",")
  }
  print(paste0(percentDiff$mgi_symbol[[j]],"_",i,"_done"))
}



for (j in 1:nrow(percentDiff)) {
  for (i in 1:length(desc_sc_down)){
    if (percentDiff$mgi_symbol[[j]] %in% desc_sc_down[[i]] == TRUE)
      
      percentDiff$`Murine_RNA-seq_cluster_DOWN`[[j]] <- paste0(percentDiff$`Murine_RNA-seq_cluster_DOWN`[[j]],"cluster_",i,",")
  }
  print(paste0(percentDiff$mgi_symbol[[j]],"_",i,"_done"))
}

write.csv(percentDiff, "percentDiff_table_Myh11_WT_mutAll_murine_scRNA-seq.csv")

percentDiff <- subset(percentDiff, subset = `Murine_RNA-seq_cluster_UP` != " " | `Murine_RNA-seq_cluster_DOWN` != " " )
#723 genes potentially modulated with murine SMC significance

write.csv(percentDiff, "FINALpercentDiff_table_Myh11_WT_mutAll_murine_scRNA-seq.csv")
