setwd("~/Desktop/") #change to directory of interest

#####make gene frequency table#####
library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

resultsTable <- read.csv("Myh11_TF_list_wMouseGenesONLY_031820222_FINAL.csv")
resultsTable <- resultsTable[,-1]

{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("mmusculus_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoMouse <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "entrezgene_id", 
                                    "mgi_symbol", "description", "external_gene_name"), mart=ensembl)
}

#already set up with known human homologues
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
homologues <- getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"), filters = "external_gene_name", values = annoMouse$mgi_symbol, mart = mouse)

percentTable <- homologues
length(unique(percentTable$mgi_symbol))
colnames(percentTable)[1] <- "mgi_symbol" #note that these are NOT exlusively protein-coding genes
percentTable["Frequency"] <- NA
percentTable["Frequency"] <- 0
percentTable["Associated_Myh11_TFs"] <- NA

dataSet <- resultsTable
geneList <- (unlist(strsplit(dataSet$Mouse_Genes_with_TF_sequence, " "))) #add TF_name
geneList <- gsub(",", "", gsub("([a-zA-Z]),", "\\1 ", geneList))
geneList <- as.data.table(geneList)
colnames(geneList)[1] <- "gene"
geneList$gene <- str_trim(geneList$gene, "right")
write.csv(geneList, "geneList_03182022.csv")

#####making a frequency table#####
{for (i in 1:nrow(percentTable)) {
  gene_val <- percentTable$mgi_symbol[[i]]
  percentTable$Frequency <- length(which(geneList$gene==gene_val))
  print(paste0(i, "_", gene_val, "_done"))
  
}}

#faster method to make table:
frequencyTable <- table(geneList$gene)
frequencyTable <- as.data.table(frequencyTable) #this is my frequency table!
colnames(frequencyTable)[1] <- "mgi_symbol"
colnames(frequencyTable)[2] <- "Frequency"

finalTable <- percentTable %>% inner_join(frequencyTable,by="mgi_symbol")
finalTable <- finalTable[,-3]
colnames(finalTable)[4] <- "Frequency"
finalTable['% of Myh11 TFs that bind Mouse gene’s PPR'] <- NA


for (i in 1:nrow(finalTable)) {
  finalTable$`% of Myh11 TFs that bind Mouse gene’s PPR`[[i]] <- finalTable$Frequency[[i]]/75532*100
}

finalTable <- percentTable %>% inner_join(frequencyTable,by="mgi_symbol")
finalTable <- finalTable[,-3]
colnames(finalTable)[4] <- "Frequency"
finalTable['% of Myh11 TFs that bind Mouse gene’s PPR'] <- NA

for (i in 1:nrow(finalTable)) {
  finalTable$`% of Myh11 TFs that bind Mouse gene’s PPR`[[i]] <- finalTable$Frequency[[i]]/75532*100
}

#OPTIONAL: determine whether human orthologs of mouse genes appear in human scRNA-seq data --> users can load their own -omics data of interest
finalTable['HumanCarotid_scRNA-seq_cluster'] <- NA
HC_SMC <- read.csv("carotid_VSMC_cluster_markers copy.csv") #has top 50 marker genes for each cluster 

for (j in 1:nrow(finalTable)) {
  for (i in 1:ncol(HC_SMC)) {
    HC_SMC_common <- length(intersect(finalTable$hsapiens_homolog_associated_gene_name[[j]], HC_SMC[[i]]))
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      finalTable$`HumanCarotid_scRNA-seq_cluster`[[j]] <- paste0(finalTable$`HumanCarotid_scRNA-seq_cluster`[[j]], clusterNum, sep=" ", collapse=", ")
    }
  }
  print(paste0(finalTable$mgi_symbol[[j]],"_",j,"_done"))
}


write.csv(finalTable, "Myh11_orthGenesbyTFpercent_HCscRNAseq_03182022_20kPPR.csv") #file is available for users to download from Myh11_Input_Files on the TRBP GitHub repo --> can continue with pipeline
finalTable <- read.csv("Myh11_orthGenesbyTFpercent_HCscRNAseq_03182022_20kPPR.csv")
finalTable <- finalTable[,-1]
#rename columns if desired
colnames(finalTable)[4] <- 'Percent_of_Myh11_TFs_that_bind_Mouse_gene_PPR' 

#adding in TF information (i.e. which Myh11-assocaited TFs bind to each gene)
finalTable['Number_Myh11_TFs'] <- NA
colnames(finalTable)[6] <- "Number_UNIQUE_Myh11_TFs"
finalTable <- finalTable[,-3]
finalTable['Associated_Myh11_TFs'] <- NA

for (i in 5746:nrow(finalTable)) {
  gene_val <- finalTable$mgi_symbol[[i]]
  finalTable$Number_Myh11_TFs[[i]] <- length(which(dataSet$Mouse_Genes_with_TF_sequence==gene_val))
  print(paste0(i, "_", gene_val, "_done"))
  
}

write.csv(finalTable, "Myh11_orthGenesbyTFpercent_andNumUNIQUE_TFs_HCscRNAseq_03182022_20kPPR.csv") #also available in the repo!
Myh11_Freq_Table <- read.csv("Myh11_orthGenesbyTFpercent_andNumUNIQUE_TFs_HCscRNAseq_03182022_20kPPR.csv")

#####making heat maps for TFs of interest #####
Myh11_TF_list <- read.csv("Myh11_TF_list_wMouseGenesONLY_031820222_FINAL.csv") #from line23
seqData <- subset(Myh11_TF_list, subset = TF_name %in% c("SRF", "KLF4", "Yy1", "ELK1", "MEF2", "MAX", "MSY-1", "TEF-1", "RTEF-1")) #Select TFs of interest
seqData <- seqData[,8:10]

library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

seqData <- separate_rows(seqData, Mouse_Genes_with_TF_sequence, convert = TRUE) #separate_rows is a tidyr command :)
finalData <- subset(seqData, subset = Mouse_Genes_with_TF_sequence %in% c("Myh11", "Macf1", "Frmd4a", "Tmcc1", "Ank2", "Foxp1", "Etl4", "Nedd4l", "Dlg2", "Pde4b", "Acta2", "Tagln")) #Select genes of interest (i.e. top 10 genes with greatest % of Myh11-specific TF binding sites)
table(finalData$TF_sequence)
length(unique(finalData$TF_sequence))
finalData['Num_TF_binding_sites'] <- NA

finalMatrix <- finalData[,-2]
finalMatrix <- distinct(finalMatrix)
length(unique(finalMatrix$Mouse_Genes_with_TF_sequence))


for (i in 1:nrow(finalMatrix)) {
  gene_val <- finalMatrix$Mouse_Genes_with_TF_sequence[[i]]
  TF_val <- finalMatrix$TF_name[[i]]
  data <- subset(finalData, subset = Mouse_Genes_with_TF_sequence == gene_val)
  data2 <- subset(data, subset = TF_name == TF_val)
  finalMatrix$Num_TF_binding_sites[[i]] <- length(data2$TF_sequence)
}

#the finalMatrix contains a raw file with number of binding sequences per gene and TF
colnames(finalMatrix)[2] <- "Gene" #modify column names to ease data presentation in heat map

ggp <- ggplot(finalMatrix, aes(TF_name, Gene)) +                          
  geom_tile(aes(fill = Num_TF_binding_sites))
ggp

ggp + theme(axis.text.x = element_text(size = 10))     
ggp + theme(axis.text.y = element_text(size = 6))  

contractile_genes <- subset(Myh11_Freq_Table, subset = mgi_symbol %in% c("Myh11", "Acta2", "Tagln", "Col18a1"))
contractile_genes <- contractile_genes[,-1]
write.csv(contractile_genes, "contractile_genes_Myh11_FreqTable.csv") #to save data for future


