setwd("~/Desktop/R Pdgfrb:Dual")
Ly6a_TF_list <- read.csv("Ly6a_TF_list_HuMoCHECKED_021112022.csv")
Ly6a_TF_list <- Ly6a_TF_list[,-1]

######repeat PPR identification for ALL genes in mm10#####
#get list of all genes in human genome from anno
library(Biostrings)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
library(biomaRt)
library(org.Mm.eg.db)
library(dplyr)

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
mouseGenome <- subset(mouseGenome, subset = entrezgene_id != "NA") #29135 "valid" entrez gene id's

library(jsonlite)
library(tidyverse)
library(httr)
library(stringr)
library(data.table)
library(purrr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)
library(dplyr)
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

#run on Rivanna
for (i in 1:nrow(WholeGenome_List)) {
  if (WholeGenome_List$strand[[i]] == "+") {
    chrom <- WholeGenome_List$seqnames[[i]]
    chromStart <- WholeGenome_List$start[[i]]-1000
    chromEnd <- WholeGenome_List$start[[i]]+500
    Gene_name <- WholeGenome_List$external_gene_name[[i]]
    my.dnastring <- as.character(Biostrings::getSeq(Mmusculus, chrom, chromStart, chromEnd))
    
  }
  
  if (WholeGenome_List$strand[[i]] == "-") {
    chrom <- WholeGenome_List$seqnames[[i]]
    chromStart <- WholeGenome_List$end[[i]]-500
    chromEnd <- WholeGenome_List$end[[i]]+1000
    Gene_name <- WholeGenome_List$external_gene_name[[i]]
    my.dnastring <- as.character(Biostrings::getSeq(Mmusculus, chrom, chromStart, chromEnd))
    
  }
  print(paste0(i, "_", Gene_name, "_done"))
  
  WholeGenome_List$PPR_sequence[[i]] <- my.dnastring
  
}
length(unique(WholeGenome_List$external_gene_name)) #19321 unique gene IDs

write.csv(WholeGenome_List, "WholeGenome_List_mm10_PPRseq_02112022.csv")
#need to clean up list

####identify mouse genes with Ly6a-specific TFs
Ly6a_TF_list <- read.csv("Ly6a_TF_list_HuMoCHECKED_021112022.csv")
Ly6a_TF_list <- Ly6a_TF_list[,-1]

Ly6a_TF_list[,12] <-NA 
colnames(Ly6a_TF_list)[12] <- "Mouse_Genes_with_TF_sequence"
Ly6a_TF_list[,13] <-NA 
colnames(Ly6a_TF_list)[13] <- "Human_Genes_with_TF_sequence"

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

j=3

for (j in 94:nrow(Ly6a_TF_list)) {
  result <- WholeGenome_List[WholeGenome_List$PPR_sequence %like% Ly6a_TF_list$TF_sequence[[j]], ]
  result_mouse_gene <- paste(as.character(unique(result[[17]])), sep="' '", collapse=", ") #col #16 has all gene names
  if (length(result$external_gene_name) > 0) {
    Ly6a_TF_list$Mouse_Genes_with_TF_sequence[[j]] <- result_mouse_gene
    #eliminate repeats in this results
    homologues = getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"), filters = "external_gene_name", values = result[[17]], mart = mouse)
    Ly6a_TF_list$Human_Genes_with_TF_sequence[[j]] <-paste(as.character(homologues[[2]]), sep="' '", collapse=", ") #this also eliminates repeats! for human set
    #convert "result" to human ortholog  
  }
  print(paste0(j, "_", Ly6a_TF_list$MouseTF_name[[j]], "_done"))
}



Ly6a_TF_list <- read.csv("Ly6a_TF_list_HuMo_predictions_02122022.csv") #2295
Ly6a_TF_list <- Ly6a_TF_list[,-1]
Ly6a_TF_list[1,13]

##########narrowing down TFs/genes##########
Ly6a_TF_list <- subset(Ly6a_TF_list, Human_Genes_with_TF_sequence != "NA") #2215
Ly6a_TF_list <- subset(Ly6a_TF_list, Human_Genes_with_TF_sequence != "NA, NA") #2095
length(unique(Ly6a_TF_list$MouseTF_name)) #408 unique TFs
#Ly6a_TF_list[,14] <- NA
#colnames(Ly6a_TF_list)[14] <- "Num_of_HumanGenes"
#Num_of_HumanGenes

#make gene frequency table
library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("mmusculus_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoMouse <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "entrezgene_id", 
                                    "mgi_symbol", "description", "external_gene_name"), mart=ensembl)
}

resultsTable <- as.data.frame(annoMouse$mgi_symbol)
resultsTable <- resultsTable[!apply(resultsTable == "", 1, all),]
resultsTable <- as.data.frame(resultsTable)
colnames(resultsTable)[1] <- "mgi_symbol" #note that these are NOT exlusively protein-coding genes
resultsTable["Frequency"] <- NA
resultsTable["Frequency"] <- 0
resultsTable["Associated_Ly6a_TFs"] <- NA

dataSet <- Ly6a_TF_list

i=1
j=1
k=19503
  
for (i in 1:nrow(dataSet)) {
  
  geneList <- unique(unlist(strsplit(dataSet$Mouse_Genes_with_TF_sequence[[i]], " ")))
  geneList <- gsub(",", "", gsub("([a-zA-Z]),", "\\1 ", geneList))
  geneList <- as.data.table(geneList)
  colnames(geneList)[1] <- "gene"
  geneList$gene <- str_trim(geneList$gene, "right") #eliminates weird spacing that interferes with if statement below   
  
  for (j in 1:nrow(geneList)) {
    for (k in 1:nrow(resultsTable)) {
      if (geneList$gene[[j]] == resultsTable$mgi_symbol[[k]]) {
        
        resultsTable$Frequency[[k]] <-  resultsTable$Frequency[[k]] + 1
        TFname <- dataSet$MouseTF_name[[i]]
        
        if (i==1) {
          resultsTable$Associated_Ly6a_TFs[[k]] <- paste0(TFname)
        }
        else {
          resultsTable$Associated_Ly6a_TFs[[k]] <- paste0(resultsTable$Associated_Ly6a_TFs[[k]], ",", TFname)
        }
      }
    }
  }
  print(paste0(i, "_", dataSet$MouseTF_name[[i]], "_done"))
}

#this loop takes the longest (13 hours) to run!

data <- read.csv("Ly6a_topGenesTFs_percents_02222022.csv")
data <- data[,-1]
data <- subset(data, subset = Frequency !=0) #55,000 --> 19474 mouse genes that bind at least one Ly6a-specific TF in PPR
data[,4] <- NA
data <- data[,-5]
colnames(data)[4] <- "% of Ly6a TFs that bind Mouse gene’s PPR"
for (i in 1:nrow(data)) {
  data$`% of Ly6a TFs that bind Mouse gene’s PPR`[[i]] <- data$Frequency[[i]]/2092*100
}

write.csv(data, "Ly6a_topGenesTFs_percents_02222022.csv")
##########working with ATAC-seq data##########
TOMTOM <- read.csv("TOMTOM_TFnames_PDGFBB_01242022.csv") #direct output of MEME Suite TOMTOM tool 
#PDGF_genes <- read.csv("PDGF-BB_HCASMC_ATACseq_genes_DONE.csv") #don't use this file
PDGFBB_genes <- read.csv("output_PDGFBB_ATAC_wSeq.csv")
PDGFBB_genes <- subset(PDGFBB_genes, subset= score >= 100) #because MEME was done on enrichment score > 100
length(unique(PDGFBB_genes$gene))
#determine which Ly6a-specific TFs are associated with PDGF-ATAC-seq

commonTFs <- as.data.frame(intersect(TOMTOM$TF_name, Ly6a_TF_list$MouseTF_name))
write.csv(commonTFs, "common_Ly6aTFs_PDGFBB_HCASMC_TOMTOM_TFs_03022022.csv")
colnames(commonTFs)[1] <- "commonTF_name"
library(data.table)
Ly6a_TF_list1 <- Ly6a_TF_list
Ly6a_TF_list1 <- setDT(Ly6a_TF_list1)[MouseTF_name %chin% TOMTOM$TF_name] #this is how you select rows to keep based on whether or not they appear in another table/list

length(intersect(Ly6a_TF_list1$MouseTF_name, TOMTOM$TF_name)) #check that there are still 81 TFs
Ly6a_TF_list1[["Num_in_PDGFBB_ATAC-seq"]] <- NA #need two brackets?
Ly6a_TF_list1[["Percent_humanOrtholog_PDGFBB_ATAC-seq"]] <- NA #need two brackets?
Ly6a_TF_list1[["Common_humanOrtholog_PDGFBB_ATAC-seq"]] <- NA

i=1
for (i in 1:nrow(Ly6a_TF_list1)) {
  geneList <- as.data.frame(unique(unlist(strsplit(Ly6a_TF_list1$Human_Genes_with_TF_sequence[[i]],","))))
  #geneList1 <- as.data.frame(unlist(strsplit(Ly6a_TFs$COPY_HumanGeneSets[[i]],",")))
  colnames(geneList)[1] <- "humanGene"
  geneList <- subset(geneList, subset = humanGene != "")
  geneList <- subset(geneList, subset = humanGene != " " )
  geneList$humanGene <- gsub('\\s+', '', geneList$humanGene) #sometimes, there are spaces before the name of each gene
  
  Ly6a_TF_list1$`Num_in_PDGFBB_ATAC-seq`[[i]] <- length(intersect(geneList$humanGene, PDGFBB_genes$gene))
  Ly6a_TF_list1$`Percent_humanOrtholog_PDGFBB_ATAC-seq`[[i]] <- ((Ly6a_TF_list1$`Num_in_PDGFBB_ATAC-seq`[[i]])/(length(geneList$humanGene)))*100
  
  if (Ly6a_TF_list1$`Percent_humanOrtholog_PDGFBB_ATAC-seq`[[i]] != 0) {
    Ly6a_TF_list1$`Common_humanOrtholog_PDGFBB_ATAC-seq`[[i]] <- paste(as.character(intersect(geneList$humanGene, PDGFBB_genes$gene)), sep="' '", collapse=", ")
  }
  print(paste0(i, "_", Ly6a_TF_list1$MouseTF_name[[i]], "_done"))
  
}

write.csv(Ly6a_TF_list1, "Ly6a_TF_list_PDGFBB_hits_allCommonTFs_03022022.csv")

Ly6a_TF_list1 <- subset(Ly6a_TF_list1, subset = `Num_in_PDGFBB_ATAC-seq` > 0)
write.csv(Ly6a_TF_list1, "Ly6a_TF_list_PDGFBB_hits_TFs_morethan1hit_03022022.csv")

Ly6a_TF_list1 <- read.csv("Ly6a_TF_list_PDGFBB_hits_TFs_morethan1hit_03022022.csv")
Ly6a_TF_list1 <- Ly6a_TF_list1[,-1]

Ly6a_TF_list2 <- read.csv("Ly6a_TF_list_PDGFBB_hits_allCommonTFs_03022022.csv")
Ly6a_TF_list2 <- Ly6a_TF_list2[,-1]

#Make simplified tables 
Table <- Ly6a_TF_list1[,c(8,16)]
colnames(Table)[1] <- "Ly6a TF name"
colnames(Table)[2] <- "Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs"

simplified_Table <- as.data.frame(unique(Table$`Ly6a TF name`))
simplified_Table[2] <- NA
colnames(simplified_Table)[1] <- "Ly6a TF name"
colnames(simplified_Table)[2] <- "Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs"

i=18
for (i in 1:nrow(simplified_Table)) {
  
  Ly6a_TF <- simplified_Table$`Ly6a TF name`[i]
  genelist <- subset(Table, subset = `Ly6a TF name`  == Ly6a_TF)
  simplified_Table$`Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs`[[i]] <- paste0(unique(genelist$`Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs`), sep=" ", collapse=",")
  
  #convert "result" to human ortholog  
  print(paste0(i, "_", simplified_Table$`Ly6a TF name`[[i]], "_done"))
  
}

write.csv(simplified_Table, "FINAL_candidateTF_gene_table_03032022.csv")
simplified_Table <- read.csv("FINAL_candidateTF_gene_table_03032022.csv")
simplified_Table <- simplified_Table[,-1]


##########working with human carotid artery RNA-seq data##########
#first compare commonTFs list to see if ANYTHING common using Ly6a_TF_list2


HC_SMC <- read.csv("carotid_VSMC_cluster_markers copy.csv") #has top 50 marker genes for each cluster 
commonTFs <- as.data.frame(unique(Ly6a_TF_list1$HumanTF_gene_name))
colnames(commonTFs) <- "commonTF_name"

#need to determine which TFs and genes are relevant
commonTFs[2] <- NA
colnames(commonTFs)[2] <- "HumanCarotid_scRNA-seq_clusterNumber"
colnames(commonTFs)[1] <- "TF_name"

Gene_List <- as.data.frame(unique(unlist(strsplit(Table$`Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs`, " "))))
colnames(Gene_List)[1] <- "Gene_name"
Gene_List$Gene_name<-as.character(gsub("\\,", "", Gene_List$Gene_name))
Gene_List <- unique(Gene_List$Gene_name)
Gene_List <- as.data.frame(Gene_List)
colnames(Gene_List)[1] <- "Gene_name"
Gene_List[2] <- NA
colnames(Gene_List)[2] <- "HumanCarotid_scRNA-seq_clusterNumber" 

for (j in 1:nrow(commonTFs)) {
  for (i in 1:ncol(HC_SMC)) {
    
    HC_SMC_common <- length(intersect(commonTFs$TF_name[[j]], HC_SMC[[i]]))
    
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      commonTFs$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(commonTFs$`HumanCarotid_scRNA-seq_clusterNumber`[[j]], clusterNum, sep=" ", collapse=", ")
    }
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}


for (j in 1:nrow(Gene_List)) {
  for (i in 1:ncol(HC_SMC)) {
    
    HC_SMC_common <- length(intersect(Gene_List$Gene_name[[j]], HC_SMC[[i]]))
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      
      Gene_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(Gene_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]], clusterNum, sep=" ", collapse=", ")
    }
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}

write.csv(commonTFs, "commonTFs_in_HCarotid_scRNA-seq_03032022.csv")
Gene_List1 <- Gene_List[!(is.na(Gene_List$`HumanCarotid_scRNA-seq_clusterNumber`) | Gene_List$`HumanCarotid_scRNA-seq_clusterNumber`==""), ]
Gene_List1$`HumanCarotid_scRNA-seq_clusterNumber`<-as.character(gsub("\\NA", "", Gene_List1$`HumanCarotid_scRNA-seq_clusterNumber`))

write.csv(Gene_List1, "commonGENES_in_HCarotid_scRNA-seq_03032022.csv")





