setwd("~/Desktop)
#####ensuring PPR sequences are determined correctly#####

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

list <- as.data.frame("Ly6a")
colnames(list)[1] <- "gene"

anno_Ly6a <-inner_join(list, anno, by=c("gene"= "external_gene_name")) #Entrez gene id is 110454

library(jsonlite)
library(tidyverse)
library(httr)
library(stringr)
library(data.table)
library(purrr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
list1 <- as.character("110454") 

txbygene <- transcriptsBy(txdb, "gene")[list1]
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
UTR_List1 <- merge(UTR_List, map1, by = "value") #5'UTR is on neg strand, so TSS is +1000, -500 BP from chrEnd
UTR_List1 <- subset(UTR_List1, exon_rank == 1)
#there are multiple 5' UTRs, so find average
mean(UTR_List1[,4])
mean(UTR_List1[,5])
#UCSC Genome uses chr15:74994878-74998031, so 74998031 is the TSS
PPR_end <- 74998031+1000
PPR_start <- 74998031-500

responseTF = fromJSON("https://api.genome.ucsc.edu/getData/track?genome=mm10;track=hub_186875_JASPAR2022_TFBS_mm10;chrom=chr15;start=74997531;end=74999031")
Ly6a_TF_list <- responseTF$hub_186875_JASPAR2022_TFBS_mm10
colnames(Ly6a_TF_list)[4] <- "MatrixID"
colnames(Ly6a_TF_list)[7] <- "TF_name"
Ly6a_TF_list <- subset(Ly6a_TF_list, subset = strand == "-") #why does a TF associated with neg strand gene bind to positive strand?
write.csv(Ly6a_TF_list, "Ly6a_TF_list_NEW_02102022.csv")

Ly6a_TF_list['TF_sequence'] <- NA

for (i in 1:nrow(Ly6a_TF_list)) {
  chrom <- Ly6a_TF_list$chrom[[i]]
  chromStart <- Ly6a_TF_list$chromStart[[i]]
  chromEnd <- Ly6a_TF_list$chromEnd[[i]]
  TF_name <- Ly6a_TF_list$MouseTF_name[[i]]
  
  my.dnastring <- as.character(Biostrings::getSeq(Mmusculus, chrom, chromStart, chromEnd))
  Ly6a_TF_list$TF_sequence[[i]] <- my.dnastring
  
  print(paste0(i, "_", TF_name, "_done"))
}
write.csv(Ly6a_TF_list, "Ly6a_TF_list_PPRseq_02102022.csv")


#make sure these TFs are MOUSE-ONLY
library(dplyr)
library(data.table)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

colnames(Ly6a_TF_list)[7] <- "MouseTF_name"
Ly6a_TF_list[8] <- NA
Ly6a_TF_list[9] <- NA
Ly6a_TF_list[10] <- NA
colnames(Ly6a_TF_list)[8] <- "MouseTF_gene_name"
colnames(Ly6a_TF_list)[9] <- "HumanTF_gene_name"
colnames(Ly6a_TF_list)[10] <- "HumanTF_name"


for (j in 1:nrow(Ly6a_TF_list)) {
  humanGenes =  getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"), filters = "external_gene_name", values = Ly6a_TF_list$MouseTF_name[[j]], mart = mouse)
  Ly6a_TF_list$MouseTF_gene_name[[j]] <-paste(as.character(humanGenes[[1]]), sep="' '", collapse=", ")
  Ly6a_TF_list$HumanTF_gene_name[[j]] <-paste(as.character(humanGenes[[2]]), sep="' '", collapse=", ")
  
  #convert "result" to human ortholog  
  print(paste0(j, "_", Ly6a_TF_list$MouseTF_name[[j]], "_done"))
  
}

write.csv(Ly6a_TF_list, "Ly6a_TF_list_HuMo_sorted_02102022.csv")
Ly6a_TF_list <- read.csv("Ly6a_TF_list_HuMo_sorted_02102022.csv")
Ly6a_TF_list <- Ly6a_TF_list[,-1]
Ly6a_TF_list <- Ly6a_TF_list[,-11]
Ly6a_TF_list <- subset(Ly6a_TF_list, subset = MouseTF_gene_name != "")
write.csv(Ly6a_TF_list, "Ly6a_TF_list_HuMoCHECKED_021112022.csv")









#####//////////////////////////////////////#####

Ly6a_mouseTFs1 <- read.csv("Ly6a_TF_list_human_mouse_wholeGenome_geneList.csv")
length(unique(Ly6a_mouseTFs$TF_name))
Ly6a_mouseTFs <- Ly6a_mouseTFs[,-c(1)]
colnames(Ly6a_mouseTFs)[5] <- "MouseTF_name"
mouse_humanTFs <- Ly6a_mouseTFs[,c(1:5)]
mouse_humanTFs[6] <- NA
mouse_humanTFs[7] <- NA
mouse_humanTFs[8] <- NA
colnames(mouse_humanTFs)[6] <- "MouseTF_gene_name"
colnames(mouse_humanTFs)[7] <- "HumanTF_gene_name"
colnames(mouse_humanTFs)[8] <- "HumanTF_name"

#JASPAR documentation: TF names are based on "standardized Entrez gene symbol"

library(dplyr)
library(data.table)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


for (j in 1:nrow(mouse_humanTFs)) {
  humanGenes =  getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"), filters = "external_gene_name", values = mouse_humanTFs$MouseTF_name[[j]], mart = mouse)
  mouse_humanTFs$MouseTF_gene_name[[j]] <-paste(as.character(humanGenes[[1]]), sep="' '", collapse=", ")
  mouse_humanTFs$HumanTF_gene_name[[j]] <-paste(as.character(humanGenes[[2]]), sep="' '", collapse=", ")
  
  
  #convert "result" to human ortholog  
  print(paste0(j, "_", mouse_humanTFs$MouseTF_name[[j]], "_done"))
  
}

step1 <- read.csv("mouse_humanTFs_step1_01212022.csv")
step1 <- step1[,-c(1)]

humanTF <- read.delim("Homo_sapiens_TF.txt")
mouseTF <- read.delim("Mus_musculus_TF.txt")

mouseTFs_noHumanOrtho <-as.data.frame(setdiff(step1$MouseTF_name, step2$MouseTF_name))


step1[8,7]

step2 <- subset(step1, subset = HumanTF_gene_name != "")
length(unique(step2$MouseTF_name))

write.csv(step2, "Ly6a_TFnames_wHumanOrth.csv")

i=2508
#repeat loop until row lengths for step2 and ly6a_TFs line up? need to optimize this... 
for (i in 1:nrow(Ly6a_mouseTFs)) {
  value <- length(intersect((Ly6a_mouseTFs$TF_name[[i]]), step2$MouseTF_name))
  
  if (value == 0) {
    Ly6a_mouseTFs <- Ly6a_mouseTFs[-c(i),]
    
  }
  
  
  #convert "result" to human ortholog  
  print(paste0(i, "_", Ly6a_mouseTFs$TF_name[[i]], "_done"))
  
}

length(setdiff(Ly6a_mouseTFs$TF_name, step2$MouseTF_name))
length(intersect(Ly6a_mouseTFs$TF_name, step2$MouseTF_name))

write.csv(Ly6a_mouseTFs, "Ly6a_TFnames_wHumanOrth.csv")
Ly6a_mouseTFs <- read.csv("Ly6a_TFnames_wHumanOrth.csv")
Ly6a_mouseTFs <- Ly6a_mouseTFs[,-c(1)]

#filtering values for PowerPoint diagrams
Ly6a_mouseTFs[2,10]

i=2
#repeat loop until same number of rows appears --> need to revise this loop

for (i in 1:nrow(Ly6a_mouseTFs)) {
  if (Ly6a_mouseTFs[i,10] == "") {
    Ly6a_mouseTFs <- Ly6a_mouseTFs[-c(i),]
  }
  print(paste0(i, "_", Ly6a_mouseTFs$TF_name[[i]], "_done"))
}

MEME_TFs <- c("ETS1", "FOS", "GATA6", "GSC2", "KLF4", "KLF5", "NFIX", "SOX10", "SP2", "TFAP2A", "ZEB1", "ZNF263")
MEME_TFs <- as.data.frame(MEME_TFs)
length(intersect(MEME_TFs$MEME_TFs, Ly6a_mouseTFs$TF_name))
setdiff(MEME_TFs$MEME_TFs, Ly6a_mouseTFs$TF_name)

setwd("~/Desktop/R Pdgfrb:Dual/12212021 Lab Meeting Data")

df <- read.csv("Ly6a_TF_list_HUMAN_PREDICTED_GENESETS_num_112321.csv")




df$COPY_HumanGeneSets

setwd("~/Desktop/R Pdgfrb:Dual/12212021 Lab Meeting Data")
Ly6a_TFs <- read.csv("CLEANED_Ly6a_TF_list_HUMAN_PREDICTED_GENESETS_num_112321 copy.csv")

for (i in 1:nrow(df)) {
  value <- length(intersect((df$TF_name[[i]]), step2$MouseTF_name))
  
  if (value == 0) {
    df <- df[-c(i),]
    
  }
  print(paste0(i, "_", Ly6a_mouseTFs$TF_name[[i]], "_done"))
}

test <- subset(df, subset = TF_name == "ZNF263") #checking it works
length(unique(df$TF_name)) #have more than 1 human ortholog

setwd("~/Desktop/R Pdgfrb:Dual")

TOMTOM <- read.table("tomtom_01242022.tsv") #direct output of MEME Suite TOMTOM tool 
names(TOMTOM) <- TOMTOM[1,]
TOMTOM <- TOMTOM[-c(1),]
TOMTOM[11] <- NA
colnames(TOMTOM)[11] <- "TF_name"



#convert matrix ids to TF family/names
library(jsonlite)
library(data.table)
library(purrr)
library(stringr)
library(curl)

#use a JASPAR API (homo sapiens=species)
responseTF = fromJSON("https://jaspar.genereg.net/api/v1/matrix/MA1596.1/")
responseTF$name

i=1
for (i in 1:nrow(TOMTOM)) {

  matrixID <- TOMTOM$Target_ID[i]
  responseTF = fromJSON(paste("https://jaspar.genereg.net/api/v1/matrix/",matrixID,"/",sep = ""))
  TOMTOM$TF_name[i] <- responseTF$name
  
  #convert "result" to human ortholog  
  print(paste0(i, "_", TOMTOM$Target_ID[[i]], "_done"))
  
}

write.csv(TOMTOM, "TOMTOM_TFnames_PDGFBB_01242022.csv")

length(unique(TOMTOM$TF_name))
df <- subset(df, subset = Num_of_HumanGenes > 0)
length(unique(df$TF_name))
length(intersect(TOMTOM$TF_name, df$TF_name))
candidate_TFs <- as.data.frame(intersect(TOMTOM$TF_name, df$TF_name))
colnames(candidate_TFs)[1] <- "humanTF_name"

write.csv(candidate_TFs, "candidate_TFs_01242022.csv")
length(unique(candidate_TFs$humanTF_name)) #check these are unique

#########what PDGF-BB genes bind these TFs?##########
df_2 <- subset(df, subset = TF_name == candidate_TFs$humanTF_name )

for (i in 1:nrow(TOMTOM)) {
  value <- length(intersect((TOMTOM$TF_name[[i]]), candidate_TFs$humanTF_name))
  
  if (value == 0) {
    TOMTOM <- TOMTOM[-c(i),]
    
  }
  print(paste0(i, "_", TOMTOM$TF_name[[i]], "_done"))
}
  
length(unique(TOMTOM$TF_name))

for (i in 1:nrow(df)) {
  value <- length(intersect((df$TF_name[[i]]), candidate_TFs$humanTF_name))
  
  if (value == 0) {
    df <- df[-c(i),]
    
  }
  print(paste0(i, "_", df$TF_name[[i]], "_done"))
}

length(unique(df$TF_name))

write.csv(df, "Ly6aTFs_mouse_human_gene_sets_87TFs.csv") #THIS IS THE FILTERED TABLE TO USE --> look at ETS1 sample

################work done on 01/25/2022################

setwd("~/Desktop/R Pdgfrb:Dual")
Ly6a_TFs <- read.csv("Ly6aTFs_mouse_human_gene_sets_87TFs.csv")
length(unique(Ly6a_TFs$TF_name))
Ly6a_TFs["Num_in_PDGFBB_ATAC-seq"] <- NA
Ly6a_TFs["Percent_humanOrtholog_PDGFBB_ATAC-seq"] <- NA
Ly6a_TFs["Common_humanOrtholog_PDGFBB_ATAC-seq"] <- NA
Ly6a_TFs <- Ly6a_TFs[,-1]
PDGFBB_genes <- read.csv("output_PDGFBB_ATAC_wSeq.csv")
PDGFBB_genes <- subset(PDGFBB_genes, subset= score >= 100) #because MEME was done on enrichment score > 100
length(unique(PDGFBB_genes$gene)) #this table contains repeats 

i=1
for (i in 1:nrow(Ly6a_TFs)) {
  
  geneList <- as.data.frame(unique(unlist(strsplit(Ly6a_TFs$COPY_HumanGeneSets[[i]],","))))
  #geneList1 <- as.data.frame(unlist(strsplit(Ly6a_TFs$COPY_HumanGeneSets[[i]],",")))
  colnames(geneList)[1] <- "humanGene"
  geneList <- subset(geneList, subset = humanGene != "")
  geneList <- subset(geneList, subset = humanGene != " " )

  Ly6a_TFs$`Num_in_PDGFBB_ATAC-seq`[[i]] <- length(intersect(geneList$humanGene, PDGFBB_genes$gene))
  Ly6a_TFs$`Percent_humanOrtholog_PDGFBB_ATAC-seq`[[i]] <- ((Ly6a_TFs$`Num_in_PDGFBB_ATAC-seq`[[i]])/(length(geneList$humanGene)))*100
  
  if (Ly6a_TFs$`Percent_humanOrtholog_PDGFBB_ATAC-seq`[[i]] != 0) {
    Ly6a_TFs$`Common_humanOrtholog_PDGFBB_ATAC-seq`[[i]] <- intersect(geneList$humanGene, PDGFBB_genes$gene)
    
  }
  print(paste0(i, "_", Ly6a_TFs$TF_name[[i]], "_done"))
  
}

write.csv(Ly6a_TFs, "Ly6a_TFs_PDGFBB_HumanGeneSets_01252022.csv")

Ly6a_TFs <- read.csv("Ly6a_TFs_PDGFBB_HumanGeneSets_01252022.csv")
copy_Ly6a_TFs <- Ly6a_TFs

filtered_Ly6a_TFs <- subset(copy_Ly6a_TFs, subset = copy_Ly6a_TFs$Num_in_PDGFBB_ATAC.seq != 0)
filtered_Ly6a_TFs <- filtered_Ly6a_TFs[,-1]
length(unique(filtered_Ly6a_TFs$TF_name))

write.csv(filtered_Ly6a_TFs, "FINAL_Ly6a_TFs_PDGFBB_HumanGeneSets_01252022.csv")

################making summary tables################
Table <- filtered_Ly6a_TFs[,c(9,18)]
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

write.csv(simplified_Table, "FINAL_candidateTF_gene_table_01252022.csv")
simplified_Table <- read.csv("FINAL_candidateTF_gene_table_01252022.csv")
simplified_Table <- simplified_Table[,-1]

################comparing to human scRNA-seq data################
HC_SMC <- read.csv("carotid_VSMC_cluster_markers copy.csv") #has top 50 marker genes for each cluster 
TF_List <- as.data.frame(simplified_Table$Ly6a.TF.name)
TF_List[2] <- NA
colnames(TF_List)[2] <- "HumanCarotid_scRNA-seq_clusterNumber"
colnames(TF_List)[1] <- "TF_name"

Gene_List <- as.data.frame(unique(Table$`Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs`))
Gene_List[2] <- NA
colnames(Gene_List)[2] <- "HumanCarotid_scRNA-seq_clusterNumber"
colnames(Gene_List)[1] <- "Gene_name"



#length(unique(Table$`Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs`)) #checking # of UNIQUE Gene (29)

for (j in 1:nrow(TF_List)) {
  for (i in 1:ncol(HC_SMC)) {
    
    HC_SMC_common <- length(intersect(TF_List$TF_name[[j]], HC_SMC[[i]]))
    
    
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      
      TF_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(TF_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]], clusterNum, sep=" ", collapse=", ")
      
    }
    
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}


for (j in 1:nrow(Gene_List)) {
  for (i in 1:ncol(HC_SMC)) {
    
    HC_SMC_common <- length(intersect(Gene_List$Gene_name[[j]], HC_SMC[[i]]))
    
    
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      
      Gene_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(TF_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]], clusterNum, sep=" ", collapse=", ")
      
    }
    
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}

#####testing other types of for loops#####
List = list()
for (i in 1:ncol(HC_SMC)) {
  for (j in 1:nrow(TF_List)) {
    
    HC_SMC_common <- length(intersect(TF_List$TF_name[[j]], HC_SMC[[i]]))
    
    if (HC_SMC_common !=0) {
      
      clusterNum <- colnames(HC_SMC)[[i]]
      List[[length(List)+1]] = clusterNum
      df_List <- as.data.frame(List)
      TF_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(df_List[1,], collapse=", ")
      
    }
    
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}

List = list()
for (i in 1:ncol(HC_SMC)) {
  for (j in 1:nrow(Gene_List)) {
    
    HC_SMC_common <- length(intersect(Gene_List$Gene_name[[j]], HC_SMC[[i]]))
    
    if (HC_SMC_common !=0) {
      
      clusterNum <- colnames(HC_SMC)[[i]]
      List[[length(List)+1]] = clusterNum
      df_List <- as.data.frame(List)
      Gene_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(df_List[1,], collapse=", ")
      
    }
    
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}

######identifying bulk conditions that up/down regulate Ly6a######

geneList <- "Ly6a"

geneList2 <- as.data.frame(geneList)
colnames(geneList2)[1] <- "gene"

setwd("~/Desktop/R Pdgfrb:Dual/2021-08-07_Merge_GFA_E080_E081_E082_WT.YFP.POS.MEDIA.cells_Gene 2/bulkRNAseq_P-0.05_Fold-1")
L_bulk_up = list.files(pattern="UP_Gene") #make sure previous versions of code/files not saved (otherwise out of bounds error)
df_bulk_up <- lapply(setNames(L_bulk_up, tools::file_path_sans_ext(basename(L_bulk_up))), read.csv) 
my_gene <- c("gene")
genes_bulk_up = lapply(df_bulk_up, "[", , my_gene)

pairs_up <- expand.grid(val1 = genes_bulk_up, val2 = geneList2)

results_up <- mapply(
  function(x, y) length(intersect(x, y)), pairs_up$val1,pairs_up$val2
)


NAMES_results_up <- mapply(
  function(x, y) intersect(x, y), pairs_up$val1,pairs_up$val2
)


######following up with K and G genes of interest######


LS_SMC <- read.csv("carotid_VSMC_late stage copy.csv") #has top 50 marker genes for each cluster 
LS_SMC <- LS_SMC[,-c(4)] #cleaning up file read in
names(LS_SMC) <- as.matrix(LS_SMC[1,])
names(murine_VSMC_MAC_markers) <- as.matrix(murine_VSMC_MAC_markers[1, ])
LS_SMC <- LS_SMC[-1,]

TF_List <- as.data.frame(simplified_Table$Ly6a.TF.name)
TF_List[2] <- NA
colnames(TF_List)[2] <- "HumanCarotid_scRNA-seq_clusterNumber"
colnames(TF_List)[1] <- "TF_name"

Gene_List <- as.data.frame(unique(Table$`Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs`))
Gene_List[2] <- NA
colnames(Gene_List)[2] <- "HumanCarotid_scRNA-seq_clusterNumber"
colnames(Gene_List)[1] <- "Gene_name"
#Ly6a gene is in mSMC6 cluster (late stage murine lesions)

#convert mouse and TF lists to mouse IDs?
library(data.table)
library(dplyr)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

Gene_List["Mouse_Ortholog"] <- NA
TF_List["Mouse_Ortholog"] <- NA


human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouseGenes =  getBM(attributes = c("external_gene_name", "mmusculus_homolog_associated_gene_name"), filters = "external_gene_name", values = geneList$humanGene, mart = human)
geneList <- mouseGenes
colnames(geneList)[1] <- "Human_Gene"
colnames(geneList)[2] <- "Mouse_Ortholog"


#length(unique(Table$`Ly6a Human Orthologous Genes in PDGF-BB-treated HCASMCs`)) #checking # of UNIQUE Gene (29)

for (j in 1:nrow(TF_List)) {
  for (i in 1:ncol(LS_SMC)) {
    
    LS_SMC_common <- length(intersect(TF_List$TF_name[[j]], LS_SMC[[i]]))
    
    
    if (LS_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      
      TF_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(TF_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]], clusterNum, sep=" ", collapse=", ")
      
    }
    
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}


for (j in 1:nrow(Gene_List)) {
  for (i in 1:ncol(HC_SMC)) {
    
    HC_SMC_common <- length(intersect(Gene_List$Gene_name[[j]], HC_SMC[[i]]))
    
    
    if (HC_SMC_common !=0) {
      clusterNum <- colnames(HC_SMC)[[i]]
      
      Gene_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]] <- paste0(TF_List$`HumanCarotid_scRNA-seq_clusterNumber`[[j]], clusterNum, sep=" ", collapse=", ")
      
    }
    
  }
  print(paste0(colnames(HC_SMC)[[i]], "_done"))
  
}

#automate nanodrop steps



