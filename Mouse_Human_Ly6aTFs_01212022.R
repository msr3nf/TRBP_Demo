setwd("~/Desktop/R Pdgfrb:Dual")
Ly6a_mouseTFs <- read.csv("Ly6a_TF_list_human_mouse_wholeGenome_geneList.csv")
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

write.csv(df, "Ly6aTFs_mouse_human_gene_sets_87TFs.csv") #THIS IS THE FILTERED TABLE TO USE



