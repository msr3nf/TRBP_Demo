#####meme suite in R#####
BiocManager::install("memes")
library(memes)
library(magrittr)
library(Biostrings)
library(universalmotif)
BiocManager::install("JASPAR2022")
library(JASPAR2022)
library(TFBSTools)
packageVersion("TFBSTools")

curated.list <- read.csv("/Users/mahimareddy/Downloads/scRNA_up_groups.2.csv")
curated.list <- subset(curated.list, subset=avg_log2FC >= log2(1.5)) #yields 30 genes 
 #identify PPR regions based on output of Step 1
WholeGenome_List <- read.csv("/Users/mahimareddy/Desktop/WholeGenome_List_mm10_20kPPRseq_03222025.trimmed.csv")
curated.genomes <- subset(WholeGenome_List, external_gene_name %in% curated.list$gene)
curated.genomes$PPR_sequence <- as.character(curated.genomes$PPR_sequence)
curated.genomes <- curated.genomes[!duplicated(curated.genomes$PPR_sequence),] #remove any duplicate sequences
curated.genomes$fasta <- paste(">", curated.genomes$X, "\n", curated.genomes$PPR_sequence, sep="")
colnames(curated.genomes)

#converting PPR sequences to BStringSet for runMeme
fasta_sequences <- BStringSet(curated.genomes$PPR_sequence) #has 89 sequences
names(fasta_sequences) <- paste0("seq_", seq_along(fasta_sequences)) #need unique identifier
#now identify the motifs
meme_results <- runMeme(fasta_sequences, alph = "dna", parse_genomic_coord = FALSE, nmotifs = 10, minw=6, minw=50) #based on the recommended parameters from the website

#now use tomtom to get predicted TFs (vignette: "Motif Comparison using TomTom")
motifs <- lapply(meme_results$name, function(x) create_motif(x))
motifs
#specify the database to use for tomtom
opts <- list(species = 10090) #this is for mouse; you can change this for species of interest
jaspar_motifs <- getMatrixSet(JASPAR2022, opts)
options(meme_db = jaspar_motifs)
jaspar_file <- tempfile(fileext = ".meme")
export_meme(jaspar_motifs, jaspar_file)
class(motifs)
sapply(motifs, class)

tomtom_results <- runTomTom(motifs, database = jaspar_motifs)
