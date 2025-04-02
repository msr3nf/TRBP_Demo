#####meme suite in R#####
BiocManager::install("memes")
BiocManager::install("JASPAR2022")
{
 library(memes)
library(magrittr)
library(Biostrings)
library(universalmotif)
library(JASPAR2022)
library(TFBSTools)

curated.list <- read.csv("/Users/Downloads/scRNA_up_groups.2.csv")
curated.list <- subset(curated.list, subset=avg_log2FC >= log2(1.5)) #yields 30 genes 
 #identify PPR regions based on output of Step 1
WholeGenome_List <- read.csv("/Users/Desktop/WholeGenome_List_mm10_20kPPRseq_03222025.trimmed.csv")
curated.genomes <- subset(WholeGenome_List, external_gene_name %in% curated.list$gene)
curated.genomes$PPR_sequence <- as.character(curated.genomes$PPR_sequence)
curated.genomes <- curated.genomes[!duplicated(curated.genomes$PPR_sequence),] #remove any duplicate sequences
curated.genomes$fasta <- paste(">", curated.genomes$X, "\n", curated.genomes$PPR_sequence, sep="")
colnames(curated.genomes)

#converting PPR sequences to BStringSet for runMeme
fasta_sequences <- BStringSet(as.character(curated.genomes$PPR_sequence)) #has 89 sequences
names(fasta_sequences) <- paste0("seq_", seq_along(fasta_sequences)) #need unique identifier
#now identify the motifs
meme_results <- runMeme(fasta_sequences, alph = "dna", parse_genomic_coord = FALSE, nmotifs = 10, minw=6, minw=50) #based on the recommended parameters from the website

#now use tomtom to get predicted TFs (vignette: "Motif Comparison using TomTom")
motifs <- lapply(meme_results$name, function(x) create_motif(x))
motifs
sapply(motifs, class)
#specify the database to use for tomtom
opts <- list(species = 10090) #this is for mouse; you can change this for species of interest
jaspar_motifs <- getMatrixSet(JASPAR2022, opts)

jaspar_motifs_file <- tempfile(fileext = ".meme")
write_meme(jaspar_motifs, jaspar_motifs_file)
options(meme_db = jaspar_motifs_file)


tomtom_results <- runTomTom(motifs)

#Extract "best match name" (the TF) for all motifs in tomtom_results --> this is the final TF list
best_match_TF_names <- sapply(tomtom_results, function(x) x$best_match_name)
best_match_consensus <- sapply(tomtom_results, function(x) x$consensus) #note the IUPAC nucleotide code
best_match_TF_df <- data.frame(best_match_name = unique(unlist(best_match_TF_names)))
best_match_consensus_df <- data.frame(best_consensus = unique(unlist(best_match_consensus)))
final.result <- cbind(best_match_TF_df, best_match_consensus_df)
                              
}
