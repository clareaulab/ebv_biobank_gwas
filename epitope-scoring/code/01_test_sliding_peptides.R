library(Biostrings)
library(data.table)
library(dplyr)

generate_sliding_windows <- function(seq, window_size = 9) {
  n <- nchar(seq)
  if (n < window_size) return(character(0))
  sapply(1:(n - window_size + 1), function(i) substr(seq, i, i + window_size - 1))
}

# Read in fasta file
fasta_file <- "../data/EBV_proteins_genename.fasta" 
proteins <- readAAStringSet(fasta_file)

# Create a dataframe of 9-mers
window_size <- 9
result <- data.frame(Protein = character(), Peptide = character(), stringsAsFactors = FALSE)

for (i in seq_along(proteins)) {
  prot_id <- names(proteins)[i]
  sequence <- as.character(proteins[[i]])
  peptides <- generate_sliding_windows(sequence, window_size)
  if (length(peptides) > 0) {
    temp_df <- data.frame(Protein = rep(prot_id, length(peptides)), Peptide = peptides, stringsAsFactors = FALSE)
    result <- rbind(result, temp_df)
  }
}

result <- result %>% unique()


