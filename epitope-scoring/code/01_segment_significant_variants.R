library(Biostrings)
library(dplyr)
library(data.table)
library(BuenColors)
library(seqinr)
options(warn=-1)

make_wt_and_mutant <- function(sequence, sequence_name = "", pos = 1, wt = "X", mut = "Y", kmer = 8){
  
  middle_seq <- toupper(as.character(substr(sequence, pos - kmer+1, pos + kmer-1)))
  lapply(1:(kmer), function(i){
    mutant_name = paste0(">", sequence_name, "_", as.character(kmer),"mer_MUT_", as.character(kmer - i + 1), "_",  mut)
    mutant = paste0(substr(middle_seq, i, kmer-1), mut, substr(middle_seq, kmer+1, kmer+i-1))
    wildtype_name = paste0(">", sequence_name,"_", as.character(kmer),"mer_WT_", as.character(kmer - i + 1), "_", wt)
    wildtype = paste0(substr(middle_seq, i, kmer-1), wt, substr(middle_seq, kmer+1, kmer+i-1))
    data.frame(out = c(mutant_name, mutant, wildtype_name, wildtype))
  }) %>% rbindlist
  
}

# read in the last four VUS from viral-sequences
vus_keep <- fread("../../viral-sequences/output/VUS_table.tsv") %>% tail(4)

si <- seqinr::read.fasta("../data/EBV_proteins_genename.fasta", as.string = TRUE)

# Verify WT alleles match
substr(si[["BALF2"]], 613, 613)
substr(si[["RPMS1"]], 51, 51)
substr(si[["BNRF1"]], 694, 694)
substr(si[["BALF2"]], 317, 317)

write_xmer <- function(kmer){
  vusdf <- rbind(
    make_wt_and_mutant(si[["BALF2"]], "BALF2", 613, "I", "V", kmer = kmer),
    make_wt_and_mutant(si[["RPMS1"]], "RPMS1", 51, "D", "N", kmer = kmer),
    make_wt_and_mutant(si[["BNRF1"]], "BNRF1", 694, "P", "H", kmer = kmer),
    make_wt_and_mutant(si[["BALF2"]], "BALF2", 317, "V", "M", kmer = kmer)
  )
  write.table(vusdf, file = paste0("../output/vus_peptides/wt_mut_peptides_kmer", kmer, ".fasta"),
              sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

lapply(c(8:11, 15), write_xmer)

# these fasta files can be used as input to NetMHCpan/NetMHCIIpan


