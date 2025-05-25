library(Biostrings)
library(dplyr)
library(data.table)
library(BuenColors)
options(warn=-1)


# import and munge variant data
aou_df <- fread("../data/AOU_chrEBV.tsv.gz")
ukb_df <- fread("../data/UKB_chrEBV.tsv.gz")

# import VEP annotations 
ebv_dt_vep <- readRDS("../output/full_EBV_annotation.rds") %>%
  mutate(id = paste0("chrEBV", Location, REF, ">", ALT))

transform_df <- function(df){
  df$REF <- factor(df$REF, levels = c("A", "C", "G", "T"))
  smol <- df[,c("BP", "REF", "TOTAL")]
  smol_melt = smol %>% mutate(
    A = df$A / df$TOTAL,
    C = df$C / df$TOTAL,
    G = df$G / df$TOTAL,
    T = df$T / df$TOTAL
  ) %>% reshape2::melt(id.vars = c("BP", "REF", "TOTAL")) %>% filter(BP != variable) %>%
    mutate(id = paste0("chrEBV", BP, REF, ">", variable))
  
}
aou_dft <- merge(transform_df(aou_df), ebv_dt_vep, by = "id")
ukb_dft <- merge(transform_df(ukb_df), ebv_dt_vep, by = "id")

# Import reference genome and look for peptides
dna <- unname(readDNAStringSet("../data/ebv_reference.fa"))

## removing the nucleotides prior to the reading frame start position.
dna6 <-c(
  subseq(c(dna, reverseComplement(dna)), start=1),
  subseq(c(dna, reverseComplement(dna)), start=2),
  subseq(c(dna, reverseComplement(dna)), start=3)
)[c(1,3,5,2,4,6)]

protein6 <- Biostrings::translate(dna6)

parse_coordinates_peptide <- function(ebv_peptide){
  # Query all reading frames for peptide
  hit <-data.frame(vmatchPattern(ebv_peptide, protein6))
  #return(length(hit$group))
  if((length(hit$group) == 0)){ # no match
    return("no hits")
  }
  if((length(hit$group) == 2)){ # multiple matches
    return("multiple hits")
  }
  if((hit$group %in% c(1,2,3))){
    # No reverse complment
    ff <- case_when(
      hit$group == 1 ~ 1,
      hit$group == 2 ~ 2,
      hit$group == 3 ~ 3
    )
    idx1 <-  hit$start*3 + ff -3 # kinda manual to find these
    idx2 <- idx1 + 3*nchar(ebv_peptide)
    
    return(data.frame(
      query = ebv_peptide, 
      result = Biostrings::translate(DNAStringSet(substr(dna, idx1, idx2)), no.init.codon = TRUE), 
      idx1, idx2, strand = "+", DNAseq = substr(dna, idx1, idx2)
    ))
    
  }
  
  if((hit$group %in% c(4,5,6))){
    
    # yes reverse complement
    ff <- case_when(
      hit$group == 6 ~ 2,
      hit$group == 5 ~ 1,
      hit$group == 4 ~ 0
    )
    idx2 <- nchar(dna) - hit$start*3 - ff + 3 # kinda manual to find these
    idx1 <- idx2 - 3*nchar(ebv_peptide)
    return(data.frame(
      query = ebv_peptide, 
      result = Biostrings::translate(reverseComplement(DNAStringSet(substr(dna, idx1, idx2))), no.init.codon = TRUE), 
      idx1, idx2, strand = "-", DNAseq = substr(dna, idx1, idx2)
    ))
  }
}

# import all peptides presented on class 1 / 2
all_peptides <- fread("../../epitope-scoring/output/C12_allPeptides.tsv") # %>% filter(known)

# Filter for peptides with unique mapping
all_peptides_ss <- all_peptides %>% filter(Peptide != "APGPGGGAAVPSGAT") %>% # 12 hits
  filter(!(Peptide %in% c("QAFLQGVKDSEDASR", "SEEPLAQAF", "RPRPPARSL"))) # 3 hits

coordinates_df <- lapply(all_peptides_ss$Peptide, parse_coordinates_peptide) %>% rbindlist() %>%
  mutate(class = ifelse(nchar(query) == 15, "classII", "classI"))
classII_df <- coordinates_df %>% filter(class == "classII")
classI_df <- coordinates_df %>% filter(class == "classI")

class1coords <- unique(unlist(Map(`:`, classI_df$idx1, classI_df$idx2))) %>% sort()
class2coords <- unique(unlist(Map(`:`, classII_df$idx1, classII_df$idx2)))%>% sort()

ukb_dft %>%
  filter(TOTAL > 100) %>%
  filter(Consequence %in% c("missense_variant")) %>%
  mutate(in_class1 = Location %in% class2coords) %>%
  ggplot(aes(x = value + 0.001, color = in_class1)) +
  stat_ecdf() +
  scale_x_log10()

ukb_dft %>%
  filter(Consequence %in% c("missense_variant", "synonymous_variant")) %>%
  filter(Location %in% class2coords) %>%
  arrange(desc(value)) %>%
  filter(Consequence == "missense_variant") %>%
  mutate(rank = 1:n()) %>%
  filter(value > 0) %>% 
  ggplot(aes(x = rank, y = value, label = BP, color = Feature)) + 
  geom_text()


