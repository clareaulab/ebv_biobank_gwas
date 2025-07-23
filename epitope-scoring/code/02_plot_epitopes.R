library(BuenColors)
library(rtracklayer)
library(dplyr)
library(data.table)
library(seqinr)
library(stringr)
library(readxl)

# known epitopes 
known_EBV <- read_xlsx("../data/epitope_table_export_1745338611.xlsx")

### Compute statistics
class1_peptidome <- rbind(
  fread("../intermediate/EBV_proteins_unique_8mers.csv"),
  fread("../intermediate/EBV_proteins_unique_9mers.csv"),
  fread("../intermediate/EBV_proteins_unique_10mers.csv"),
  fread("../intermediate/EBV_proteins_unique_11mers.csv")
)

class2_peptidome <- rbind(
  fread("../intermediate/EBV_proteins_unique_15mers.csv")
)
class1_candidate_peptides <- unique(class1_peptidome$Peptide)
class2_candidate_peptides <- unique(class2_peptidome$Peptide)

observed_MHC1 <- intersect(known_EBV[[3]], class1_candidate_peptides)
observed_MHC2 <- intersect(known_EBV[[3]], class2_candidate_peptides)

# helper function
fix_gene_id <- function(gene){
   gsub("[.]", "_", gsub("/", "_", gene))
 }

# import protein fasta
fasta <- import("../data/EBV_proteins_genename.fasta", type = "AA")
protein_df <- data.frame(
  gene = fix_gene_id(names(fasta)),
  seq = unname(as.character(fasta))
)

# import hits
hit_df_c1 <- fread("../intermediate/peptides.best.rank.alleles.csv"); hit_df_c1 <- hit_df_c1[complete.cases(hit_df_c1),]
hit_df_c1 <- hit_df_c1[!duplicated(hit_df_c1),]
hit_df_c1$rel_position <- sapply(1:dim(hit_df_c1)[1], function(i){
  gene1 <- hit_df_c1[i,"ID"][[1]]
  peptide1 <- hit_df_c1[i,"Peptide"][[1]]
  seq1 <- protein_df %>% filter(gene == gene1) %>% pull(seq)
  ((str_locate(seq1,peptide1))[[1]]+5)/nchar(seq1)
})
hit_df_c1$HLA_class <- substr(hit_df_c1$allele, 1, 1)

# 9 of 83 predicted epitopes are known
intersect(hit_df_c1$Peptide, known_EBV[[3]]) %>% length()
length(unique(hit_df_c1$Peptide))

# Import coordinates
gff3 <- data.frame(import("../data/sequence.gff3")) %>% filter(type == "CDS") %>%
  mutate(gene = fix_gene_id(gene)) %>%
  group_by(gene) %>% slice_head(n=1)

full_merge_df <- merge(hit_df_c1, gff3[,c("gene","start", "end", "width")], by.x = "ID" ,by.y = "gene") 
  
full_merge_df %>%
  filter(Peptide == "ATIGTAMYK")

px <- full_merge_df %>%
  mutate(coord_midpoint = width*rel_position + start) %>%
  ggplot(aes(x = coord_midpoint, y = -log10(ELRank), color = HLA_class)) +
  geom_point() + scale_color_manual(values = jdb_palette("corona")) +
  scale_y_continuous(limits = c(1.5, 3.6)) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none")
cowplot::ggsave2(px, file = "../plots/rank_score_MHC1.pdf", width = 2.5, height = 1)


# Import class 2 
class2df <- rbind(
  fread("../intermediate/best_peptide_DPAB.csv") %>%
  mutate(HLA_class = "DP"),
  fread("../intermediate/best_peptide_DQAB.csv") %>%
    mutate(HLA_class = "DQ"),
  fread("../intermediate/best_peptide_DRB.csv") %>%
    mutate(HLA_class = "DR")
)
class2df <- class2df[!duplicated(class2df),]
colnames(class2df) <- c("allele", "Peptide", "ID", "ELRank",  "HLA_class")

intersect(class2df$Peptide, known_EBV[[3]]) %>% length()
length(unique(class2df$Peptide))


# class 1 binom test
binom.test(x = sum(hit_df_c1$Peptide %in% known_EBV[[3]]), n = dim(hit_df_c1)[1], p = (length(observed_MHC1)/length(class1_candidate_peptides))) %>% str()

# class 2 binom test
binom.test(x = sum(class2df$Peptide %in% known_EBV[[3]]), n = dim(class2df)[1], p = (length(observed_MHC2)/ length(class2_candidate_peptides))) %>% str()

#  7 / 106 total
class2df$rel_position <- sapply(1:dim(class2df)[1], function(i){
  gene1 <- class2df[i,"ID"][[1]]
  peptide1 <- class2df[i,"Peptide"][[1]]
  seq1 <- protein_df %>% filter(gene == gene1) %>% pull(seq)
  x <- ((str_locate(seq1[[length(seq1)]],peptide1))[[1]]+5)/nchar(seq1[[length(seq1)]])
  x
})

full_merge_df_class2 <- merge(class2df, gff3[,c("gene","start", "end", "width")], by.x = "ID" ,by.y = "gene") 
cl2 <- full_merge_df_class2 %>%
  mutate(coord_midpoint = width*rel_position + start) %>%
  ggplot(aes(x = coord_midpoint, y = -log10(ELRank), color = HLA_class)) +
  geom_point() + scale_color_manual(values = jdb_palette("corona")[c(4,5,6)]) +
  scale_y_continuous(limits = c(1, 3.6)) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none")
cowplot::ggsave2(cl2, file = "../plots/rank_score_MHC2.pdf", width = 2.5, height = 1)


total_c12 <- rbind(full_merge_df_class2 %>% mutate(class = "class2"), 
                   full_merge_df  %>% mutate(class = "class1"))
total_c12$known <- total_c12$Peptide %in% ( known_EBV[[3]])

write.table(total_c12, file = "../output/C12_allPeptides.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
redgrey <- total_c12 %>%
  arrange(known) %>% 
  mutate(coord_midpoint = width*rel_position + start) %>%
  ggplot(aes(x = coord_midpoint, y = -log10(ELRank), color = known, shape = class)) +
  geom_point() +
  scale_y_continuous(limits = c(1, 3.6)) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("lightgrey", "firebrick"))
cowplot::ggsave2(redgrey, file = "../plots/IE.pdf", width = 2.5, height = 1)

