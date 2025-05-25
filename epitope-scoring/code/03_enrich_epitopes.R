library(Biostrings)
library(data.table)
library(dplyr)
library(BuenColors)

# Read in fasta file
fasta_file <- "../data/EBV_proteins_genename.fasta" 
proteins <- readAAStringSet(fasta_file)

# annotations
gene_anno <- fread("../data/EBV_annotation.tsv", header = FALSE)
gene_anno_vec <- gene_anno[["V2"]];  names(gene_anno_vec) <- gsub("[.]", "_", gsub("/", "_", gene_anno[["V1"]]))

all_class1 <- rbind(
  fread("../intermediate/EBV_proteins_unique_8mers.csv"),
  fread("../intermediate/EBV_proteins_unique_9mers.csv"),
  fread("../intermediate/EBV_proteins_unique_10mers.csv"),
  fread("../intermediate/EBV_proteins_unique_11mers.csv")
) %>% mutate(what = unname(gene_anno_vec[ID]))

all_class2 <- rbind(
  fread("../intermediate/EBV_proteins_unique_15mers.csv")
)%>% mutate(what = unname(gene_anno_vec[ID]))

# see what is best ranked
br_class1 <- fread("../intermediate/peptides.best.rank.alleles.csv") %>% pull(Peptide)
br_class2 <- rbind(
  fread("../intermediate/best_peptide_DPAB.csv"),
  fread("../intermediate/best_peptide_DQAB.csv"),
  fread("../intermediate/best_peptide_DRB.csv")
  )%>% pull(Peptide)

all_class1$is_br <- all_class1$Peptide %in% br_class1
all_class2$is_br <- all_class2$Peptide %in% br_class2

class1_summary <- all_class1 %>%
  group_by(what) %>%
  summarize(total_possible = n(), total_br = sum(is_br)) %>%
  mutate(prop_BR = total_br/sum(total_br), prop_total = total_possible/sum(total_possible)) %>%
  mutate(log2_fold_change = log2(prop_BR/prop_total))

class2_summary <- all_class2 %>%
  group_by(what) %>%
  summarize(total_possible = n(), total_br = sum(is_br)) %>%
  mutate(prop_BR = total_br/sum(total_br), prop_total = total_possible/sum(total_possible)) %>%
  mutate(log2_fold_change = log2(prop_BR/prop_total))


class1_summary$pvalue <- sapply(1:4, function(i){
  prop.test(class1_summary$total_br[i], sum(class1_summary$total_br), p = class1_summary$prop_total[i])$p.value
  })

class2_summary$pvalue <- sapply(1:4, function(i){
  prop.test(class2_summary$total_br[i], sum(class2_summary$total_br), p = class2_summary$prop_total[i])$p.value
})

rbind(
  class1_summary %>% mutate(class = "I"),
  class2_summary %>% mutate(class = "II")
) %>%
  ggplot(aes(x = log2_fold_change, y = -log10(pvalue), color = what, shape = class)) + 
  geom_point() + geom_hline(yintercept = -log10(0.05/4), linetype = 2) +
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = jdb_palette("brewer_spectra")[c(9,3,7,1)]) +
  scale_x_continuous(limits = c(-2, 2)) -> pv
cowplot::ggsave2(pv, file = "../plots/volcano_what_genes.pdf", width = 1.4, height = 1)
