library(Biostrings)
library(dplyr)
library(data.table)
library(BuenColors)

# Look for it in AOU
aou_df <- fread("../data/AOU_chrEBV.tsv.gz")
ukb_df <- fread("../data/UKB_chrEBV.tsv.gz")

aou_df$REF <- factor(aou_df$REF, levels = c("A", "C", "G", "T"))
ukb_df$REF <- factor(ukb_df$REF, levels = c("A", "C", "G", "T"))

#https://stackoverflow.com/questions/20036255/get-the-vector-of-values-from-different-columns-of-a-matrix
count_mat <- matrix(data = c(aou_df[["A"]], aou_df[["C"]], aou_df[["G"]], aou_df[["T"]]),ncol = 4)
aou_df$ref_count <- count_mat[cbind(1:dim(aou_df)[1],as.numeric(aou_df$REF))]

count_mat2 <- matrix(data = c(ukb_df[["A"]], ukb_df[["C"]], ukb_df[["G"]], ukb_df[["T"]]),ncol = 4)
ukb_df$ref_count <- count_mat2[cbind(1:dim(ukb_df)[1],as.numeric(ukb_df$REF))]

mdf <- merge(
  (aou_df %>% mutate(minor_af_aou = 1- ref_count/TOTAL))[,c("BP", "minor_af_aou", "TOTAL")],
  (ukb_df %>% mutate(minor_af_ukb = 1- ref_count/TOTAL))[,c("BP", "minor_af_ukb", "TOTAL")], by = "BP"
) 

well_covered_df <- mdf %>%
  filter(TOTAL.x > 100 & TOTAL.y > 100)

p1 <- well_covered_df %>%
  filter(minor_af_aou > 0 | minor_af_ukb > 0) %>%
  ggplot(aes(x = minor_af_aou*100, y = minor_af_ukb*100)) +
  geom_point() +
  pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "non-reference allele frequency, AoU", y = "non-reference allele frequency, UKB") 
cowplot::ggsave2(p1, file = "../plots/EBV_af_af.pdf", width = 2, height = 2)

cowplot::ggsave2(p1+ theme_void(), file = "../plots/EBV_af_af.png", width = 2*2, height = 2*2, dpi = 400)

cor(well_covered_df$minor_af_aou, well_covered_df$minor_af_ukb)

mdf %>% 
  filter(TOTAL.x > 100 & TOTAL.y > 100) %>%
  filter(abs(minor_af_aou-minor_af_ukb) > 0.3)

sum(well_covered_df$minor_af_aou > 0.01)
# 5471
 sum(well_covered_df$minor_af_ukb > 0.01)
# 13568
