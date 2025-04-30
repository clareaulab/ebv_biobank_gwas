library(data.table)
library(dplyr)
library(reshape2)
library(BuenColors)

aou_df <- fread("../data/AOU_chrEBV.tsv.gz")
ukb_df <- fread("../data/UKB_chrEBV.tsv.gz")

positions_t1t2 <- c(36209, 36226, 36251, 36252, 36258, 36275,  36302, 36312, 36320) 

process_data_frame <- function(df){
  reshape2::melt(df[,c("BP", "REF","TOTAL", "A","C","G", "T")], id.vars = c("BP", "REF", "TOTAL")) %>%
    mutate(var = value / TOTAL*100) %>%
    filter(var > 1 & var < 50) %>%
    mutate(variant = paste0("chrEBV:", BP, REF, ">", variable))
}

library(ggbeeswarm)
rbind(
  process_data_frame(aou_df %>% filter(BP %in% positions_t1t2)) %>% mutate(cohort = "AoU"),
  process_data_frame(ukb_df %>% filter(BP %in% positions_t1t2)) %>% mutate(cohort = "UKB")
) -> full_df

p1 <- full_df %>%
  ggplot(aes(x = cohort, y = 100-var, color = variant)) + 
  geom_bar(data = full_df %>% group_by(cohort) %>% summarize(var = mean(var)), fill = "lightgrey", color = "black", stat = "identity") +
  geom_quasirandom() +
  scale_color_manual(values = jdb_palette("corona")) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "cohort", y = "Estimated % Type 1 EBV genomes")
cowplot::ggsave2(p1, file = "../plots/EBV_bargraphs.pdf", width = 2.5, height = 2)

