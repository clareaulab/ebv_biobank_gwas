library(data.table)
library(BuenColors)
library(dplyr)

dt <- fread("../data/NPC_VUS.csv")
vus <- dt[!is.na(dt$variant_nt),]
ordered_dt <- vus %>%
  mutate(MAFavg = (MAF_aou + MAF_ukb)/2) %>%
  arrange(desc(MAFavg)) %>% mutate(rank = 1:n())
reshape2::melt(ordered_dt[,c("gene", "variant_aa", "variant_nt","rank","MAFavg", "MAF_aou", "MAF_ukb")], 
              id.vars = c("gene", "variant_aa", "variant_nt", "rank","MAFavg")) -> in_df
  
in_df %>% 
  ggplot(aes(x = rank, color = variable)) +
  geom_bar(data = in_df %>% filter(variable == "MAF_aou"), stat = "identity", aes(y = MAFavg*100), fill = "lightgrey", color = "black")+
  geom_point(aes(y = value*100)) +
  pretty_plot(fontsize = 8) + L_border() + scale_color_manual(values = c('firebrick', "dodgerblue3")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "NPC VUS ranked by mean allele frequency", y = "Biobank allele frequency") +
  theme(legend.position = "none") -> p1
p1
cowplot::ggsave2(p1, file = "../plots/VUS_waterfall.pdf", width = 2.8, height = 2)

write.table(in_df, file = "../output/VUS_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
