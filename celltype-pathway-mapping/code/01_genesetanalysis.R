library(BuenColors)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(org.Hs.eg.db)
library(annotables)
library(viridis)

# Function to compute fold enrichments
calc_fold_enrichments <- function(a_df){
  new_df <- a_df %>% dplyr::select(GeneRatio, BgRatio) %>% 
    tidyr::separate(GeneRatio, into = c("GA","GB"), sep = "/", remove = TRUE, convert = TRUE) %>% 
    tidyr::separate(BgRatio, into = c("BA","BB"), sep = "/", remove = TRUE, convert = TRUE) 
  
  new_df <- new_df %>% mutate(gene_ratio = GA/GB, bg_ratio = BA/BB, fold_enrichment= gene_ratio/bg_ratio)
  a_df <- cbind(a_df,new_df)
  return(a_df)
}

# Options for cluster profiler
options(clusterProfiler.download.method = 'curl')
options(download.file.extra = '-k')
options(clusterProfiler.download.method = 'wget')
options(download.file.extra = '--no-check-certificate')

# Import enriched hits
non_chr6_hits <- readLines("../data/EBV_0.0018_Ch6rem.txt")
non_hla_hits <- readLines("../data/EBV_0.0018_Ch6inc_HLArem.txt")
all_hits <- readLines("../data/ebv_hits.txt")

# Define universe to do chromosome 6 enrichment
universe_df <- grch38 %>% filter(biotype == "protein_coding" & chr %in% c(1:22, "X"))
universe_df_noHLA <- universe_df[!grepl("HLA", universe_df$symbol),]
universe_df_no6 <- universe_df %>% filter(chr != 6)

# Run enrichment now
map_genes <- function(vec){
  bitr(vec, fromType="SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"), 
       OrgDb="org.Hs.eg.db")
}
enriched_genes_no6 <- map_genes(non_chr6_hits)
enriched_genes_noHLA <- map_genes(non_hla_hits)
enriched_genes <- map_genes(all_hits)

# Define background
background_no6 <- map_genes(universe_df_no6$symbol)
background_noHLA <- map_genes(universe_df_noHLA$symbol)
background_allgenes <- map_genes(universe_df$symbol)

# Do enrichment
goBP_enrich_no6 <- enrichGO(gene = enriched_genes_no6$ENSEMBL,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            universe = background_no6$ENSEMBL, 
                            qvalueCutoff  = 0.05) %>%  simplify(. , cutoff = 0.7) %>% as.data.frame() %>% calc_fold_enrichments()

goBP_enrich_noHLA <- enrichGO(gene = enriched_genes_noHLA$ENSEMBL,
                              OrgDb         = org.Hs.eg.db,
                              keyType       = 'ENSEMBL',
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              universe = background_noHLA$ENSEMBL, 
                              qvalueCutoff  = 0.05) %>%  simplify(. , cutoff = 0.7) %>% as.data.frame() %>% calc_fold_enrichments()

goBP_enrich_all <- enrichGO(gene = enriched_genes$ENSEMBL,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            universe = background_allgenes$ENSEMBL, 
                            qvalueCutoff  = 0.05)  %>%  simplify(. , cutoff = 0.7) %>% as.data.frame() %>% calc_fold_enrichments()

# Function to make enrichment plot automatically
make_enrichment_plot <- function(in_df){
  
  in_df %>% 
    arrange((zScore)) %>%
    mutate(rank = 1:n()) %>% 
    tail(5) -> plot_df 
  
  plot_df %>%
    ggplot(aes(x = rank, y = zScore, fill = -log10(qvalue))) + 
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) + 
    pretty_plot(fontsize = 8) + L_border() +
    scale_fill_gradientn(limits = c(0, max(-log10(in_df$qvalue))), oob = scales::squish, colors = jdb_palette("solar_rojos")) +
    theme(legend.position = "bottom") + labs( x= "", y = "enrichment z-score", fill = "")
}

p1 <- make_enrichment_plot(goBP_enrich_all)
p2 <- make_enrichment_plot(goBP_enrich_noHLA)
p3 <- make_enrichment_plot(goBP_enrich_no6)
cowplot::ggsave2(
  cowplot::plot_grid(p1, p2, p3, nrow = 1), 
  width = 6, height = 1.8, file = "../plots/enrich_BP.pdf"
)

head(goBP_enrich_noHLA[,!(colnames(goBP_enrich_noHLA) %in% c("geneID"))])
head(goBP_enrich_all[,!(colnames(goBP_enrich_all) %in% c("geneID"))])
head(goBP_enrich_no6[,!(colnames(goBP_enrich_no6) %in% c("geneID"))])

library(ggrepel)
p0 <- ggplot(goBP_enrich_all, aes(x = fold_enrichment, y = -log10(qvalue), label = Description)) +
  geom_point(size = 0.5) + geom_text_repel(max.overlaps = 5,  size = 1) +
  pretty_plot(fontsize = 8) + L_border()


kegg_enrich <- enrichKEGG(gene = enriched_genes$ENTREZID,
                            pAdjustMethod = "BH",
                            universe = background_allgenes$ENTREZID, 
                            qvalueCutoff  = 0.05) %>% as.data.frame() %>% calc_fold_enrichments()

p1 <- ggplot(kegg_enrich, aes(x = fold_enrichment, y = -log10(qvalue), label = Description)) +
  geom_point(size = 0.5) + geom_text_repel(max.overlaps = 5, size = 1) +
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(p0, file = "../plots/go_bp_volcano.pdf", width = 3.75, height = 2.5)
cowplot::ggsave2(p1, file = "../plots/kegg_volcano.pdf", width = 3.75, height = 2.5)
