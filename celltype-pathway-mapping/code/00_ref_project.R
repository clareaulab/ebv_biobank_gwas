library(Seurat) # should be V4
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)
library(SeuratDisk)
library(viridis)
library(annotables)

# Import reference
## NOTE: The human PBMC dataset can be downloaded here: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
## Additional instructions on loading the dataset can be found here: https://satijalab.org/seurat/archive/v4.3/multimodal_reference_mapping
ref_path <- "~/Downloads/pbmc_multimodal.h5seurat" # change path according to where you download the Seurat object form the URL above

options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat(ref_path)

# Load hit genes
ebv_hits <- (fread("../data/ebv_hits.txt", header = FALSE)[[1]])
ss <- (rownames(reference@assays$SCT))[(rownames(reference@assays$SCT) %in% ebv_hits) & !grepl("HLA", rownames(reference@assays$SCT))]
ss <- ss[ss != "RPS18"]
reference <- AddModuleScore(reference, features = list(ebv_hits, ss), name = "EBV_hits")

pu <- FeaturePlot(reference, c("EBV_hits1"), raster = FALSE, reduction = "wnn.umap", max.cutoff = "q95", min.cutoff = "q05") &
  scale_color_viridis() & theme_void() & theme(legend.position = "none") 
cowplot::ggsave2(pu, file = "../plots/umap_ebv_score.png", width = 6, height = 6.1, dpi = 400)

pu2 <- FeaturePlot(reference, c("EBV_hits1"), raster = FALSE, reduction = "wnn.umap", max.cutoff = "q95", min.cutoff = "q05") &
  scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)])) & 
  theme_void() & theme(legend.position = "none") 
cowplot::ggsave2(pu2, file = "../plots/umap_ebv_score_reds.png", width = 6, height = 6.1, dpi = 400)



px <- ggplot(reference@meta.data, aes(x = celltype.l2, y = EBV_hits1)) +
  geom_violin(aes(reorder(celltype.l2, EBV_hits1, mean, decreasing = TRUE))) + 
  geom_boxplot(outlier.shape = NA, fill = NA, color = "firebrick", width = 0.5) +
  pretty_plot(fontsize = 8) + labs(x = "Annotated cell type", y = "EBV hit module score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  L_border()
cowplot::ggsave2(px, file = "../plots/box_ebv_azimuth_L2.pdf", width = 5, height = 2.5)


pb <- ggplot(reference@meta.data %>% filter(!(celltype.l1 %in% c("other", "other T"))), aes(x = celltype.l1, y = EBV_hits1)) +
  geom_violin() + 
  geom_boxplot(outlier.shape = NA, fill = NA, color = "firebrick", width = 0.5) +
  pretty_plot(fontsize = 8) + labs(x = "Annotated cell type", y = "EBV hit module score") +
  scale_x_discrete(limits=rev) +
  L_border() + coord_flip()
cowplot::ggsave2(pb, file = "../plots/box_ebv_module.pdf", width = 1.8, height = 2.2)

pu2 <- DimPlot(reference, reduction = "wnn.umap", group.by = "celltype.l2", raster = FALSE) +
  scale_color_manual(values = (c(jdb_palette("corona"), "salmon"))) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(pu2, file = "../plots/seurat_baseline.png", width = 6, height = 6.1, dpi = 400)


pl1 <- DimPlot(reference, reduction = "wnn.umap", group.by = "celltype.l1", raster = FALSE) +
  scale_color_manual(values = (c(jdb_palette("corona")[c(1:3,6,4,5)],"lightgrey", "lightgrey"))) +
  theme_void() + theme(legend.position = "none")

cowplot::ggsave2(pl1, file = "../plots/seurat_baseline_L1.png", width = 6, height = 6.1, dpi = 400)
