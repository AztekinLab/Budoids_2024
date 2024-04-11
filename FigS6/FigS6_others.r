library(Seurat)
library(ggplot2)
library(patchwork)

setwd("/work/gr-aztekin/3.project/culture_AER_final/FigS6")
combined <- readRDS("../Fig1/mESC_2D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"

# subset
clusters <- c("Surface Ectoderm 1", "Surface Ectoderm 2", "AER-like")
ect_only <- subset(combined, cells = WhichCells(combined, idents = clusters, expression = orig.ident == "day7"))
mes_only <- subset(combined, cells = WhichCells(combined, idents = "Mesodermal", expression = orig.ident == "day7"))
ect_only$celltype <- as.character(ect_only$celltype)
mes_only$celltype <- as.character(mes_only$celltype)

markers <- readxl::read_xlsx("../marker.xlsx")
markers <- as.data.frame(markers)


col.ct <- c(
    "#dcc3fd", "#9f96f1",
    "#fbcea3", "#f99e4a", "#db995a", "#ad733c",
    "#8E0068", "#00BF93",
    "#94b44f", "#2176ff",
    "grey"
)
names(col.ct) <- levels(combined$celltype)


ect <- c("Pre-placodal_ectoderm", "Trophectoderm", "Otic_cell", "Neural_crest")
mes <- c("Mesothelium", "Cardiac_mesoderm", "Cranial_mesoderm")
marker.NOT.ect <- unique(na.omit(unlist(c(markers[, ect]))))
marker.NOT.mes <- unique(na.omit(unlist(c(markers[, mes]))))


p1 <- DoHeatmap(ect_only,
    features = rev(marker.NOT.ect),
    group.by = "celltype", assay = "RNA", slot = "data",
    group.colors = col.ct[clusters],
    disp.max = 3, angle = 0, hjust = 0.5
) & guides(colour = "none") &
    scale_fill_viridis_c(option = "E", limits = c(0, 3))


p2 <- DoHeatmap(mes_only,
    features = rev(marker.NOT.mes),
    group.by = "celltype", assay = "RNA", slot = "data",
    group.colors = col.ct["Mesodermal"],
    disp.max = 3, angle = 0, hjust = 0.5
) & guides(colour = "none") &
    scale_fill_viridis_c(option = "E", limits = c(0, 3)) &
    theme(plot.background = element_blank())


pdf("annotion_marker_NOT_heatmap.pdf")
wrap_plots(p1, p2, ncol = 1, guides = "collect", heights = c(2, 1))
dev.off()
