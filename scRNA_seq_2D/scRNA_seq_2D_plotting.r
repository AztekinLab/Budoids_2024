library(Seurat)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

setwd('scRNA_seq_2D')
source("../0.scripts/dotplot_fromTOM.r")
source("../0.scripts/cell_type_composition.r")

combined <- readRDS("mESC_2D.rds")
DefaultAssay(combined) <- "RNA"

col.ct <- c(
    "#dcc3fd", "#9f96f1",
    "#fbcea3", "#f99e4a", "#db995a", "#ad733c",
    "#8E0068", "#94c5e5",
    "#94b44f", "#2176ff",
    "grey"
)
names(col.ct) <- levels(combined$celltype)


pdf("FigS2A_clustering_celltype.pdf")
DimPlot(combined,
    cols = col.ct,
    group.by = "celltype", label = F, label.box = F
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

dev.off()


pdf("Fig1G_clustering_celltype_day7.pdf")
DimPlot(combined,
    cells = WhichCells(combined, expression = orig.ident == "day7"),
    cols = col.ct,
    group.by = "celltype", label = F, label.box = F
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

dev.off()



## time and batch
pdf("FigS2BC_clustering_time_batch.pdf")

DimPlot(combined,
    cols = c("#bad9e7", "#7bb3cd", "#548ba5", "#46606b"),
    group.by = "orig.ident", label = F, label.box = F
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())


DimPlot(combined,
    cols = c("#E41A1C", "#377EB8"),
    group.by = "batch", label = F, label.box = F
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

dev.off()


## marker expression
markers <- readxl::read_xlsx("../marker.xlsx")
markers <- as.data.frame(markers)

# all
cols <- c("Stem_cell", "Surface_ectoderm", "AER", "Limb_mesoderm", "Distal_limb_mesoderm", "FL_HL", "Neural_ectoderm", "Neuron", "Endoderm")
genes.to.plot <- unlist(lapply(cols, function(x) markers[, x]))
genes.to.plot <- na.omit(unique(genes.to.plot))

p1 <- dotplot(combined@assays$RNA@data,
    genes = rev(c(genes.to.plot, "Ctnna2", "Crabp1", "Crabp2")),
    norm = "max", condition = combined$celltype, title = "Cell type markers",
    collow = "lightgrey", colhigh = "red4",
    aspect.ratio = 0.25, ySize = 7, xSize = 7, dot.scale = 3,
    plot = F, return = T
)

pdf("FigS2D_marker_all.pdf", width = 10, height = 4)
p1 + coord_flip()
dev.off()


# selected: AER, Surface ectoderm and mesoderm in day7
genes.to.plot <- c(
    "Epcam", "Cdh1", "Wnt6", "Tfap2c,Tfap2a", "Krt8", "Krt18", "Krt5", "Krt14",
    "Trp63", "Fgf8", "Msx2", "Wnt5a", "Bmp4", "Sp8", "Sp6", "Dlx5",
    "Prrx1", "Meis1", "Meis2", "Hoxb5", "Hoxb6", "Irx3", "Irx5",
    "Etv3", "Etv4", "Etv5", "Msx1", "Ptch1", "Gli3", "Hand1", "Hand2",
    "Fgf10", "Shh", 'Grem1', 'Tbx5',' Tbx4', 'Pitx1'
)
clusters <- c(
    "Surface Ectoderm 1", "Surface Ectoderm 2",
    "AER-like", "Mesodermal"
)

Idents(combined) <- "celltype"
p2 <- dotplot(combined@assays$RNA@data,
    genes = genes.to.plot,
    cells = WhichCells(combined, idents = clusters, expression = orig.ident == "day7"),
    norm = "max", condition = combined$celltype, title = "Selected cell types",
    collow = "lightgrey", colhigh = "red4",
    aspect.ratio = 4, ySize = 7, xSize = 7, dot.scale = 4, plot = F, return = T
)

pdf("Fig1J_marker_selected.pdf", width = 3, height = 6)
p2
dev.off()




# perc of cluster across time & batch
combined$time_batch <- paste(combined$orig.ident, combined$batch, sep = "_")

p3 <- composition_plot(
    df = combined@meta.data, group = time_batch,
    fill = celltype, colours = col.ct
) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) &
    labs(title = "Cell type compositions across timepoints and batches")

pdf("FigS2E_cell_type_composition.pdf", width = 7, height = 3)
p3
dev.off()
# save as table
prop.table(table(combined$time_batch, combined$celltype), margin = 1) %>% write.csv('FigS2E_cell_type_composition.csv', quote = F)


df <- combined@meta.data[combined$orig.ident == "day7", ]
df$need <- "Ectoderm"
df$need[df$celltype %in% c("Naïve mESC", "Primed mESC", "Neuroectoderm", "Sox17+ cells")] <- "Non-targeted"
df$need[df$celltype == "Mesodermal"] <- "Mesoderm"
df$need <- factor(df$need, levels = c("Ectoderm", "Mesoderm", "Non-targeted"))


p3.5 <- composition_plot(
    df = df,
    group = time_batch,
    fill = need, colours = c("#ffc58e", "#9ABACC", "grey")
) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) &
    labs(title = "Cell type compositions in day7")
ggsave("Fig1H_cell_type_composition_day7_simplified.pdf", width = 4, height = 4)
# save as table
prop.table(table(df$time_batch, df$need), margin = 1) %>% write.csv('FigS1H_cell_type_composition_day7_simplified.csv', quote = F)



p4 <- composition_plot(
    df = combined@meta.data, group = celltype,
    fill = Phase, colours = c("#808080", "#4daf4a", "#d9e021")
) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) &
    labs(title = "Cell cycle phases across cell types")

pdf("FigS2F_cell_cycle_composition.pdf", width = 7, height = 3)
p4
dev.off()
# save as table
prop.table(table(combined$celltype, combined$Phase), margin = 1) %>% write.csv('FigS2F_cell_cycle_composition.csv', quote = F)


# EMT D5 analysis
library(ComplexHeatmap)
EMT <- c("Cdh1", "Cdh2", "Epcam", "Zeb1", "Zeb2", "Vim", "Snai1", "Snai2", "Twist1")
seu_d5 <- subset(combined, orig.ident == "day5" & celltype %in% c("Surface Ectoderm 3", "Surface Ectoderm 4"))

pheat <- Heatmap(as.matrix(seu_d5@assays$RNA@data[EMT, ]),
    col = circlize::colorRamp2(c(0, 1, 2, 3), viridis::cividis(4)), # cividis(100),
    cluster_columns = T,
    cluster_rows = T,
    show_row_dend = F,
    show_column_dend = F,
    show_column_names = F,
    column_split = seu_d5$celltype,
    # column_title = "EMT signature in Day5 cells",
    column_title = c("Surface Ectoderm 4", "Surface Ectoderm 3"),
    use_raster = F
)

pdf("FigS2G_EMT_heatmap.pdf", width = 8, height = 4)
pheat
dev.off()


pdf("FigS2H_EMT_umap.pdf", width = 12, height = 12)
FeaturePlot(seu_d5, EMT,
    ncol = 3,
    cols = c("lightgrey", "red4"),
    order = T, raster = T, pt.size = 2
)
dev.off()



# checking non-limb ectoderm/mesodermal markers

Idents(combined) <- "celltype"

# subset
clusters <- c("Surface Ectoderm 1", "Surface Ectoderm 2", "AER-like")
ect_only <- subset(combined, cells = WhichCells(combined, idents = clusters, expression = orig.ident == "day7"))
mes_only <- subset(combined, cells = WhichCells(combined, idents = "Mesodermal", expression = orig.ident == "day7"))
ect_only$celltype <- as.character(ect_only$celltype)
mes_only$celltype <- as.character(mes_only$celltype)


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


pdf("FigS3_CE_NOTmarker_heatmap.pdf")
wrap_plots(p1, p2, ncol = 1, guides = "collect", heights = c(2, 1))
dev.off()



## Mesodermal core specific genes
genes <- c("Prrx1", "Lef1", "Tfap2c", "Trp63")
p5 <- VlnPlot(subset(combined, celltype %in% c("Surface Ectoderm 1", "Surface Ectoderm 2", "AER-like", "Mesodermal") & orig.ident == 'day7'),
    features = genes,
    cols = col.ct,
    pt.size = 0, ncol = 4
) & theme(axis.title.x = element_blank())

pdf("FigS4D_meso_gene_selected.pdf", width = 16, height = 5)
p5
dev.off()


pdf("FigS9A_D7only_epcam_cd9_cd44_ridge.pdf", width = 10, height = 3)
RidgePlot(subset(combined, orig.ident == "day7"), c("Epcam", "Cd9"),
    idents = c("Surface Ectoderm 1", "Surface Ectoderm 2", "AER-like", "Mesodermal"), # main populations, % > 5%
    same.y.lims = T, cols = col.ct
)
dev.off()
