library(Seurat)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

setwd('/work/gr-aztekin/3.project/culture_AER_final/Fig1')
source("../../axolotl_AER_final/0.scripts/dotplot_fromTOM.r")
source("../../axolotl_AER_final/0.scripts/cell_type_composition.r")

combined <- readRDS("mESC_2D.rds")
DefaultAssay(combined) <- "RNA"

## cluster color
# col.lb_pro <- '#00BF93'
# col.se <- '#FFC58E'
# col.aer <- '#8E0068'

col.ct <- c(
    "#dcc3fd", "#9f96f1",
    "#fbcea3", "#f99e4a", "#db995a", "#ad733c",
    "#8E0068", "#94c5e5",
    "#94b44f", "#2176ff",
    "grey"
)
names(col.ct) <- levels(combined$celltype)


pdf("clustering_celltype.pdf")
DimPlot(combined,
    cols = col.ct,
    group.by = "celltype", label = F, label.box = F
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

dev.off()


pdf("clustering_celltype_day7.pdf")
DimPlot(combined,
    cells = WhichCells(combined, expression = orig.ident == "day7"),
    cols = col.ct,
    group.by = "celltype", label = F, label.box = F
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

dev.off()



## time and batch
pdf("clustering_time_batch.pdf")

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
    plot = T, return = T
)

pdf("marker_all2.pdf", width = 10, height = 4)
p1 + coord_flip()
dev.off()


# selected: AER, Surface ectoderm and mesoderm in day7
genes.to.plot <- c(
    "Epcam", "Cdh1", "Wnt6", "Tfap2c,Tfap2a", "Krt8", "Krt18", "Krt5", "Krt14",
    "Trp63", "Fgf8", "Msx2", "Wnt5a", "Bmp4", "Sp8", "Sp6", "Dlx5",
    "Prrx1", "Meis1", "Meis2", "Hoxb5", "Hoxb6", "Irx3", "Irx5",
    "Etv3", "Etv4", "Etv5", "Msx1", "Ptch1", "Gli3", "Hand1", "Hand2",
    "Fgf10", "Shh"
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
    aspect.ratio = 4, ySize = 7, xSize = 7, dot.scale = 4, plot = T, return = T
)

pdf("marker_selected2.pdf", width = 3, height = 6)
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

pdf("cell_type_composition.pdf", width = 7, height = 3)
p3
dev.off()
# save as table
prop.table(table(combined$time_batch, combined$celltype), margin = 1) %>% write.csv('cell_type_composition.csv', quote = F)


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
p3.5
ggsave("cell_type_composition_day7_simplified.pdf", width = 4, height = 4)
# save as table
prop.table(table(df$time_batch, df$need), margin = 1) %>% write.csv('cell_type_composition_day7_simplified.csv', quote = F)



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

pdf("cell_cycle_composition.pdf", width = 7, height = 3)
p4
dev.off()
# save as table
prop.table(table(combined$time_batch, combined$Phase), margin = 1) %>% write.csv('cell_cycle_composition.csv', quote = F)



## Mesodermal core specific genes
genes <- c("Prrx1", "Lef1", "Tfap2c", "Trp63")
p5 <- VlnPlot(subset(combined, celltype != "Unknown"),
    features = genes,
    cols = col.ct,
    pt.size = 0, ncol = 4
) & theme(axis.title.x = element_blank())

pdf("meso_gene_selected.pdf", width = 16, height = 5)
p5
dev.off()
