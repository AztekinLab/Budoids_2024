library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrastr)
source("../../axolotl_AER_final/0.scripts/cell_type_composition.r")
source("../../axolotl_AER_final/0.scripts/dotplot_fromTOM.r")

setwd("/work/gr-aztekin/3.project/culture_AER_final/Fig3")
combined <- readRDS("mESC_3D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"


meso <- c("#94c5e5", "#6fb6e5", "#0087e1", "#078eba", "#046484") # green
ect <- c("#fbcea3", "#f99e4a", "#ad420d") # orange
aer <- "#8E0068"

neu_cre <- "#cedb9c"
neu_ect <- "#94b44f"
neuron <- "#637939"
chp <- "#405020"

endo <- "#2176ff"
pluri <- "#9f96f1"
unknow <- "grey"

ct_col <- c(ect, aer, meso, neu_ect, neuron, chp, neu_cre, endo, pluri, unknow)
names(ct_col) <- levels(combined$celltype)
time_col <- c("#c2e1ef", "#68a5c2", "#46606b")
names(time_col) <- levels(combined$orig.ident)
batch_col <- c("#E41A1C", "#377EB8", "#4DAF4A")



# final clustering UMAP
DimPlot(combined, group.by = "celltype", cols = ct_col)
ggsave("clustering_celltype2.pdf")


# time and batch
set.seed(1234)
pdf("clustering_time.pdf", width = 9, height = 8)
# DimPlot(combined,
#     cols = time_col,
#     group.by = "orig.ident", label = F, label.box = F,
#     pt.size = 2, shuffle = TRUE,
#     raster = T, raster.dpi = c(1024, 1024)
# ) &
#     theme(legend.position = c(.05, .3)) & labs(title = element_blank())

DimPlot(combined,
    cols = time_col,
    group.by = "orig.ident", label = F, label.box = F,
    pt.size = 0.1, shuffle = TRUE
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

dev.off()


pdf("clustering_batch.pdf", width = 24, height = 8)
DimPlot(combined,
    cols = batch_col,
    group.by = "batch", split.by = "orig.ident", label = F, label.box = F,
    pt.size = 2, raster = T, raster.dpi = c(1024, 1024), shuffle = TRUE
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

DimPlot(combined,
    cols = batch_col,
    group.by = "batch", split.by = "orig.ident", label = F, label.box = F,
    pt.size = 0.1, shuffle = TRUE
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())
dev.off()


pdf("clustering_batch2.pdf", width = 9, height = 8)
DimPlot(combined,
    # cols = batch_col,
    group.by = "time_batch", label = F, label.box = F,
    pt.size = 0.1, shuffle = TRUE
) &
    theme(legend.position = c(.05, .3)) & labs(title = element_blank())

dev.off()



# percentage
combined$time_batch <- paste(combined$orig.ident, combined$batch, sep = "_")
combined$time_batch <- factor(combined$time_batch,
    levels = c("day7_Rep1", "day7_Rep2", "day9_Rep3", "day12_Rep2", "day12_Rep3")
)

p1 <- composition_plot(combined@meta.data,
    group = time_batch, fill = celltype,
    colours = ct_col
)

p1
ggsave("celltype_composition.pdf")


df <- combined@meta.data
df$need <- "Non-targeted"
df$need[df$celltype %in% c("Surface Ectoderm 1", "Surface Ectoderm 2", "BasalEctoderm-like", "AER-like")] <- "Ectoderm"
df$need[df$celltype %in% c("Mesodermal-like", "Mes_Hoxb", "Mes_Hox13", "Mes_Sox9", "Mes_Col2a1")] <- "Mesoderm"
df$need <- factor(df$need, levels = c("Ectoderm", "Mesoderm", "Non-targeted"))

df$need2 <- as.character(df$need)
df$need2[df$celltype %in% c("Mesodermal-like", "Mes_Hoxb", "Mes_Hox13")] <- "Mesoderm(soft)"
df$need2[df$celltype %in% c("Mes_Sox9", "Mes_Col2a1")] <- "Mesoderm(hard)"
df$need2 <- factor(df$need2, levels = c("Ectoderm", "Mesoderm(soft)", "Mesoderm(hard)", "Non-targeted"))


p2 <- composition_plot(
    df = df,
    group = time_batch,
    fill = need, colours = c("#ffdaa9", "#6ba96b", "#757680")
) &
    theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) &
    labs(title = "Cell type compositions")
ggsave("celltype_composition_simplified.pdf", width = 6, height = 6)



p3 <- composition_plot(
    df = df,
    group = time_batch,
    fill = need2, colours = c("#ffdaa9", "#4f994f", "#046484", "#757680")
) &
    theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) &
    labs(title = "Cell type compositions")
ggsave("celltype_composition_simplified2.pdf", width = 6, height = 6)


# freq <- table(combined$orig.ident, combined$celltype)
# perc <- round(100 * freq / rowSums(freq), 2)
# write.table(perc, "percent_celltype.txt", sep = "\t", quote = F)



# cell cycle percentage
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

pdf("cell_cycle_composition.pdf", width = 9, height = 3)
p4
dev.off()



## marker expression
markers <- readxl::read_xlsx("../marker.xlsx")
markers <- as.data.frame(markers)

cols <- c(
    "Surface_ectoderm", "Basasl_ectoderm", "AER", "Limb_mesoderm", "Distal_limb_mesoderm",
    "Fibro", "Tenocyte", "Pericyte", "Chondrocyte",
    "Neural_ectoderm", "Neuron", "Choroid_plexus", "Neural_crest",
    "Endoderm", "Stem_cell"
)
genes.to.plot <- unlist(lapply(cols, function(x) markers[, x]))
genes.to.plot <- na.omit(unique(genes.to.plot))

p1 <- dotplot(combined@assays$RNA@data,
    genes = c(genes.to.plot, "Anxa4", "Irs3"), # add DEG for "unknown"
    norm = "max", condition = combined$celltype, title = "Cell type markers",
    collow = "lightgrey", colhigh = "red4",
    aspect.ratio = 4, ySize = 7, xSize = 7, dot.scale = 3,
    plot = T, return = T
)

pdf("marker_all.pdf", width = 4, height = 10)
p1
dev.off()

# selected genes
# Msx1, Grem1, Fgf10
# Hox
library(ggrastr)

mes <- c("Prrx1", "Msx1", "Grem1", "Fgf10")
# hox <- c("Hoxa9", "Hoxa11", "Hoxa13", "Hoxd9", "Hoxd11", "Hoxd13")
hox <- c("Hoxa9", "Hoxa11os", "Hoxa13", "Hoxd9", "Hoxd11", "Hoxd13")
carti <- c("Sox9", "Col2a1", "Col9a1", "Sox6")
ect <- c("Epcam", "Krt8", "Krt18", "Krt5", "Krt14", "Trp63")
diffect <- c("Krt10", "Krtdap", "Zfp750", "Lor", "Ivl") # K10 (early diff), filaggrin, involucrin and loricrin (late diff)
aer <- c("Fgf8", "Sp6", "Dlx5", "Msx1")

fp_col <- c("grey", "red4")

p1 <- FeaturePlot(combined,
    features = ect, order = T, ncol = 3,
    pt.size = 0.1, combine = F, cols = fp_col,
    cells = WhichCells(combined, idents = c("Surface Ectoderm 1", "Surface Ectoderm 2", "BasalEctoderm-like", "AER-like"))
)
p1 <- lapply(p1, function(p) {
    rasterize(p, layers = "Point", dpi = 600)
})


p1.1 <- FeaturePlot(combined,
    features = diffect, order = T, ncol = 3,
    pt.size = 0.1, combine = F, cols = fp_col,
    cells = WhichCells(combined, idents = c("Surface Ectoderm 1", "Surface Ectoderm 2", "BasalEctoderm-like", "AER-like"))
)
p1.1 <- lapply(p1.1, function(p) {
    rasterize(p, layers = "Point", dpi = 600)
})



p2 <- FeaturePlot(combined,
    features = aer, order = T, ncol = 3,
    pt.size = 0.1, combine = F, cols = fp_col,
    cells = WhichCells(combined, idents = c("Surface Ectoderm 1", "Surface Ectoderm 2", "BasalEctoderm-like", "AER-like"))
)
p2 <- lapply(p2, function(p) {
    rasterize(p, layers = "Point", dpi = 600)
})

p3 <- FeaturePlot(combined,
    features = mes, order = T, ncol = 3,
    pt.size = 0.1, combine = F, cols = fp_col,
    cells = WhichCells(combined, idents = c("Mesodermal-like", "Mes_Hoxb", "Mes_Hox13"))
)
p3 <- lapply(p3, function(p) {
    rasterize(p, layers = "Point", dpi = 600)
})

p4 <- FeaturePlot(combined,
    features = hox, order = T, ncol = 3,
    pt.size = 0.1, combine = F, cols = fp_col,
    cells = WhichCells(combined, idents = c("Mes_Hoxb", "Mes_Hox13"))
)
p4 <- lapply(p4, function(p) {
    rasterize(p, layers = "Point", dpi = 600)
})

p5 <- FeaturePlot(combined,
    features = carti, order = T, ncol = 3,
    pt.size = 0.1, combine = F, cols = fp_col,
    cells = WhichCells(combined, idents = c("Mes_Sox9", "Mes_Col2a1"))
)
p5 <- lapply(p5, function(p) {
    rasterize(p, layers = "Point", dpi = 600)
})



pdf("marker_selected.pdf", width = 20, height = 12)
wrap_plots(p1, ncol = 3)
wrap_plots(p1.1, ncol = 3)
wrap_plots(p2, ncol = 3)
wrap_plots(p3, ncol = 3)
wrap_plots(p4, ncol = 3)
wrap_plots(p5, ncol = 3)
dev.off()


# hox sequential activation
pdf("hox_exp_sequential.pdf", width = 8, height = 6)
# VlnPlot(subset(combined, subset = celltype %in% c("Mes_Hoxb")), features = hox, split.by = "orig.ident", pt.size = 0.1, col = time_col) +
#     theme(legend.position = "right")

VlnPlot(subset(combined, subset = celltype %in% c("Mes_Hox13")), features = hox, split.by = "orig.ident", pt.size = 0.1, col = time_col) +
    theme(legend.position = "right")
dev.off()





# FL/HL identity
## co-expression scatterplot
library(ggridges)

genes <- c("Tbx5", "Pitx1", "Tbx4")
mat <- combined@assays$RNA@data[genes, ]
df <- as.data.frame(t(as.matrix(mat)))
df$stage <- combined$orig.ident

p <- ggplot(df, aes(Tbx5, Pitx1, color = stage)) +
    scale_color_manual(values = time_col) +
    geom_point() +
    theme_bw() +
    guides(shape = "none", size = "none")

xdensity <- ggplot(df, aes(Tbx5, stage, fill = stage)) +
    # geom_density() +
    geom_density_ridges(scale = 2) +
    scale_fill_manual(values = time_col) +
    theme_classic() +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        # axis.title = element_blank(),
    )

ydensity <- ggplot(df, aes(Pitx1, stage, fill = stage)) +
    # geom_density() +
    geom_density_ridges(scale = 2) +
    scale_fill_manual(values = time_col) +
    theme_classic() +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        # axis.title = element_blank(),
    ) +
    coord_flip()


design <- "
AAAD
BBBC
BBBC
BBBC"
p_all <- wrap_plots(list(xdensity, rasterize(p, layers = "Point", dpi = 600), ydensity, guide_area()), design = design, guides = "collect")
ggsave("coexpr_pitx1.pdf", width = 4, height = 4)



p1 <- ggplot(df, aes(Tbx5, Tbx4, color = stage)) +
    scale_color_manual(values = time_col) +
    geom_point() +
    theme_bw() +
    guides(shape = "none", size = "none")

ydensity1 <- ggplot(df, aes(Tbx4, stage, fill = stage)) +
    # geom_density() +
    geom_density_ridges(scale = 2) +
    scale_fill_manual(values = time_col) +
    theme_classic() +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        # axis.title = element_blank(),
    ) +
    coord_flip()

p_all <- wrap_plots(list(xdensity, rasterize(p1, layers = "Point", dpi = 600), ydensity1, guide_area()), design = design, guides = "collect")
ggsave("coexpr_tbx4.pdf", width = 4, height = 4)
