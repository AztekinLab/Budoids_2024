library(Seurat)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

library(ComplexHeatmap)
library(circlize)
library(viridis)

setwd("/work/gr-aztekin/3.project/culture_AER_final/FigS4_and_S5")
combined <- readRDS("../Fig1/mESC_2D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"

# prepare in vitro dataset
aer_vitro <- subset(combined, cells = WhichCells(combined, idents = "AER-like", expression = orig.ident == "day7"))
mes_vitro <- subset(combined, cells = WhichCells(combined, idents = "Mesodermal", expression = orig.ident == "day7"))
aer_vitro$celltype <- as.character(aer_vitro$celltype)
mes_vitro$celltype <- as.character(mes_vitro$celltype)

# prepare in vivo dataset
fvivo <- Sys.glob("/work/gr-aztekin/3.project/culture_AER_final/1.data/processed/mouse*")
vivo <- lapply(fvivo, readRDS)

aer_vivo <- lapply(vivo, function(x) {
    x$celltype <- x$CellType2_sep
    subset(x, subset = CellType2_sep == "AER")
})
aer_vivo <- merge(aer_vivo[[1]], aer_vivo[-1])
Idents(aer_vivo) <- "orig.ident"


be_vivo <- lapply(vivo, function(x) {
    x$celltype <- x$CellType2_sep
    subset(x, subset = CellType2_sep == "BasalEctoderm")
})
be_vivo <- merge(be_vivo[[1]], be_vivo[-1])
Idents(be_vivo) <- "orig.ident"


mes_vivo <- lapply(vivo, function(x) {
    x$celltype <- x$CellType2_sep
    subset(x, subset = CellType2_sep == "CT")
})
mes_vivo <- merge(mes_vivo[[1]], mes_vivo[-1])
Idents(mes_vivo) <- "orig.ident"


muscle_vivo <- lapply(vivo, function(x) {
    x$celltype <- x$CellType2_sep
    subset(x, subset = CellType2_sep == "Muscle")
})
muscle_vivo <- merge(muscle_vivo[[1]], muscle_vivo[-1])
Idents(muscle_vivo) <- "orig.ident"


# prepare genes to plot
markers <- readxl::read_xlsx("../marker.xlsx")
markers <- as.data.frame(markers)
AER <- sort(c(
    "Fgf8", "Fgf4", "Fgf9", "Fgf17", "Fn1", "Msx2", "Sp6", "Sp8",
    "Dlx1", "Dlx2", "Dlx3", "Dlx4", "Dlx5", "Dlx6", "Gja1", "Lrp6", "Mdk",
    "Wnt5a", "Rspo2", "Bmp2", "Bmp4", "Bmp7", "Tgfb2", "Trp63", "Cd44", "Epcam"
))

mes <- unique(c(
    na.omit(markers$Limb_mesoderm), na.omit(markers$Distal_limb_mesoderm),
    "Tbx4", "Tbx5", "Pitx1", "Hoxd9", "Hoxa9"
))
mes <- mes[stringr::str_order(mes, numeric = T)]



# plot Mesodermal markers
set.seed(123)
cell_idx1 <- WhichCells(mes_vitro, downsample = 500)
mes_vitro_mat <- mes_vitro@assays$RNA@data[mes, cell_idx1]
mes_vitro_mat <- as.matrix(mes_vitro_mat)

cell_idx2 <- WhichCells(mes_vivo, downsample = 500)
mes_vivo_mat <- mes_vivo@assays$RNA@data[mes, cell_idx2]
mes_vivo_mat <- as.matrix(mes_vivo_mat)
mes_mat <- cbind(mes_vitro_mat, mes_vivo_mat)


cellid <- c(
    as.character(mes_vitro$orig.ident[cell_idx1]),
    as.character(mes_vivo$orig.ident[cell_idx2])
)
cellid <- factor(cellid, levels = c("day7", "Mouse_E95", "Mouse_E105", "Mouse_E115", "Mouse_E125"))

col.cond <- brewer.pal(5, "Accent")
names(col.cond) <- levels(cellid)

top_anno <- HeatmapAnnotation(
    Condition = cellid,
    col = list(Condition = col.cond)
)

heat <- Heatmap(mes_mat,
    col = colorRamp2(c(0, 0.5, 2, 3), cividis(4)), # cividis(100),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    row_order = mes,
    top_annotation = top_anno,
    column_split = cellid,
    column_title = "Mesoderm marker expression (Downsampling 500 cells)",
    use_raster = F
)

pdf("mes_markers.pdf", width = 10, height = 6)
heat
dev.off()

# plot AER markers
set.seed(123)
cell_idx1 <- WhichCells(aer_vitro, downsample = 200)
aer_vitro_mat <- aer_vitro@assays$RNA@data[AER, cell_idx1]
aer_vitro_mat <- as.matrix(aer_vitro_mat)

cell_idx2 <- WhichCells(aer_vivo, downsample = 200)
aer_vivo_mat <- aer_vivo@assays$RNA@data[AER, cell_idx2]
aer_vivo_mat <- as.matrix(aer_vivo_mat)

aer_mat <- cbind(aer_vitro_mat, aer_vivo_mat)


cellid <- c(
    as.character(aer_vitro$orig.ident[cell_idx1]),
    as.character(aer_vivo$orig.ident[cell_idx2])
)
cellid <- factor(cellid, levels = c("day7", "Mouse_E95", "Mouse_E105", "Mouse_E115", "Mouse_E125"))

top_anno_a <- HeatmapAnnotation(
    Condition = cellid,
    col = list(Condition = col.cond)
)

# order cells by Fgf8
df_split <- split.data.frame(as.data.frame(t(aer_mat)), f = cellid)
cn_list <- lapply(df_split, function(df) {
    rownames(df)[order(df[, "Fgf8"], decreasing = F)]
})

heat_a <- Heatmap(aer_mat,
    col = colorRamp2(c(0, 0.5, 2, 3), cividis(4)), # cividis(100),
    cluster_columns = F,
    cluster_rows = F,
    row_order = AER,
    show_row_dend = F,
    show_column_names = F,
    column_order = unlist(cn_list),
    top_annotation = top_anno_a,
    column_split = cellid,
    column_title = "AER marker expression (Downsampling 200 cells)",
    use_raster = F
)

pdf("AER_markers.pdf", width = 10, height = 5)
heat_a
dev.off()


## ================= MetaNeighbor ==========
library(MetaNeighbor)
library(SummarizedExperiment)
source("./pl_helper.r")

# AER, Mes
seu_sub <- merge(aer_vitro, list(aer_vivo, mes_vitro, mes_vivo))
seu_sub$celltype[grepl("^AER", seu_sub$celltype)] <- "AER"
seu_sub$celltype[grepl("^Mes|CT", seu_sub$celltype)] <- "Mesoderm"
unique(seu_sub$celltype)

se <- SummarizedExperiment(seu_sub@assays$RNA@counts, colData = seu_sub@meta.data)
var_genes <- variableGenes(dat = se, exp_labels = se$orig.ident)
length(var_genes)

AUROC.scores <- MetaNeighborUS(
    var_genes = var_genes, # 4082
    dat = se,
    study_id = se$orig.ident,
    cell_type = se$celltype,
    fast_version = T
)
# topHitsByStudy(AUROC.scores)

write.csv(AUROC.scores, "MetaNeigbor_onlyd7_mes_ect.csv")

pdf("MetaNeigbor_onlyd7_mes_ect_0.95.pdf")
plot_AUROC_heatmap(AUROC.scores)
dev.off()



# AER, BE, Mes, Muscle
seu_sub <- seu_sub <- merge(aer_vitro, list(aer_vivo, mes_vitro, mes_vivo, be_vivo, muscle_vivo))
seu_sub$celltype[grepl("^AER", seu_sub$celltype)] <- "AER"
seu_sub$celltype[grepl("^Mes|CT", seu_sub$celltype)] <- "Mesoderm"
unique(seu_sub$celltype)

se <- SummarizedExperiment(seu_sub@assays$RNA@counts, colData = seu_sub@meta.data)
var_genes <- variableGenes(dat = se, exp_labels = se$orig.ident)


AUROC.scores <- MetaNeighborUS(
    var_genes = var_genes, # 4124
    dat = se,
    study_id = se$orig.ident,
    cell_type = se$celltype,
    fast_version = T
)
topHitsByStudy(AUROC.scores)

write.csv(AUROC.scores, "MetaNeigbor_onlyd7_mes_ect_muscle.csv")

pdf("MetaNeigbor_onlyd7_mes_ect_muscle_0.95.pdf")
plot_AUROC_heatmap(AUROC.scores)
dev.off()
