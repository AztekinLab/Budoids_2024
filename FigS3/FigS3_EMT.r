library(Seurat)
library(ComplexHeatmap)
library(viridis)
library(circlize)

setwd("/work/gr-aztekin/3.project/culture_AER_final/FigS3")
combined <- readRDS("../Fig1/mESC_2D.rds")
DefaultAssay(combined) <- "RNA"


# EMT D5
EMT <- c("Cdh1", "Cdh2", "Epcam", "Zeb1", "Zeb2", "Vim", "Snai1", "Snai2", "Twist1")
seu_d5 <- subset(combined, orig.ident == "day5" & celltype %in% c("Surface Ectoderm 3", "Surface Ectoderm 4"))


pheat <- Heatmap(as.matrix(seu_d5@assays$RNA@data[EMT, ]),
    col = colorRamp2(c(0, 1, 2, 3), cividis(4)), # cividis(100),
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

pdf("EMT_heatmap.pdf", width = 8, height = 4)
pheat
dev.off()


pdf("EMT_umap.pdf", width = 12, height = 12)
FeaturePlot(seu_d5, EMT,
    ncol = 3,
    cols = c("lightgrey", "red4"),
    order = T, raster = T, pt.size = 2
)
dev.off()
