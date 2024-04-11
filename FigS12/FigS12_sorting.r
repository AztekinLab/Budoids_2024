
# supp12
# cell sorting
# day7 only
library(Seurat)
library(patchwork)

setwd("/work/gr-aztekin/3.project/culture_AER_final/FigS12")
combined <- readRDS("../Fig1/mESC_2D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"

combined2 <- subset(combined, subset = orig.ident == "day7")

perc <- 100*table(combined2$celltype)/dim(combined2)[2]
perc[perc >5]
# Surface Ectoderm 1 Surface Ectoderm 2           AER-like         Mesodermal
#           34.92420           28.01797           20.13850           13.98091

pdf("D7only_epcam_cd9_cd44.pdf", width = 12, height = 3)

cd9 <- FeaturePlot(combined, c("Epcam", "Cd9"),
    order = T,
    blend = T, blend.threshold = 0, pt.size = 0.1,
    cells = WhichCells(combined, expression = orig.ident == "day7")
)

# cd44 <- FeaturePlot(combined, c("Epcam", "Cd44"),
#     order = T,
#     blend = T, blend.threshold = 0, pt.size = 0.1,
#     cells = WhichCells(combined, expression = orig.ident == "day7")
# )

# cd9 / cd44
cd9
dev.off()




col.ct <- c(
    "#dcc3fd", "#9f96f1",
    "#fbcea3", "#f99e4a", "#db995a", "#ad733c",
    "#8E0068", "#94c5e5",
    "#94b44f", "#2176ff",
    "grey"
)
names(col.ct) <- levels(combined$celltype)


pdf("D7only_epcam_cd9_cd44_vln.pdf", width = 8, height = 5)
VlnPlot(subset(combined, orig.ident == "day7"), c("Epcam", "Cd9"),
    idents = c("Surface Ectoderm 1", "Surface Ectoderm 2", "AER-like", "Mesodermal"), same.y.lims = T, pt.size = 0, cols = col.ct
)
dev.off()



pdf("D7only_epcam_cd9_cd44_ridge.pdf", width = 10, height = 3)
RidgePlot(subset(combined, orig.ident == "day7"), c("Epcam", "Cd9"),
    idents = c("Surface Ectoderm 1", "Surface Ectoderm 2", "AER-like", "Mesodermal"), same.y.lims = T, cols = col.ct
)
dev.off()



pdf("D7only_epcam_cd9_cd44_ridge_all.pdf", width = 10, height = 5)
RidgePlot(subset(combined, orig.ident == "day7"), c("Epcam", "Cd9"), same.y.lims = T, cols = col.ct)
dev.off()
