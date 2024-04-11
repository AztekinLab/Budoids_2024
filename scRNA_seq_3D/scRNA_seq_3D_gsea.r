library(Seurat)
library(ggpubr)
library(patchwork)

setwd("/work/gr-aztekin/3.project/culture_AER_final/scRNA_seq_3D")
combined <- readRDS("mESC_3D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"

# scGESA
mes_cluster <- c("Mesodermal-like", "Mes_Sox9", "Mes_Col2a1", "Mes_Hoxb", "Mes_Hox13")
aer_cluster <- "AER-like"

aer <- subset(combined, idents = aer_cluster)
mes <- subset(combined, idents = mes_cluster)

## get genes
genelist <- read.csv( 'genelist_gsea.csv')
aer_list <- na.omit(genelist$AER)
chon_list <- na.omit(genelist$Chondrocyte)
fibro_list <- na.omit(genelist$Fibroblast)


input_list <- list(
    chon = chon_list, fibro = fibro_list
)

aer <- AddModuleScore(aer,
    features = list(aer = aer_list),
    assay = "RNA", name = "aer"
)

mes <- AddModuleScore(mes,
    features = input_list,
    assay = "RNA", name = names(input_list)
)


col_time <- c("#7bb3cd", "#548ba5", "#46606b")

p1 <- VlnPlot(aer, features = "aer1", pt.size = 0, group.by = "orig.ident", cols = col_time, y.max = max(aer$aer1) + 0.1) &
    stat_compare_means(
        comparisons = list(c("day7", "day12")),
        label = "p.signif",
        method = "wilcox.test"
    )

p2 <- VlnPlot(mes, features = "chon1", pt.size = 0, group.by = "orig.ident", cols = col_time, y.max = max(mes$chon1) + 0.1) &
    stat_compare_means(
        comparisons = list(c("day7", "day12")),
        label = "p.signif",
        method = "wilcox.test"
    )
p3 <- VlnPlot(mes, features = "fibro2", pt.size = 0, group.by = "orig.ident", cols = col_time, y.max = max(mes$fibro2) + 0.1) &
    stat_compare_means(
        comparisons = list(c("day7", "day12")),
        label = "p.signif",
        method = "wilcox.test"
    )

pdf("Fig3F_cell_fate.pdf", width = 10, height = 3)
wrap_plots(p1, p2, p3, nrow = 1, guides = "collect")
dev.off()
