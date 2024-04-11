library(Seurat)
library(ggpubr)
library(patchwork)

setwd("/work/gr-aztekin/3.project/culture_AER_final/Fig3")
combined <- readRDS("mESC_3D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"

# scGESA
mes_cluster <- c("Mesodermal-like", "Mes_Sox9", "Mes_Col2a1", "Mes_Hoxb", "Mes_Hox13")
aer_cluster <- "AER-like"

aer <- subset(combined, idents = aer_cluster)
mes <- subset(combined, idents = mes_cluster)

## marker expression
markers <- readxl::read_xlsx("../marker.xlsx")
markers <- as.data.frame(markers)


genelist <- readRDS("/work/gr-aztekin/3.project/axolotl_AER_final/functional_list_updateAxolotl_updateEctoderm_updateNames.rds")
aer_list <- na.omit(unique(c(genelist$Ectoderm$marker$Mouse$gene_name[-c(15, 24, 25, 32, 39)], markers$AER))) # "Fgf15" "Lgr5"  "Lgr6"  "Rspo3" "Wnt3a", not mice AER marker
chon_list <- na.omit(unique(c(genelist$differentiation$chondro$Mouse$gene_name[-16], markers$Chondrocyte))) # sox10
fibro_list <- na.omit(unique(c(genelist$differentiation$fibro$Mouse$gene_name, markers$Fibro)))
teno_list <- na.omit(unique(c(genelist$differentiation$teno$Mouse$gene_name, markers$Tenocyte)))
multipotent_list <- na.omit(unique(c(genelist$CT$stemness$Mouse$gene_name, markers$Distal_limb_mesoderm)))
# matureBE_list <- c('Krt10','Krt1')

multipotent_list <- c("Igfbp3", "Adamts15", "Alb", "C4bp", "Col23a1", "Col8a1", "Cxcl12", "Dact1", "Dlk1", "Ebf2", "Emx2", "Gsc", "H19", "Hoxc4", "Hs3st3b1", "Igf1", "Igf2", "Igfbp5", "Itm2a", "Meg3", "Meis2", "Meox1", "Mir483", "Mir675", "Nr2f1", "Ntn1", "Ntrk2", "Nxnl2", "Osr2", "Pdzrn3", "Pitx2", "Pkdcc", "Rgcc", "Ror1", "Scx", "Zcchc5", "Eno3", "Gria2", "Hmga2", "Hoxb3", "Hoxb4", "Hoxb5", "Hoxc6", "Irx3", "Lbx1", "Meis1", "Peg3", "Stfa2", "Tcf15", "Tshz3", "Zim1")

input_list <- list(
    chon = chon_list, fibro = fibro_list, teno = teno_list,
    multi = multipotent_list
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

p4 <- VlnPlot(mes, features = "multi4", pt.size = 0, group.by = "orig.ident", cols = col_time, y.max = max(mes$multi4) + 0.1) &
    stat_compare_means(
        comparisons = list(c("day7", "day12")),
        label = "p.signif",
        method = "wilcox.test"
    )

pdf("cell_fate2.pdf", width = 12, height = 3)
wrap_plots(p1, p2, p3, p4, nrow = 1, guides = "collect")
dev.off()
