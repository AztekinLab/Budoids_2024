library(Seurat)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

library(ComplexHeatmap)
library(circlize)
library(viridis)

setwd("/work/gr-aztekin/3.project/culture_AER_final/FigSx")
combined <- readRDS("../Fig3/mESC_3D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"


# prepare in vitro dataset
aer_vitro <- subset(combined, cells = WhichCells(combined, idents = "AER-like"))
be_vitro <- subset(combined, cells = WhichCells(combined, idents = "BasalEctoderm-like"))
mes_vitro <- subset(combined, cells = WhichCells(combined, idents = c("Mesodermal-like", "Mes_Hoxb", "Mes_Hox13", "Mes_Sox9", "Mes_Col2a1")))

aer_vitro$celltype <- as.character(aer_vitro$celltype)
be_vitro$celltype <- as.character(be_vitro$celltype)
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



## ================= MetaNeighbor ==========
library(MetaNeighbor)
library(SummarizedExperiment)
source("../FigS4_and_S5/pl_helper.r")

# AER, Mes
seu_sub <- merge(aer_vitro, list(be_vitro, mes_vitro, aer_vivo, be_vivo, mes_vivo))
# seu_sub$celltype[grepl("^AER", seu_sub$celltype)] <- "AER"
seu_sub$celltype[grepl("CT", seu_sub$celltype)] <- "Mesoderm"
unique(seu_sub$celltype)

se <- SummarizedExperiment(seu_sub@assays$RNA@counts, colData = seu_sub@meta.data)
var_genes <- variableGenes(dat = se, exp_labels = se$orig.ident)
length(var_genes)

AUROC.scores <- MetaNeighborUS(
    var_genes = var_genes, # 2705
    dat = se,
    study_id = se$orig.ident,
    cell_type = se$celltype,
    fast_version = T
)
# topHitsByStudy(AUROC.scores)

write.csv(AUROC.scores, "MetaNeigbor_mes_ect.csv")

rs <- stringi::stri_replace_all_regex(rownames(AUROC.scores),
    pattern = c(".*\\|", "-.*"), replacement = c("", ""), vectorize = FALSE
)
rs[grep("^Mes", rs)] <- "Mesoderm"


pdf("MetaNeigbor_mes_ect_0.95.pdf", width = 10, height = 9)
plot_AUROC_heatmap(AUROC.scores, Thres = 0.95, column_split = rs, row_split = rs)
dev.off()
