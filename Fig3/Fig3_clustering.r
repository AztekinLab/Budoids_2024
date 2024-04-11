library(Seurat)
library(dplyr)
source("../../axolotl_AER_final/0.scripts/cell_cycle_removal.r")

setwd("/work/gr-aztekin/3.project/culture_AER_final/Fig3")
seu1 <- readRDS("../1.data/processed/batch1.rds")
seu2 <- readRDS("../1.data/processed/batch2.rds")
seu31 <- readRDS("../1.data/processed/batch3_d9.rds")
seu32 <- readRDS("../1.data/processed/batch3_d12.rds")

seu1 <- subset(seu1, orig.ident == "day7")
seu2 <- subset(seu2, orig.ident %in% c("day7", "day12"))

seu_list <- list(seu1, seu2, seu31, seu32)



# rPCA based integration
seu_list <- lapply(seu_list, FUN = SCTransform, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"))

features <- SelectIntegrationFeatures(seu_list, nfeatures = 3000)
seu_list <- PrepSCTIntegration(seu_list, anchor.features = features)
seu_list <- lapply(X = seu_list, FUN = RunPCA, features = features)


anchors <- FindIntegrationAnchors(seu_list, anchor.features = features, reduction = "rpca", normalization.method = "SCT")

sample.tree <- matrix(c(-3, -4, 1, -2, 2, -1), byrow = T, ncol = 2)
combined <- IntegrateData(anchorset = anchors, sample.tree = sample.tree, normalization.method = "SCT")


combined <- combined %>%
    RunPCA() %>%
    RunUMAP(dims = 1:40, min.dist = .1, spread = 1.25, n.neighbors = 100) %>%
    FindNeighbors(dims = 1:40, k.param = 100) %>%
    FindClusters(resolution = 0.8)

# fine tune clustering
DefaultAssay(combined) <- "integrated"
combined <- FindSubCluster(combined, cluster = 0, resolution = 0.3, graph.name = "integrated_snn", subcluster.name = "C0_sub")
combined <- FindSubCluster(combined, cluster = 1, resolution = 0.3, graph.name = "integrated_snn", subcluster.name = "C1_sub")
combined <- FindSubCluster(combined, cluster = 6, resolution = 0.3, graph.name = "integrated_snn", subcluster.name = "C6_sub")


# meta <- readRDS("../Fig1/metadata.rds")
# meta <- subset(meta, orig.ident == "day7")
# rownames(meta) <- paste(rownames(meta), gsub("Rep", "", meta$batch), sep = "_")
# all(rownames(meta) %in% rownames(combined@meta.data))

# cells1 <- rownames(meta)[meta$celltype == "Surface Ectoderm 1"]
# cells2 <- rownames(meta)[meta$celltype == "Surface Ectoderm 2"]
# cells3 <- rownames(meta)[meta$celltype == "AER-like"]
# cells4 <- rownames(meta)[meta$celltype == "Mesodermal"]

# pdf("clustering150.pdf", width = 8, height = 8)
# DimPlot(combined, label = TRUE, group.by = "batch")
# DimPlot(combined, label = TRUE, group.by = "Phase")
# DimPlot(combined, label = TRUE, group.by = "orig.ident")
# DimPlot(combined, label = TRUE)

# DimPlot(combined, label = TRUE, cells.highlight = cells1)
# DimPlot(combined, label = TRUE, cells.highlight = cells2)
# DimPlot(combined, label = TRUE, cells.highlight = cells3)
# DimPlot(combined, label = TRUE, cells.highlight = cells4)

# DimPlot(combined, label = TRUE, group.by = "C0_sub")
# DimPlot(combined, label = TRUE, group.by = "C1_sub")
# DimPlot(combined, label = TRUE, group.by = "C6_sub")

# DimPlot(combined, label = TRUE, group.by = "celltype")

# FeaturePlot(combined, c("Fgf8", "Prrx1"), order = T)
# FeaturePlot(combined, c("Hoxa13", "Hoxd13"), order = T)
# FeaturePlot(combined, c("nFeature_RNA"), order = T)
# dev.off()



# DefaultAssay(combined) <- "RNA"
# pdf("marker.pdf")
# FeaturePlot(combined, c("Fgf8", "Epcam", "Trp63", "Krtdap"), order = T, label = T)
# FeaturePlot(combined, c("Krt14", "Krt5"), order = T, label = T)
# FeaturePlot(combined, c("Krt1", "Krt10"), order = T, label = T)
# FeaturePlot(combined, c("Krt8", "Krt18"), order = T, label = T)

# FeaturePlot(combined, c("Prrx1", "Sox9", "Twist2", "Hoxd13"), order = T)
# FeaturePlot(combined, c("Hoxb9", "Hoxd13"), order = T)
# FeaturePlot(combined, c("Sox9", "Col2a1", "Col9a1"), order = T)

# FeaturePlot(combined, c("Nanog", "Sox2"), order = T)
# FeaturePlot(combined, "Sox17", order = T) # endoderm
# FeaturePlot(combined, c("Ttr", "Otx2", "Bmp7"), order = T) # choroid plexus
# FeaturePlot(combined, c("Sox2", "Pax6"), order = T) # neuralectoderm
# FeaturePlot(combined, c("Sox10", "Pax3"), order = T) # neural crest
# FeaturePlot(combined, c("Neurod1", "Tubb3"), order = T) # neuron
# dev.off()

# mar <- FindAllMarkers(combined, only.pos = T)
# mar <- subset(mar, p_val_adj < 0.05)
# write.csv(mar, file = "findallmarker.csv")


# cell type annotation
combined$celltype <- as.character(combined$seurat_clusters)
combined$celltype[combined$celltype %in% c(2, 8)] <- "Surface Ectoderm 1"
combined$celltype[combined$C0_sub == "0_0"] <- "Surface Ectoderm 2"
combined$celltype[combined$C0_sub == "0_1"] <- "AER-like"
combined$celltype[combined$celltype == 6] <- "AER-like"
combined$celltype[combined$C6_sub == "6_0"] <- "BasalEctoderm-like"


combined$celltype[combined$celltype %in% c(1, 4)] <- "Mes_Hoxb"
combined$celltype[combined$C1_sub == "1_1"] <- "Mes_Hox13"
combined$celltype[combined$celltype == 3] <- "Mes_Sox9"
combined$celltype[combined$celltype == 5] <- "Mes_Col2a1"
combined$celltype[combined$celltype %in% c(12, 13)] <- "Mesodermal-like"

combined$celltype[combined$celltype == 7] <- "Pluripotent"
combined$celltype[combined$celltype == 9] <- "Endoderm"
combined$celltype[combined$celltype == 10] <- "NeuralEctoderm"
combined$celltype[combined$celltype == 11] <- "ChoroidPlexus"
combined$celltype[combined$celltype == 14] <- "Unknown"
combined$celltype[combined$celltype == 15] <- "Neuron"
combined$celltype[combined$celltype == 16] <- "NeuralCrest"


le_ct <- c(
    "Surface Ectoderm 1", "Surface Ectoderm 2", "BasalEctoderm-like", "AER-like",
    "Mesodermal-like", "Mes_Hoxb", "Mes_Hox13", "Mes_Sox9", "Mes_Col2a1",
    "NeuralEctoderm", "Neuron", "ChoroidPlexus", "NeuralCrest",
    "Endoderm", "Pluripotent", "Unknown"
)
combined$celltype <- factor(combined$celltype, levels = le_ct)
combined$orig.ident <- factor(combined$orig.ident, levels = c("day7","day9","day12" ))


saveRDS(combined, "mESC_3D.rds")
