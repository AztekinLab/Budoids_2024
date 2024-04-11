library(Seurat)
library(dplyr)
source("../../axolotl_AER_final/0.scripts/cell_cycle_removal.r")

setwd("/work/gr-aztekin/3.project/culture_AER_final/scRNA_seq_3D")
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
