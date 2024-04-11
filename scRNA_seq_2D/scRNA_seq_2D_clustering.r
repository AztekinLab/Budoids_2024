library(Seurat)
library(dplyr)

setwd("/work/gr-aztekin/3.project/culture_AER_final/scRNA_seq_2D")
seu1 <- readRDS("../1.data/batch1.rds")
seu2 <- readRDS("../1.data/batch2.rds")
seu2 <- subset(seu2, orig.ident != "day12")

seu_list <- list(seu1, seu2)

# rPCA based integration
features <- SelectIntegrationFeatures(object.list = seu_list)
seu_list <- lapply(X = seu_list, FUN = function(seu) {
    seu %>%
        ScaleData(features = features, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score")) %>%
        RunPCA(features = features)
})


anchors <- FindIntegrationAnchors(seu_list, anchor.features = features, reduction = "rpca")
combined <- IntegrateData(anchorset = anchors)

# DefaultAssay(combined) <- "integrated"

combined <- combined %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 0.5)


combined <- FindSubCluster(combined, cluster = 6, resolution = 0.3, graph.name = "integrated_snn")


# Cell Type annotation
combined$orig.ident <- factor(combined$orig.ident, levels = sort(unique(combined$orig.ident)))
combined$celltype <- as.character(combined$integrated_snn_res.0.5)
combined$celltype[combined$celltype %in% c(0, 8)] <- "Surface Ectoderm 1" # ffc58e
combined$celltype[combined$celltype == 1] <- "Surface Ectoderm 2" # b38a63
combined$celltype[combined$celltype == 2] <- "AER-like" # 8E0068
combined$celltype[combined$celltype == 3] <- "Surface Ectoderm 3" # ffd1a5
combined$celltype[combined$celltype == 6] <- "Primed mESC" # cc9e72
combined$celltype[combined$celltype == 5] <- "Surface Ectoderm 4" # 00BF93
combined$celltype[combined$celltype == 4] <- "Mesodermal" # 7669EB
combined$celltype[combined$celltype == 7] <- "Naïve mESC" # 9f96f1
combined$celltype[combined$sub.cluster == "6_2"] <- "Neuroectoderm" # f88e36
combined$celltype[combined$celltype == 9] <- "Unknown" #
combined$celltype[combined$celltype == 10] <- "Sox17+ cells" # 4d1a25

combined$celltype <- factor(combined$celltype,
    levels = c(
        "Naïve mESC", "Primed mESC",
        "Surface Ectoderm 1", "Surface Ectoderm 2", "Surface Ectoderm 3", "Surface Ectoderm 4",
        "AER-like", "Mesodermal", "Neuroectoderm", "Sox17+ cells",
        "Unknown"
    )
)

saveRDS(combined, file = "mESC_2D.rds")
