library(Seurat)
library(dplyr)

# get all mat
file_list1 <- Sys.glob("./batch1/*/count/sample_filtered_feature_bc_matrix")
file_list2 <- Sys.glob("./batch2/*/count/sample_filtered_feature_bc_matrix")
file_list3 <- Sys.glob("./Day*/outs/filtered_feature_bc_matrix")

names(file_list1) <-c("day0", "day3", "day5", "day7")
names(file_list2) <- c("day5", "day7", "day12")
names(file_list3) <- c("day9", "day12")


# batch1 data preprocessing
seu_list1 <- lapply(1:length(file_list1), function(x) {
    dat <- Read10X(file_list1[x])
    seu <- CreateSeuratObject(dat$`Gene Expression`, min.cells = 3)

    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
    seu$orig.ident <- names(file_list1)[x]
    seu$batch <- "Rep1"
    seu
})


seu_list1[[1]] <- subset(seu_list1[[1]], subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & nCount_RNA > 12500 & nCount_RNA < 75000 & percent.mt < 10)
seu_list1[[2]] <- subset(seu_list1[[2]], subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & nCount_RNA > 12500 & nCount_RNA < 75000 & percent.mt < 10)
seu_list1[[3]] <- subset(seu_list1[[3]], subset = nFeature_RNA > 2000 & nFeature_RNA < 9500 & nCount_RNA > 10000 & nCount_RNA < 140000 & percent.mt < 10)
seu_list1[[4]] <- subset(seu_list1[[4]], subset = nFeature_RNA > 1250 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 10)




# batch2 data preprocessing
seu_list2 <- lapply(1:length(file_list2), function(x) {
    dat <- Read10X(file_list2[x])
    seu <- CreateSeuratObject(dat$`Gene Expression`)

    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
    seu$orig.ident <- names(file_list2)[x]
    seu$batch <- "Rep2"
    seu
})

seu_list2[[1]] <- subset(seu_list2[[1]], subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & nCount_RNA > 200 & nCount_RNA < 100000 & percent.mt < 10)
seu_list2[[2]] <- subset(seu_list2[[2]], subset = nFeature_RNA > 1500 & nFeature_RNA < 9500 & nCount_RNA > 200 & nCount_RNA < 100000 & percent.mt < 10)
seu_list2[[3]] <- subset(seu_list2[[3]], subset = nFeature_RNA > 1500 & nFeature_RNA < 7500 & nCount_RNA > 200 & nCount_RNA < 100000 & percent.mt < 10)



# batch3 data preprocessing
seu_list3 <- lapply(1:length(file_list3), function(x) {
    dat <- Read10X(file_list3[x])
    seu <- CreateSeuratObject(dat)

    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
    seu$orig.ident <- names(file_list3)[x]
    seu$batch <- "Rep3"
    seu
})


seu_list3[[1]] <- subset(seu_list3[[1]], subset = nFeature_RNA > 2500 & nFeature_RNA < 11000 & nCount_RNA > 2000 & nCount_RNA < 150000 & percent.mt < 10)
seu_list3[[2]] <- subset(seu_list3[[2]], subset = nFeature_RNA > 2500 & nFeature_RNA < 11000 & nCount_RNA > 2000 & nCount_RNA < 150000 & percent.mt < 10)



# Cell Cycle Scoring
s.genes <- Hmisc::capitalize(tolower(cc.genes.updated.2019$s.genes))
g2m.genes <- Hmisc::capitalize(tolower(cc.genes.updated.2019$g2m.genes))


seu1 <- merge(seu_list1[[1]], seu_list1[-1])
seu2 <- merge(seu_list2[[1]], seu_list2[-1])
seu_list <- lapply(X = c(list(seu1, seu2), seu_list3), FUN = function(x) {
    x %>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})


# save into rds
dir.create('processed')
setwd('processed')
mapply(saveRDS, seu_list, c("batch1.rds", "batch2.rds", "batch3_d9.rds", "batch3_d12.rds"))
