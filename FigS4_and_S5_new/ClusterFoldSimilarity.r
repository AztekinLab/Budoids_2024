library(Seurat)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

library(ComplexHeatmap)
library(circlize)
library(viridis)

setwd("/work/gr-aztekin/3.project/culture_AER_final/FigS4_and_S5_new")
combined <- readRDS("../Fig1/mESC_2D.rds")
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"

# prepare in vitro dataset
aer_vitro <- subset(combined, cells = WhichCells(combined, idents = c("AER-like"), expression = orig.ident == "day7"))
mes_vitro <- subset(combined, cells = WhichCells(combined, idents = "Mesodermal", expression = orig.ident == "day7"))
aer_vitro$celltype <- as.character(aer_vitro$celltype)
mes_vitro$celltype <- as.character(mes_vitro$celltype)





# prepare in vitro dataset
combined <- readRDS("../Fig3/mESC_3D.rds")
DefaultAssay(combined) <- 'RNA'
combined[['integrated']] <- NULL
combined[['SCT']] <- NULL
Idents(combined) <- "celltype"

aer_vitro2 <- subset(combined, cells = WhichCells(combined, idents = "AER-like"))
be_vitro2 <- subset(combined, cells = WhichCells(combined, idents = "BasalEctoderm-like"))
mes_vitro2 <- subset(combined, cells = WhichCells(combined, idents = c("Mesodermal-like", "Mes_Hoxb", "Mes_Hox13", "Mes_Sox9", "Mes_Col2a1")))

aer_vitro2$celltype <- as.character(aer_vitro2$celltype)
be_vitro2$celltype <- as.character(be_vitro2$celltype)
mes_vitro2$celltype <- as.character(mes_vitro2$celltype)



# prepare in vivo dataset
vivo <- readRDS('/work/gr-aztekin/3.project/limb_dev/1.data/human_E-MTAB-8813/processed/MouseAllData_annotated_filtered_subset.rds')
vivo$orig.ident <- vivo$stage

pdf('umap.pdf')
DimPlot(vivo2[, WhichCells(vivo, expression =  celltype %in% c('AER-Basal','Basal'))], group.by = 'celltype', split.by = 'stage', ncol = 3, pt.size =0.1)

dev.off()


freq <- table(vivo$celltype, vivo$orig.ident)
cells <- c()
for(i in colnames(freq)){
    subdf <- freq[,colnames(freq) == i]
    need <- names(subdf)[subdf > 15]
    cells <- c(cells, WhichCells(vivo, expression = orig.ident == i & celltype %in% need))
}
vivo2 <- subset(vivo, cells = cells)
table(vivo2$celltype, vivo2$orig.ident)


aer_vivo <- subset(vivo2, subset = celltype == "AER-Basal")
Idents(aer_vivo) <- "stage"

be_vivo <- subset(vivo2, subset = celltype == "Basal")
Idents(be_vivo) <- "stage"


meso_cluster <- c( 'Adh+Fibro', 'DistalMes', 'EarlyDistalMes',
       'EarlyProxMes', 'Hoxc5+DermFibroProg', 'InterZone', 'Meox2+Mes',
       'MesCond', 'PrehyperChon', 'ProlifChon', 'ProxMes',
       'Rdh10+DistalMes', 'RestingChon', 'TransMes')
mes_vivo <- subset(vivo2, subset = celltype %in% meso_cluster)
Idents(mes_vivo) <- "stage"



# remotes::install_github("OscarGVelasco/ClusterFoldSimilarity")
library(ClusterFoldSimilarity)

sc1 <- merge(aer_vivo, list(be_vivo, mes_vivo))
sc2 <- merge(aer_vitro, mes_vitro)
sc3 <- merge(aer_vitro2, list(be_vitro2, mes_vitro2))

freq <- table(sc3$celltype, sc3$orig.ident)
cells <- c()
for(i in colnames(freq)){
    subdf <- freq[,colnames(freq) == i]
    need <- names(subdf)[subdf > 35]
    cells <- c(cells, WhichCells(sc3, expression = orig.ident == i & celltype %in% need))
}
sc32 <- subset(sc3, cells = cells)
table(sc32$celltype, sc32$orig.ident)




scList <- list(sc2, sc1)
scList2 <- list(sc32, sc1)

scList <- lapply(scList, function(x){
    x$comp <- paste(x$orig.ident, x$celltype, sep = '|')
    Idents(x) <- factor(x$comp)
    x <- FindVariableFeatures(x, nfeatures = 3000)
    x
})


scList2 <- lapply(scList2, function(x){
    x$comp <- paste(x$orig.ident, x$celltype, sep = '|')
    Idents(x) <- factor(x$comp)
    x <- FindVariableFeatures(x, nfeatures = 3000)
    x
})

sc_vf <- lapply(scList, function(x){x[VariableFeatures(x),]})
sc_vf2 <- lapply(scList2, function(x){x[VariableFeatures(x),]})

# lapply(scList, function(x){print(length(VariableFeatures(x)))})

# BiocParallel::register(BPPARAM =  BiocParallel::MulticoreParam(workers = 6))
pdf('topinf.pdf', width = 25)
similarityTable <- clusterFoldSimilarity(scList=sc_vf, nSubsampling=25, sampleNames=c("vitro", "vivo"),topN=Inf, parallel=F)
dev.off()



# BiocParallel::register(BPPARAM =  BiocParallel::MulticoreParam(workers = 6))
pdf('topinf_3d.pdf', width = 25)
similarityTable <- clusterFoldSimilarity(scList=sc_vf2, nSubsampling=25, sampleNames=c("vitro", "vivo"),topN=Inf, parallel=F)
dev.off()



df <- subset(similarityTable, datasetL == 'vitro')
df$comp1 <- paste(df$datasetL,df$clusterL, sep = "|")
df$comp2 <- paste(df$datasetR,df$clusterR, sep = "|")
df <- df[,c(1, 9,10)]

df_wide <- reshape(df, idvar = "comp1", timevar = "comp2",  v.names="similarityValue", direction = "wide")

rownames(df_wide) <- df_wide[,1]
df_wide <- df_wide[,-1]


r_keep <- rowSums(is.na(df_wide)) != nrow(df_wide)
c_keep <- colSums(is.na(df_wide)) != ncol(df_wide)
df_wide <- df_wide[r_keep, c_keep]


## row/col labels
r_lab <- stringi::stri_replace_all_regex(rownames(df_wide),
    pattern = c(".*\\|"), replacement = c(""), vectorize = FALSE
)
c_lab <- stringi::stri_replace_all_regex(colnames(df_wide),
    pattern = c(".*\\|"), replacement = c(""), vectorize = FALSE
)


r_lab1 <- unlist(lapply(strsplit(rownames(df_wide), "\\|"), `[`, 2))
c_lab1 <- unlist(lapply(strsplit(colnames(df_wide), "\\|"), `[`, 2))


# row/col annotation
num_col <- length(unique(c_lab1)) + 1
cols <- brewer.pal(max(3, num_col), "Greys")[2:num_col]
names(cols) <- sort(as.numeric((unique(c_lab1))))
# "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#005A32"

c_anno <- HeatmapAnnotation(
    cond = c_lab1,
    col = list(cond = cols),
    show_annotation_name = FALSE
)


num_col <- length(unique(r_lab1)) + 1
cols <- brewer.pal(max(3, num_col), "Oranges")[2:num_col]
names(cols) <- stringr::str_sort(unique(r_lab1), numeric = TRUE)
r_anno <- rowAnnotation(
    cond = r_lab1,
    col = list(cond = cols),
    show_annotation_name = FALSE
)


maxs <- apply(df_wide, 1, max)

set.seed(123) # reproducible column_km
p <- Heatmap(df_wide, col = c("steelblue", "white", "red3"),
    width = ncol(df_wide)*unit(5, "mm"),
    height = nrow(df_wide)*unit(5, "mm"),
column_dend_height = unit(3, "cm"),
show_row_dend = F, column_km = 2, row_km = 2,
        # annotation
        row_labels = r_lab, column_labels = c_lab,
        bottom_annotation = c_anno,
        right_annotation = r_anno,

        rect_gp = gpar(col = "white", lwd = 1),
         column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),

        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid

        if((df_wide[i, j] == max(df_wide[i,]))){
grid.text('+', x, y)
        }

        }

        )




pdf('test_3d.pdf', width = 20, height = 7)
# draw(p, annotation_legend_list = list(lgd), merge_legend = TRUE)
draw(p)
dev.off()


pdf('test.pdf', width = 15, height = 4)
# draw(p, annotation_legend_list = list(lgd), merge_legend = TRUE)
draw(p)
dev.off()
