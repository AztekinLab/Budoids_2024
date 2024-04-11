
library(RColorBrewer)
library(ComplexHeatmap)

plot_AUROC_heatmap <- function(AUROC.scores, Thres = 0.95, row_split = NULL, column_split = NULL) {
    r_keep <- rowSums(is.na(AUROC.scores)) != nrow(AUROC.scores)
    c_keep <- colSums(is.na(AUROC.scores)) != ncol(AUROC.scores)
    AUROC.scores <- AUROC.scores[r_keep, c_keep]


    ## row/col labels
    r_lab <- stringi::stri_replace_all_regex(rownames(AUROC.scores),
        pattern = c(".*\\|"), replacement = c(""), vectorize = FALSE
    )
    c_lab <- r_lab

    r_lab1 <- stringi::stri_replace_all_regex(rownames(AUROC.scores),
        pattern = c("\\|.*"), replacement = c(""), vectorize = FALSE
    )
    # r_lab1 <- r_lab
    c_lab1 <- r_lab1


    # row/col annotation
    num_col <- length(unique(r_lab1))
    cols <- brewer.pal(max(3, num_col), "Accent")[1:num_col]
    names(cols) <- unique(r_lab1)

    c_anno <- HeatmapAnnotation(
        cond = anno_points(rep(1, times = nrow(AUROC.scores)),
            ylim = c(0, 1),
            size = unit(5, "mm"),
            height = unit(1, "mm"),
            border = FALSE,
            axis = FALSE,
            gp = gpar(col = cols[r_lab1])
        ),
        show_annotation_name = FALSE
    )


    r_anno <- rowAnnotation(
        cond = anno_points(rep(1, times = nrow(AUROC.scores)),
            ylim = c(0, 1),
            size = unit(5, "mm"),
            width = unit(1, "mm"),
            border = FALSE,
            axis = FALSE,
            gp = gpar(col = cols[r_lab1])
        ),
        show_annotation_name = FALSE
    )

    lgd <- Legend(
        labels = names(cols), type = "points", pch = 19,
        background = NULL, size = unit(3, "mm"),
        legend_gp = gpar(fill = cols, col = cols),
    )
    lgd_sig <- Legend(
        labels = paste(">", Thres), type = "points", pch = "*",
        background = NULL, size = unit(5, "mm")
    )


    p <- Heatmap(AUROC.scores,
        col = c("steelblue", "white", "red3"), border = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (AUROC.scores[i, j] > Thres) {
                grid.text("*", x, y, gp = gpar(fontsize = 10))
            }
        },

        # for clustering
        row_split = if (is.null(row_split)) sub(".*\\|", "", rownames(AUROC.scores)) else row_split,
        column_split = if (is.null(column_split)) sub(".*\\|", "", colnames(AUROC.scores)) else column_split,
        clustering_method_rows = "ward.D2",
        clustering_distance_rows = "spearman",
        clustering_method_columns = "ward.D2",
        clustering_distance_columns = "spearman",
        show_row_dend = TRUE, show_column_dend = TRUE,
        row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
        # row_names_side = 'left', column_names_side = 'top',
        rect_gp = gpar(col = "white", lwd = 1),
        width = ncol(AUROC.scores) * unit(5, "mm"),
        height = nrow(AUROC.scores) * unit(5, "mm"),

        # annotation
        row_labels = r_lab, column_labels = c_lab,
        bottom_annotation = c_anno,
        right_annotation = r_anno
    )

    return(draw(p, annotation_legend_list = list(lgd, lgd_sig), merge_legend = TRUE))
}
