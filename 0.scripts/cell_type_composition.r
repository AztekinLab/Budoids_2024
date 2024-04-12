library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

cal_perc <- function(df, group, fill, enquo_var = TRUE) {
  if (enquo_var) {
    group <- enquo(group)
    fill <- enquo(fill)
  }

  cal_freq(df, group, fill, enquo_var = FALSE) %>%
    mutate(percent = n * 100 / sum(n))
}

cal_freq <- function(df, group, fill, enquo_var = TRUE) {
  if (enquo_var) {
    group <- enquo(group)
    fill <- enquo(fill)
  }
  df %>%
    dplyr::count(!!group, !!fill) %>%
    group_by(!!group)
}


composition_plot <- function(df, seu, group, fill, colours = NULL, palette = "Set1") {
  if (is.null(df)) {
    df <- seu@meta.data
  }
  group <- enquo(group)
  fill <- enquo(fill)
  cell_perc <- cal_perc(df, group, fill, enquo_var = FALSE)

  # cell_perc[,!!fill] <- factor(cell_perc[,!!fill], levels = levels(df[,!!fill]))

  p <- ggplot(cell_perc, aes(!!group, percent, fill = !!fill)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    theme_bw()

  if (!is.null(colours)) {
    p <- p + scale_fill_manual(values = colours)
  } else {
    p <- p + scale_fill_brewer(palette = palette)
  }

  return(p)
}

# USAGE
# p <- composition_plot(seu_list[[1]], CellType, Phase)
