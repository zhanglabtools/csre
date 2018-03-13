library(magrittr)
library(tidyverse)
library(matrixStats)
library(ggpubr)
library(ggtree)
library(pheatmap)
library(RColorBrewer)
library(GenomicRanges)
library(rtracklayer)
source("R/functions/getMetadata.R")
whichData <- "roadmap"
getMetadata(whichData)

# load signal
regionHeight <- 0.1
regionWidth <- 60
fileCsre <- file.path(dirDataHome, "csre",
                      paste0("csre_nlp_c_z_ssz_2n_z_pd_sh_",
                             regionHeight, "_", regionWidth, ".txt"))
source("R/roadmap/getRoadmapMetadata.R")

nameCsre <- paste("csre_s", regionHeight, regionWidth, sep = "_")
dirFigureCsre <- file.path(dirFigure, nameCsre)
if (!file.exists(dirFigureCsre)) dir.create(dirFigureCsre, recursive = TRUE)
dirGr <- file.path(dirDataHome, "gr")
gr <- readRDS(file.path(dirGr, paste0("gr_", regionHeight, "_",
                                      regionWidth, ".rds")))

grl <- splitByCellType(gr)

disLen <- data.frame(chr = seqnames(gr),
                     cellType = gr$cellType,
                     length = width(gr))
lenCsreOfCellType <- sum(width(grl))
lenChrCellType <- disLen %>%
  group_by(chr, cellType) %>%
  summarise(length = sum(length)) %>%
  ungroup() %>%
  mutate(lenChr = sizeChrs[chr],
         fc = length / lenChr / (lenCsreOfCellType[cellType] / lenAllChrs)) %>%
  mutate(cellType = as.character(cellType)) %>%
  left_join(ref, by = c("cellType" = "EID"))

matLenChrCellType <-
  reshape2::acast(lenChrCellType, chr ~ cellType, value.var = "fc")

colnames(matLenChrCellType) <- ref[match(colnames(matLenChrCellType), ref$EID), "SHORT_NAME"]
annRow <- ref[match(colnames(matLenChrCellType), ref$SHORT_NAME), c("GROUP", "TYPE", "SEX")]
rownames(annRow) <- colnames(matLenChrCellType)
annRow$Log10Length <- log10(sum(width(grl)))
# color <- colorRampPalette((c("green","black", "red")))(100)
color <- colorRampPalette(rev(brewer.pal(n = 10, name = "Spectral")))(100)

#### draw_colnames_old is used to restore original pheatmap
draw_colnames_old <- function(coln, gaps, ...) {
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") -
                   unit(3, "bigpts"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
  return(res)
}
## set enclosing env to let it find "find_coordinates", "textGrob" and so on
## when is it called by pheatmap
environment(draw_colnames_old) <- asNamespace("pheatmap")
#### draw_colnames_angle is used to set preferred angle
draw_colnames_angle <- function(coln, gaps, ...) {
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") -
                   unit(3, "bigpts"),
                 vjust = 0.7, hjust = 1,
                 rot = 90, gp = gpar(...))
  return(res)
}
environment(draw_colnames_angle) <- asNamespace("pheatmap")
#### use preferred angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_angle,
                  ns = asNamespace("pheatmap"))
#### use original angle
# assignInNamespace(x = "draw_colnames",
#                   value = draw_colnames_old,
#                   ns = asNamespace("pheatmap"))

phtChr <-
  pheatmap(log2(matLenChrCellType + 1),
           # clustering_distance_cols = "correlation",
           # clustering_distance_rows = "correlation",
           clustering_method = "ward.D2",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           treeheight_row = 16,
           annotation_col = annRow,
           annotation_colors = list(GROUP = groupColors,
                                    TYPE = typeColors),
           color = color,
           fontsize = 6)$gtable
png(file.path(dirFigureCsre, "figS4_chr.png"),
    width = 12, height = 6, units = "in",
    res = 1200)
grid::grid.draw(phtChr)
dev.off()
#### use original angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_old,
                  ns = asNamespace("pheatmap"))
