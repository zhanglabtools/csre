library(magrittr)
library(tidyverse)
library(ggpubr)
library(ggtree)
library(GenomicRanges)
library(pheatmap)
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


######## config pheatmap
#### draw_colnames_old is used to restore original pheatmap
draw_colnames_old <- function(coln, gaps, ...) {
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") -
                   unit(3, "bigpts"), vjust = 0.5, hjust = 0, rot = 270,
                 gp = gpar(...))
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
                 vjust = 0.5, hjust = 0,
                 rot = 305, gp = gpar(...))
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


#### GO heatmap
fileGO <- file.path(dirFigureCsre, "resGO.rds")
resGO <- readRDS(fileGO)
keptOntology <- c("BP")
minSizeTerm <- 5
maxSizeTerm <- 500
# maxSizeTerm <- 1000
# maxSizeTerm <- Inf
numTopGO <- 1
thresQGO <- 1e-6
thresFcNlp <- 0
keptGO <- resGO %>%
  filter(sizeTerm >= minSizeTerm,
         sizeTerm <= maxSizeTerm,
         ONTOLOGY %in% keptOntology,
         fcMultiplyNlp >= thresFcNlp,
         logqByCellType >= -log(thresQGO, base = 10)) %>%
  group_by(cellType) %>%
  top_n(n = numTopGO, wt = logqByCellType) %>%
  ungroup() %>%
  mutate(cellType = as.character(cellType)) %>%
  left_join(ref, by = c("cellType" = "EID")) %>%
  arrange(SHORT_NAME, desc(logqByCellType))

longPlotGO <- resGO %>%
  filter(TERM %in% keptGO$TERM,
         cellType %in% keptGO$cellType) %>%
  mutate(cellType = as.character(cellType)) %>%
  left_join(ref, by = c("cellType" = "EID")) %>%
  mutate(TERM = factor(TERM, unique(keptGO$TERM))) %>%
  complete(cellType, TERM,
           fill = list(logqByCellType = 0,
                       logp = 0))

arGO <- reshape2::acast(longPlotGO, SHORT_NAME ~ TERM, value.var = "logp")
annRowGO <- ref[match(rownames(arGO), ref$SHORT_NAME),
              c("GROUP"), drop = FALSE]
rownames(annRowGO) <- rownames(arGO)
library(RColorBrewer)
# dGO <- dist(arGO)
# hcGO <- hclust()
# rich.colors(100)
colorGO <- colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(15)
# colorGO <- c(colorGO, rep(colorGO[length(colorGO)], 60))
# colorGO <- colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(25)
# colorGO <- c(colorGO, rep(colorGO[length(colorGO)], 100))
dGOCol <- as.dist(1 - cor(arGO))
# dGOCol <- dist(t(arGO))
hcGOCol <- hclust(dGOCol, "ward.D2")
dGORow <- as.dist(1 - cor(t(arGO)))
# dGORow <- dist(arGO)
hcGORow <- hclust(dGORow, "ward.D2")
htGO <-
  pheatmap(pmin(t(arGO), 15),
           cluster_rows = hcGOCol,
           cluster_cols = hcGORow,
           # clustering_distance_rows = "correlation",
           # clustering_distance_cols = "correlation",
           # clustering_method = "ward.D2",
           # cutree_rows = 11,
           # cutree_cols = 11,
           color = colorGO,
           # color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(20),
           # color = rich.colors(100),
           # color = blue2green2red(20),
           # breaks = c(seq(0, 12, 0.5), max(arGO)),
           annotation_col = annRowGO,
           annotation_colors = list(GROUP = groupColors[levels(ref$GROUP)]),
           annotation_legend = TRUE,
           annotation_names_row = FALSE,
           annotation_names_col = FALSE,
           cellwidth = 6.4,
           cellheight = 6.4,
           treeheight_row = 20,
           treeheight_col = 24,
           # border_color = NA,
           fontsize = 6.4
           # legend_breaks = c(seq(0, 12, 0.5), max(arGO)),
           # legend_labels = c(seq(0, 12, 0.5), max(arGO)),
  )$gtable
png(file.path(dirFigureCsre, "fig2_GO_heatmap_new.png"),
    width = 13, height = 8, units = "in",
    res = 600)
grid::grid.draw(htGO)
dev.off()

#### GWAS heatmap
fileGwas <- file.path(dirFigureCsre, "resGwas.rds")
resGwas <- readRDS(fileGwas)
library(gwascat)
data(ebicat37)
gwasGR <- GRanges(ebicat37)
genome(gwasGR) <- genome(seqInfo)
seqlevels(gwasGR) <- sortSeqlevels(seqlevels(gwasGR))
seqinfo(gwasGR) <- intersect(seqInfo, seqinfo(gwasGR))
names(mcols(gwasGR)) <- make.names(names(mcols(gwasGR)))
gwasIdName <- c("PUBMEDID", "DISEASE.TRAIT")
minNumUniqueSNPsOfStudy <- 0
maxNumUniqueSNPsOfStudy <- Inf
minNumHits <- 0
thresQByPairs <- 1
thresQByCellType <- 0.1
thresP <- 1
thresFc <- 0
thresRecall <- 0
keptStudy <- resGwas %>%
  filter(numUniqueSNPsOfStudy >= minNumUniqueSNPsOfStudy,
         numUniqueSNPsOfStudy <= maxNumUniqueSNPsOfStudy,
         numHits >= minNumHits,
         logqByPairs >= -log10(thresQByPairs),
         logqCellType >= -log10(thresQByCellType),
         logp >= -log10(thresP),
         fc >= thresFc,
         recall >= thresRecall) %>%
  group_by(cellType) %>%
  arrange(cellType, desc(logqByPairs)) %>%
  unite_("study", gwasIdName, remove = FALSE) %>%
  ungroup(cellType) %>%
  mutate(cellType = as.character(cellType)) %>%
  left_join(ref, by = c("cellType" = "EID")) %>%
  arrange(SHORT_NAME, logqByPairs)

id <- rep(TRUE, nrow(resGwas))
for (idName in gwasIdName) {
  id <- id & resGwas[[idName]] %in% keptStudy[[idName]]
}
longPlotStudy <- resGwas[id, ] %>%
  unite_("study", gwasIdName, remove = FALSE) %>%
  filter(cellType %in% keptStudy$cellType,
         study %in% keptStudy$study) %>%
  mutate(study = factor(study, unique(keptStudy$study)),
         cellType = as.character(cellType)) %>%
  complete(cellType, study,
           fill = list(logqByPairs = 0,
                       logqCellType = 0,
                       logp = 0)) %>%
  mutate(cellType = as.character(cellType)) %>%
  left_join(ref, by = c("cellType" = "EID"))

arGwas <- reshape2::acast(longPlotStudy, SHORT_NAME ~ study, value.var = "logp")
colnames(arGwas) <- sub("_", " - ", colnames(arGwas))
annRowGwas <- ref[match(rownames(arGwas), ref$SHORT_NAME),
              c("GROUP"), drop = FALSE]
rownames(annRowGwas) <- rownames(arGwas)

colorGwas <- colorRampPalette(rev(brewer.pal(n = 5, name = "Spectral")))(5)
# colorGwas <- c(colorGwas, rep(colorGO[length(colorGO)], 1))
dGwasCol <- as.dist(1 - cor(arGwas))
hcGwasCol <- hclust(dGwasCol, "ward.D2")
dGwasRow <- as.dist(1 - cor(t(arGwas)))
hcGwasRow <- hclust(dGwasRow, "ward.D2")
htGwas <-
  pheatmap(pmin(t(arGwas), 5),
         cluster_rows = hcGwasCol,
         cluster_cols = hcGwasRow,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         color = colorGwas,
         # breaks = c(0, 1, 2, 3, 4, max(arGwas)),
         annotation_col = annRowGwas,
         annotation_colors = list(GROUP = groupColors),
         annotation_legend = TRUE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         cellwidth = 6.4,
         cellheight = 6.4,
         treeheight_row = 25,
         treeheight_col = 30,
         # border_color = NA,
         fontsize = 6.4
)$gtable
png(file.path(dirFigureCsre, "fig2_Gwas_heatmap.png"),
    width = 6, height = 6, units = "in",
    res = 600)
grid::grid.draw(htGwas)
dev.off()

# plotGWAS <-
#   ggplot(longPlotStudy) +
#   geom_tile(aes(x = SHORT_NAME, y = study, fill = logp),
#             show.legend = TRUE, color = "grey") +
#   scale_fill_gradientn(colours = c("white", rep("red", 2)),
#                        name = "-log10(p-value)") +
#   scale_y_discrete(position = "left",
#                    limits = rev(levels(longPlotStudy$study)),
#                    expand = c(0, 0),
#                    labels = function(x) stringr::str_wrap(x, width = 50)) +
#   scale_x_discrete(expand = c(0, 0)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90,
#                                    hjust = 1,
#                                    vjust = 0.37),
#         axis.text.y = element_text(lineheight = 0.7,
#                                    size = 8,
#                                    vjust = 0.35),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         legend.position = c(0, 0),
#         legend.justification = c(1, 1),
#         legend.box.margin = margin(2, 0, 0, 0),
#         legend.title.align = 1,
#         legend.direction = "horizontal") +
#   ggtitle("GWAS enrichment of CSREs")
# print(plotGWAS)

#### use original angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_old,
                  ns = asNamespace("pheatmap"))
