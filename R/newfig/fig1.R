library(magrittr)
library(tidyverse)
library(ggpubr)
library(ggtree)
library(GenomicRanges)
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


#### number of csres
numCsres <- table(gr$cellType)
numCsres <- setNames(as.integer(numCsres), names(numCsres))
dfNumCsres <- ref %>%
  mutate(number = numCsres[EID])
jitterNumCsres <-
  ggplot(dfNumCsres) +
  geom_boxplot(aes(x = factor(1), y = number), outlier.shape = NA) +
  geom_jitter(aes(x = factor(1), y = number, color = GROUP),
              width = 0.4,
              size = 0.5,
              show.legend = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank()) +
  # scale_x_discrete(limits = rev(dfNumCsres$MNEMONIC),
  #                  labels = rev(dfNumCsres$STD_NAME)) +
  scale_color_manual(values = groupColors) +
  scale_y_log10()
png(file.path(dirFigureCsre, "fig1_jitterNumCsres.png"),
    width = 1.2, height = 2, units = "in",
    res = 1200)
print(jitterNumCsres)
dev.off()


#### length of csres
lengthCsres <- sum(width(grl))
dfLengthCsres <- ref %>%
  mutate(length = lengthCsres[EID])
jitterLengthCsres <-
  ggplot(dfLengthCsres) +
  geom_boxplot(aes(x = factor(1), y = (length / 1e6)), outlier.shape = NA) +
  geom_jitter(aes(x = factor(1), y = (length / 1e6), color = GROUP),
              width = 0.4,
              size = 0.5,
              show.legend = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank()) +
  # scale_x_discrete(limits = rev(dfNumCsres$MNEMONIC),
  #                  labels = rev(dfNumCsres$STD_NAME)) +
  scale_color_manual(values = groupColors) +
  scale_y_log10() +
  ylab("length (Mb)")
png(file.path(dirFigureCsre, "fig1_jitterLengthCsres.png"),
    width = 1.2, height = 2, units = "in",
    res = 1200)
print(jitterLengthCsres)
dev.off()

#### pie chart of genome coverage
bpCoverage <- colSums(table(coverage(gr)))
symCoverage <- tibble(coverage = names(bpCoverage),
                      numBps = bpCoverage)
dfCoverage <- data.frame(coverage = c("0", "1", "2", "3 or more"),
                         numBps = c(bpCoverage[1:3],
                                    "3 or more" = sum(bpCoverage[-(1:3)]))) %>%
  mutate(prop = numBps / sum(numBps),
         pos = 1 - (cumsum(prop) - prop / 2))

pieCoverage <-
  ggplot(dfCoverage) +
  geom_bar(aes(x = 1, y = prop, fill = coverage),
           stat = "identity", width = 1) +
  geom_text(aes(x = c(1.3, 1.3, 1.3, 1.6),
                y = pos,
                label = paste0(round(prop * 100), "%")),
            size = 25.4 / 72 * 8) +
  coord_polar(theta = "y") +
  theme_void(base_size = 8) +
  scale_fill_brewer(palette = "Set1",
                    name = "number of CSREs\n covering genome")
png(file.path(dirFigureCsre, "fig1_pieCoverage.png"),
    width = 3.5, height = 3, units = "in",
    res = 1200)
print(pieCoverage)
dev.off()

#### heatmap
ov <- findOverlaps(gr, gr)
csreGRRow <- gr[queryHits(ov)]
csreGRCol <- gr[subjectHits(ov)]
mcols(ov) <- DataFrame(cellTypeRow = csreGRRow$cellType,
                       cellTypeCol = csreGRCol$cellType,
                       idRow = queryHits(ov),
                       idCol = subjectHits(ov),
                       lengthRow = width(csreGRRow),
                       lengthCol = width(csreGRCol),
                       lengthOv = width(pintersect(csreGRRow, csreGRCol)))
len <- sum(width(grl))
lenAllChrs <- sum(as.numeric(sizeChrs))
fcOv <- data.frame(mcols(ov)) %>%
  group_by(cellTypeRow, cellTypeCol) %>%
  summarise(lengthOv = sum(lengthOv)) %>%
  mutate(lengthRow = len[cellTypeRow],
         lengthCol = len[cellTypeCol],
         lengthAll = lenAllChrs,
         fc = lengthOv * lengthAll / lengthRow / lengthCol,
         jacard = lengthOv / (lengthRow + lengthCol - lengthOv)) %>%
  ungroup() %>%
  mutate(cellTypeRow = factor(cellTypeRow, sortedEID),
         cellTypeCol = factor(cellTypeCol, sortedEID)) %>%
  tidyr::complete(cellTypeRow, cellTypeCol, fill = list(fc = 0, jacard = 0)) %>%
  mutate(fc = ifelse(cellTypeRow == cellTypeCol, NaN, fc))
matFc <- reshape2::acast(fcOv, cellTypeRow ~ cellTypeCol, value.var = "fc")
matFc[is.nan(matFc)] <- 0
rownames(matFc) <- colnames(matFc) <- ref$MNEMONIC[match(rownames(matFc), ref$EID)]
dd <- as.dist(max(log2(matFc + 1)) - log2(matFc + 1))
hc <- hclust(dd, method = "ward.D2")
sil <- sapply(2:126, function(x)
  mean(cluster::silhouette(cutree(hc, x), dd)[, 3]))
png(file.path(dirFigureCsre, "figS_tree_sil.png"),
    width = 6, height = 3, units = "in",
    res = 1200,
    pointsize = 8)
plot(2:126, sil, xlab = "Number of clusters", ylab = "Average Silhouette scores")
abline(v = 25, lty = 2, col = "red")
# abline(v = 39, lty = 2, col = "red")
dev.off()
ng <- cutree(hc, 25)
og <- as.integer(ref[match(names(ng), ref$MNEMONIC), "GROUP"])
# tissue <- sapply(strsplit(as.character(ref$MNEMONIC), "\\."), `[`, 1)
# og <- as.integer(ref[match(names(ng), ref$MNEMONIC), "GROUP"])
ct <- table(og, ng)
annRow <- data.frame(GROUP = levels(ref$GROUP))
rownames(annRow) <- rownames(ct)
draw_colnames_angle <- function(coln, gaps, ...) {
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") -
                   unit(3, "bigpts"),
                 vjust = 0.7, hjust = 0.5,
                 rot = 0, gp = gpar(...))
  return(res)
}
environment(draw_colnames_angle) <- asNamespace("pheatmap")
#### use preferred angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_angle,
                  ns = asNamespace("pheatmap"))
cont <- pheatmap(ct[1:19, ],
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         clustering_method = "complete",
         display_numbers = TRUE,
         number_format = "%d",
         annotation_row = annRow[1:19, , drop = FALSE],
         annotation_colors = list(GROUP = groupColors[levels(ref$GROUP)]),
         fontsize = 8,
         fontsize_number = 8,
         show_rownames = FALSE,
         annotation_names_row = FALSE)$gtable
png(file.path(dirFigureCsre, "figS_tree_cont_25.png"),
    width = 6, height = 3.5, units = "in",
    res = 1200,
    pointsize = 8)
grid::grid.draw(cont)
dev.off()
split(ref$MNEMONIC, ng) %>% lapply(as.character)

ari <- sapply(2:126, function(x) {
  phyclust::RRand(as.integer(as.factor(tissue)), cutree(hc, x))[["adjRand"]]
})
plot(2:126, ari)
table(og, ng)
library(pheatmap)
library(RColorBrewer)
# annRow <- ref[match(rownames(matFc), ref$MNEMONIC),
#               c("GROUP", "TYPE", "SEX", "ANATOMY")]
annRow <- ref[match(rownames(matFc), ref$MNEMONIC),
              c("GROUP"), drop = FALSE]
annCol <- ref[match(colnames(matFc), ref$MNEMONIC),
              c("GROUP"), drop = FALSE]
rownames(annRow) <- rownames(matFc)
rownames(annCol) <- colnames(matFc)
mat <- log2(matFc + 1)
diag(mat) <- max(mat) * 1

#### draw_colnames_old is used to restore original pheatmap
draw_colnames_old <- function(coln, gaps, ...) {
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"),
                 vjust = 0.5, hjust = 0,
                 rot = 270, gp = gpar(...))
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
                 vjust = 0.5, hjust = 1,
                 rot = 0, gp = gpar(...))
  return(res)
}
environment(draw_colnames_angle) <- asNamespace("pheatmap")
#### use preferred angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_angle,
                  ns = asNamespace("pheatmap"))
#### use original angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_old,
                  ns = asNamespace("pheatmap"))
colorMat <- colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100)
# colorMat <- c(colorMat, rep(colorMat[length(colorMat)], 100))
phtOv <- pheatmap(pmin(mat, 5),
         cluster_rows = hc,
         cluster_cols = hc,
         color = colorMat,
         # color = colorRampPalette(c(rep("white", 2), rep("yellow", 3), rep("red", 3)))(30),
         # color = colorRampPalette(rev(brewer.pal(n = 7, name =
         #                                           "Oranges")))(100),
         # kmeans_k = NA,
         # cutree_rows = 10,
         # clustering_distance_rows = "correlation",
         # clustering_method = "single",
         # gaps_row = which(diff(as.numeric(an$cluster)) == 1),
         # scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         # color = colorRampPalette(c("blue", "white", "yellow", "red"))(100),
         # breaks = c(seq(0, 200, 3), max(toPlot)),
         # annotation_row = annotation_row,
         annotation_row = annRow,
         annotation_col = annCol,
         annotation_colors = list(GROUP = groupColors[levels(ref$GROUP)]),
         # annotation_colors = ann_colors,
         # color = colorRampPalette(c("white","red"))(100),
         annotation_legend = TRUE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize = 8
)$gtable
png(file.path(dirFigureCsre, "fig1_heatmap.png"),
    width = 5.2, height = 4, units = "in",
    res = 1200)
grid::grid.draw(phtOv)
dev.off()

png(file.path(dirFigureCsre, "fig1_heatmap_legend.png"),
    width = 5.2, height = 4, units = "in",
    res = 1200)
pheatmap(pmin(mat, 5),
         cluster_rows = hc,
         cluster_cols = hc,
         color = colorMat,,
         # color = colorRampPalette(c(rep("white", 2), rep("yellow", 3), rep("red", 3)))(30),
         # color = colorRampPalette(rev(brewer.pal(n = 7, name =
         #                                           "Oranges")))(100),
         # kmeans_k = NA,
         # cutree_rows = 10,
         # clustering_distance_rows = "correlation",
         # clustering_method = "single",
         # gaps_row = which(diff(as.numeric(an$cluster)) == 1),
         # scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         # color = colorRampPalette(c("blue", "white", "yellow", "red"))(100),
         # breaks = c(seq(0, 200, 3), max(toPlot)),
         # annotation_row = annotation_row,
         annotation_row = annRow,
         # annotation_col = annCol,
         annotation_colors = list(GROUP = groupColors[levels(ref$GROUP)]),
         # annotation_colors = ann_colors,
         # color = colorRampPalette(c("white","red"))(100),
         annotation_legend = TRUE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize = 8
)
dev.off()

#### circle cluster showing names
matFc <- reshape2::acast(fcOv, cellTypeRow ~ cellTypeCol, value.var = "fc")
matFc[is.nan(matFc)] <- 0
rownames(matFc) <- colnames(matFc) <- ref$SHORT_NAME[match(rownames(matFc), ref$EID)]
dd <- as.dist(max(log2(matFc + 1)) - log2(matFc + 1))
hc <- hclust(dd, method = "ward.D2")
tree <- as.phylo(hc)
# treeWithMeta <- ggtree(tree, layout = "circular") %<+%
#   (ref %>% dplyr::select(MNEMONIC, everything()))
# treeWithMeta +
#   geom_tiplab(aes(angle = angle, color = GROUP), size = 2, offset = 0.1) +
#   geom_tippoint(aes(color = GROUP), size = 2) +
#   scale_color_manual(values = groupColors) +
#   xlim(0, 7)
ht <- ref
rownames(ht) <- ht$SHORT_NAME
circleTree <-
  gheatmap(ggtree(tree, layout = "circular") +
             geom_tiplab(aes(angle = angle), size = 2, offset = 1.8) +
             xlim(0, 35),
           ht[, c("GROUP"), drop = FALSE],
           offset = -0.2, width = 0.08, colnames = FALSE) +
  scale_fill_manual(values = groupColors) +
  # coord_polar(theta = "y") +
  theme(legend.position = "none")
png(file.path(dirFigureCsre, "fig1_circleTree.png"),
    width = 8000, height = 8000, res = 1200)
print(circleTree)
dev.off()
# ggsave("fig1_tree.png", p, "png", width = 6, height = 6, limitsize = FALSE)


#### config ggpubr
.test_pairwise_exact <- function(data, formula, method = "wilcox.test", ...)
{

  x <- deparse(formula[[2]])
  group <- attr(stats::terms(formula), "term.labels")

  # One sample test
  if(.is_empty(group)){
    res <- .test(data, formula, method = method,  ...)
    return(res)
  }

  # Pairwise test
  method <- switch(method,
                   t.test = "pairwise.t.test",
                   wilcox.test = "pairwise.wilcox.test")
  test <- match.fun(method)

  test.opts <- list(x = .select_vec(data, x),
                    g = .select_vec(data, group),  ...)
  # if(method == "pairwise.wilcox.test") test.opts$exact <- FALSE

  pvalues <- do.call(test, test.opts)$p.value %>%
    as.data.frame()
  group1 <- group2 <- p <- NULL
  pvalues$group2 <- rownames(pvalues)
  pvalues <- pvalues %>%
    tidyr::gather(key = "group1", value = "p", -group2) %>%
    dplyr::select(group1, group2, p) %>%
    dplyr::filter(!is.na(p))
  pvalues
}
environment(.test_pairwise_exact) <- asNamespace("ggpubr")
assignInNamespace(x = ".test_pairwise",
                  value = .test_pairwise_exact,
                  ns = asNamespace("ggpubr"))

stat_compare_means_exact <- function (mapping = NULL, data = NULL, method = NULL, paired = FALSE,
                                      ref.group = NULL, comparisons = NULL, hide.ns = FALSE, label.sep = ", ",
                                      label = NULL, label.x.npc = "left", label.y.npc = "top",
                                      label.x = NULL, label.y = NULL, geom = "text", position = "identity",
                                      na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...)
{
  if (!is.null(comparisons)) {
    method.info <- .method_info(method)
    method <- method.info$method
    test.args <- list(paired = paired)
    # if (method == "wilcox.test")
    #   test.args$exact <- FALSE
    pms <- list(...)
    size <- ifelse(is.null(pms$size), 0.3, pms$size)
    color <- ifelse(is.null(pms$color), "black", pms$color)
    map_signif_level <- FALSE
    if (.is_p.signif_in_mapping(mapping) | !.is_empty(label %in%
                                                      "p.signif")) {
      map_signif_level <- c(`****` = 1e-04, `***` = 0.001,
                            `**` = 0.01, `*` = 0.05)
      if (hide.ns)
        map_signif_level[4] <- " "
    }
    step_increase <- ifelse(is.null(label.y), 0.12, 0)
    ggsignif::geom_signif(comparisons = comparisons, y_position = label.y,
                          test = method, test.args = test.args, step_increase = step_increase,
                          size = size, color = color, map_signif_level = map_signif_level)
  }
  else {
    mapping <- .update_mapping(mapping, label)
    layer(stat = StatCompareMeans, data = data, mapping = mapping,
          geom = geom, position = position, show.legend = show.legend,
          inherit.aes = inherit.aes, params = list(label.x.npc = label.x.npc,
                                                   label.y.npc = label.y.npc, label.x = label.x,
                                                   label.y = label.y, label.sep = label.sep, method = method,
                                                   paired = paired, ref.group = ref.group, hide.ns = hide.ns,
                                                   na.rm = na.rm, ...))
  }
}
environment(stat_compare_means_exact) <- asNamespace("ggpubr")
unlockBinding("stat_compare_means", as.environment("package:ggpubr"))
assign("stat_compare_means", stat_compare_means_exact,
       envir = as.environment("package:ggpubr"))
lockBinding("stat_compare_means", as.environment("package:ggpubr"))
rm(stat_compare_means_exact)


#### fold change of chrX by sex
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


esMap <- setNames(nm = unique(as.character(ref$GROUP)))
esMap[!names(esMap) %in% c("ESC", "iPSC")] <- "Others"
sexMap <- setNames(nm = unique(as.character(ref$SEX)))
sexMap[!names(sexMap) %in% c("Female", "Male")] <- "Unsure"
dfXSexEs <-
  lenChrCellType %>%
  filter(chr == "chrX") %>%
  mutate(es = esMap[GROUP],
         sex = sexMap[SEX])

plotLog2FcOnChrXByEs <-
  ggplot(dfXSexEs, aes(x = sex, y = log2(fc), fill = sex)) +
  geom_boxplot(outlier.shape = NA, notch = FALSE,
               show.legend = FALSE) +
  # geom_violin() +
  geom_jitter(width = 0.15, size = 0.5,
              show.legend = FALSE) +
  # geom_hline(yintercept = 0, linetype = 2) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25),
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()) +
  ylab("log2(Fold Change)") +
  ylim(-8, 5)


scaleFun <- function(x) sprintf("%.1f", x)

SexComparisons <- list(c("Female", "Male"))
plotLog2FcOnChrXByEs <-
  ggboxplot(dfXSexEs, x = "sex", y = "log2(fc)",
            color = "sex",# palette = "jco",
            add = "jitter", #legend = "none",
            outlier.shape = NA,
            size = 0.5) +
  stat_compare_means(comparisons = SexComparisons) +
    # stat_compare_means() +
  ylab("log2(Fold Change)") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25),
                                   angle = 90,
                                   hjust = 1,
                                   vjust = 0.37),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        legend.position = "none") +
  scale_y_continuous(labels = scaleFun,
                     limits = c(NA, 5))
png(file.path(dirFigureCsre, "fig1_lenChrXSex.png"),
    width = 1.8, height = 2, units = "in",
    res = 1200)
print(plotLog2FcOnChrXByEs)
dev.off()

#### fold change of chrX of female by group
plotLog2FcOnChrXByGroup <-
  ggplot(dfXSexEs %>% filter(SEX == "Female"),
         aes(x = es, y = log2(fc), fill = es)) +
  geom_boxplot(outlier.shape = NA, notch = FALSE,
               show.legend = FALSE) +
  # geom_violin() +
  geom_jitter(width = 0.15, size = 0.5,
              show.legend = FALSE) +
  # geom_hline(yintercept = 0, linetype = 2) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25),
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()) +
  ylab("log2(Fold Change)") +
  xlab("groups") +
  scale_fill_manual(values = groupColors) +
  ylim(-3, 4.5)
esComparisons <- list(c("ESC", "Others"))
esColors <- c(groupColors[c("ESC", "iPSC")], "Others" = "#000000")
plotLog2FcOnChrXByGroup <-
  ggboxplot(dfXSexEs %>%
              filter(SEX == "Female") %>%
              mutate(es = factor(es, c("ESC", "iPSC", "Others"))),
            x = "es", y = "log2(fc)",
            color = "es",# palette = "jco",
            add = "jitter", #legend = "none",
            outlier.shape = NA) +
  stat_compare_means(comparisons = esComparisons,
                     method = "wilcox.test") +
  ylab("log2(Fold Change)") +
  xlab("groups") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.37),
        legend.position = "none") +
  scale_color_manual(values = esColors) +
  scale_y_continuous(labels = scaleFun,
                     limits = c(NA, 5))

png(file.path(dirFigureCsre, "fig1_lenChrXEsInSex.png"),
    width = 1.8, height = 2, units = "in",
    res = 1200)
print(plotLog2FcOnChrXByGroup)
dev.off()
