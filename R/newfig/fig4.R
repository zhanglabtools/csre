library(magrittr)
library(tidyverse)
library(matrixStats)
library(LOLA)
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


resNb <- readRDS(file.path(dirFigureCsre,
                           paste0("resNb_",
                                  regionHeight,
                                  "_",
                                  regionWidth,
                                  ".rds")))
curCellType <- "E082"
cluster <- resNb[[curCellType]][["Best.partition"]]
grCur <- subset(gr, mcols(gr)[["cellType"]] == curCellType)
pz <- grCur$profile_zscore
pn <- grCur$profile_nlp

#### figure 4A --------------------------------------------------------
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
draw_colnames_angle <- function(coln, gaps, ...) {
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") -
                   unit(3, "bigpts"),
                 vjust = 0.5, hjust = 1,
                 rot = 90, gp = gpar(...))
  return(res)
}
environment(draw_colnames_angle) <- asNamespace("pheatmap")
#### use preferred angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_angle,
                  ns = asNamespace("pheatmap"))
getColorAndBreaks <- function(pz) {
  library(RColorBrewer)
  if (min(pz) < 0) {
    color <- colorRampPalette(c(rev(brewer.pal(n = 7, name = "Blues")),
                                brewer.pal(n = 7, name = "Reds")))(100)
    breaks <- unique(c(seq(min(pz), 0, length.out = 50),
                       seq(0, quantile(pz, 0.998), length.out = 50),
                       max(pz)))
  } else {
    color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    breaks <- unname(c(seq(min(pz), quantile(pz, 0.95), length.out = 99),
                       quantile(pz, c(0.975, 1))))
  }
  return(list(color = color, breaks = breaks))
}
clBr <- getColorAndBreaks(pz)
color <- clBr$color
breaks <- clBr$breaks

htPz <- pheatmap(pz[order(cluster), ],
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 color = color,
                 breaks = breaks,
                 fontsize = 8)$gtable
# htPn <- pheatmap(pmin(pn[order(cluster), ], 20),
#                  cluster_cols = FALSE,
#                  cluster_rows = FALSE,
#                  color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(100),
#                  fontsize = 8)$gtable
png(file.path(dirFigureCsre, "fig4_heatmap_zscore.png"),
    width = 1.6, height = 6, units = "in",
    res = 1200)
grid::grid.draw(htPz)
dev.off()
#### use original angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_old,
                  ns = asNamespace("pheatmap"))

#### figure 4B --------------------------------------------------------
dfClusterPz <- aggregate(pz, list(cluster = cluster), median) %>%
  gather(mark, zscore, -cluster) %>%
  mutate(mark = factor(mark, nameMarks))
plotPzMedian <- ggplot(dfClusterPz) +
  geom_bar(aes(x = mark, y = zscore, fill = mark),
           stat = "identity", show.legend = FALSE) +
  # geom_hline(yintercept = 0, linetype = 2, color = 2) +
  facet_grid(cluster ~ .) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.37,
                                   size = rel(1.25),
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = rel(1.25),
                                   color = "black"),
        panel.grid.minor = element_blank()) +
  ylab("z-score")
png(file.path(dirFigureCsre, "fig4_prof_z.png"),
    width = 1.6, height = 6, units = "in",
    res = 1200)
print(plotPzMedian)
dev.off()

#### figure 4C --------------------------------------------------------

chm <- import.bed("extData/roadmap/E082_15_coreMarks_mnemonics.bed.gz")
DNase <- readRDS("extData/roadmap/grlDNase.rds")[[curCellType]]
H3K27ac <- readRDS("extData/roadmap/grlH3K27ac.rds")[[curCellType]]
extFt <- c(DNase = DNase, H3K27ac = H3K27ac)
nameStates <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk",
                "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv",
                "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk",
                "15_Quies")
chm$name <- factor(chm$name, nameStates)

grlCur <- split(grCur, cluster)
source("R/functions/fc.R")
resExtFt <- fcGrlGrl(grlCur, GRangesList(extFt), lenAllChrs)
tpExtFt <- t(log2(resExtFt + 1))

chm <- split(chm, chm$name)
resChm <- fcGrlGrl(grlCur, GRangesList(chm), lenAllChrs)
tpResChm <- log2(resChm + 1)

grlBg <- readRDS("extData/roadmap/grlBg.rds")
regionDB <- loadRegionDB("extData/LOLA/chromHMM_E082")
resChmLola <- runLOLA(grlCur, grlBg[[curCellType]], regionDB)
oddsLola <- reshape2::acast(resChmLola, userSet ~ antibody, value.var = "oddsRatio")
nlpLola <- reshape2::acast(resChmLola, userSet ~ antibody, value.var = "pValueLog")
nlpLola[nlpLola == Inf] <- max(nlpLola[nlpLola != Inf])
tpResChm <- t(log2(oddsLola + 1))

breaks <- unique(c(seq(0, max(tpResChm), 1),
                   ceiling(max(tpResChm))))
color <- colorRampPalette(rev(brewer.pal(n = 8, name = "Spectral")))(length(breaks) - 1)

#### draw_colnames_angle is used to set preferred angle
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
phtChm <- pheatmap(tpResChm,
                   breaks = breaks,
                   color = color,
                   fontsize = 8,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   treeheight_row = 7,
                   treeheight_col = 7,
                   cellwidth = 9,
                   cellheight = 9,
                   legend = TRUE)$gtable
png(file.path(dirFigureCsre, "fig4_heatmap_chm.png"),
    width = 1.8, height = 2.3, units = "in",
    res = 1200)
grid::grid.draw(phtChm)
dev.off()
#### use original angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_old,
                  ns = asNamespace("pheatmap"))

breaks <- unique(c(seq(0, max(tpResChm), 1),
                   ceiling(max(tpResChm))))
color <- colorRampPalette(rev(brewer.pal(n = 10, name = "Spectral")))(length(breaks) - 1)
# color <- colorRampPalette(rev(brewer.pal(n = 10, name = "Spectral")))(100)
rn <- c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk",
        "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv",
        "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk",
        "Quies")
gghtChm <-
  ggplot(resChmLola, aes(x = as.character(userSet), y = antibody)) +
  geom_tile(aes(fill = log2(oddsRatio + 1)),
            color = "grey70") +
  geom_text(aes(label = ifelse(qValue < 0.001, "*", "")),
            nudge_y = -0.25) +
  # geom_text(aes(label = gtools::stars.pval(qValue)),
  #           nudge_y = -0.3) +
  scale_fill_gradientn(colors = color) +
  scale_x_discrete(limits = as.character(seq_along(grlCur)),
                   expand = c(0, 0)) +
  scale_y_discrete(limits = rev(rn),
                   expand = c(0, 0),
                   position = "left") +
  theme_minimal(base_size = 8) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = rel(1)),
        axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = 'pt')) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4.5))
png(file.path(dirFigureCsre, "fig4_heatmap_chm_gg.png"),
    width = 1.8, height = 2.3, units = "in",
    res = 1200)
print(gghtChm)
dev.off()

#### figure 4D --------------------------------------------------------
source("R/functions/getSpec.R")
rpkm <- readRDS("extData/roadmap/rpkm56.rds")
rpkm <- rpkm[rowMaxs(rpkm) >= 0.5, ]
logRpkm <- log2(rpkm + 1)
zLogRpkm <- getZscore(logRpkm)
# zLogRpkm <- (rpkm - rowMedians(rpkm)) / rowSds(rpkm)

tx <- transcripts(txdb, filter = list(tx_chrom = nameChrs))
NM <- tx[grep("^NM_", tx$tx_name), ]
prNM <- promoters(NM, 2000, 2000)
NMWithPr <- punion(NM, prNM)
mcols(NMWithPr) <- mcols(NM)

grlCur <- split(grCur, resNb[[curCellType]][["Best.partition"]])

source("R/functions/getGenesForCsres.R")
# genesUsed <- getNearestGenesOld(grlCur, NM)
genesUsed <- getNearestGenes(grlCur, NM)
# genesUsed <- getOvGenes(grlCur, NM)

#### overlap of genes of each cluster
matOv <- pOvGenes(genesUsed)

dfOv <- reshape2::melt(matOv, value.name = "numGenes") %>%
  dplyr::rename(row = Var1, col = Var2)
ggOv <-
  ggplot(dfOv) +
  geom_tile(aes(x = row, y = col, fill = numGenes),
            color = "black", show.legend = FALSE) +
  geom_text(aes(x = row, y = col, label = numGenes),
            size = 25.4 / 72 * 8) + # unit of size here is mm
  scale_fill_gradientn(colours = c("white", rep("red", 2)),
                       name = "-log10(p-value)") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_text(color = "black",
                                   lineheight = 0.7,
                                   size = rel(1.25)),
        # axis.title.y = element_blank(),
        axis.ticks = element_blank()) +
  scale_y_continuous(position = "left",
                     trans = "reverse",
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "group", y = "group")
png(file.path(dirFigureCsre, "fig4_heatmap_ov_gg_new.png"),
    width = 1.7, height = 1.7, units = "in",
    res = 1200)
print(ggOv)
dev.off()

#### figure 4E --------------------------------------------------------
source("R/functions/GOGWAS.R")
listEnsOfGrl <- convertListGeneName(genesUsed,
                                    from = "REFSEQ",
                                    to = "ENSEMBL")
rpkmUsed <- zLogRpkm
expOfGrl <- lapply(listEnsOfGrl,
                   function(ens)
                     rpkmUsed[rownames(rpkm) %in% ens, curCellType])
expOfGrlRest <- lapply(listEnsOfGrl,
                       function(ens)
                         rpkmUsed[!rownames(rpkm) %in% ens, curCellType])
expOfGrlAndAll <-
  c(expOfGrl,
    list(others = rpkmUsed[!rownames(rpkm) %in%
                             unlist(listEnsOfGrl),
                           curCellType]))
expUsed <- expOfGrlAndAll
dfbox <- reshape2::melt(expUsed) %>%
  dplyr::rename(class = L1, zscore = value)
pv <- vapply(expUsed[-length(expUsed)], function(x)
  wilcox.test(x, expUsed[[length(expUsed)]])$p.value, double(1))
dfPv <- data.frame(class = names(pv), pvalue = pv) %>%
  mutate(sig = gtools::stars.pval(pvalue),
         ypos = vapply(expUsed[-length(expUsed)], median, double(1)) + 0.05)
plotBox <-
  ggplot(dfbox,
         aes(x = class, y = zscore)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(aes(fill = class),
               show.legend = FALSE,
               outlier.size = 0.2,
               varwidth = FALSE,
               notch = FALSE) +
  # geom_jitter(size = 0.2, width = 0.1) +
  geom_text(aes(x = class, y = ypos, label = sig), dfPv) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25)), # default size = rel(0.8)
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank()) +
  # geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_brewer(palette = "Set2") +
  xlab("group") +
  ylab(expression(paste("z-score of ", log[2], "(RPKM+1)")))
## theme_get()$text theme_get()$axis.text
png(file.path(dirFigureCsre, "fig4_box_exp_new.png"),
    width = 1.86, height = 1.9, units = "in",
    res = 1200)
print(plotBox)
dev.off()
wilcox.test(expUsed[[3]], expUsed[[5]])

#### figure 4F --------------------------------------------------------
listEtzOfGrl <- convertListGeneName(genesUsed,
                                    from = "REFSEQ",
                                    to = "ENTREZID")

resGO <- GOEnrichListGeneSet(listEtzOfGrl)
resGO <- resGO %>% dplyr::rename(cluster = cellType)

keptOntology <- c("BP")
minSizeTerm <- 5
maxSizeTerm <- 500
# maxSizeTerm <- 1000
# maxSizeTerm <- Inf
numTopGO <- 5
thresQGO <- 1e-2
thresFcNlp <- 0
keptGO <- resGO %>%
  filter(sizeTerm >= minSizeTerm,
         sizeTerm <= maxSizeTerm,
         ONTOLOGY %in% keptOntology,
         fcMultiplyNlp >= thresFcNlp,
         logqByCellType >= -log(thresQGO, base = 10)) %>%
  group_by(cluster) %>%
  top_n(n = numTopGO, wt = logqByCellType) %>%
  ungroup() %>%
  arrange(cluster, desc(logqByCellType))

longPlotGO <- resGO %>%
  filter(TERM %in% keptGO$TERM,
         cluster %in% keptGO$cluster) %>%
  mutate(TERM = factor(TERM, unique(keptGO$TERM))) %>%
  complete(cluster, TERM,
           fill = list(logqByCellType = 0,
                       logp = 0))

arGO <- reshape2::acast(longPlotGO, cluster ~ TERM, value.var = "logp")
color <- colorRampPalette(rev(brewer.pal(n = 10, name = "Spectral")))(10)
tpArGO <- t(arGO)
dGOCol <- dist(t(tpArGO))
dGOCol <- as.dist(1 - cor(tpArGO))
hcGOCol <- hclust(dGOCol, "ward.D2")
dGORow <- dist(tpArGO)
dGORow <- as.dist(1 - cor(t(tpArGO)))
hcGORow <- hclust(dGORow, "ward.D2")
#### draw_colnames_angle is used to set preferred angle
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
phtGO <-
  pheatmap(pmin(tpArGO, 10),
           color = color,
           cluster_rows = hcGORow,
           cluster_cols = FALSE,
           fontsize = 8,
           # breaks = breaks,
           # cluster_rows = TRUE,
           # cluster_cols = FALSE,
           treeheight_row = 16,
           treeheight_col = 4,
           cellwidth = 9,
           cellheight = 9,
           legend = TRUE)$gtable
png(file.path(dirFigureCsre, "fig4_heatmap_go_new.png"),
    width = 5.2, height = 2.1, units = "in",
    res = 1200)
grid::grid.draw(phtGO)
dev.off()
#### use original angle
assignInNamespace(x = "draw_colnames",
                  value = draw_colnames_old,
                  ns = asNamespace("pheatmap"))

