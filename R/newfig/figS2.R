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

#### number of csres
numCsres <- table(gr$cellType)
numCsres <- setNames(as.integer(numCsres), names(numCsres))
dfNumCsres <- ref %>%
  mutate(number = numCsres[EID])
plotNumCsres <-
  ggplot(dfNumCsres) +
  geom_bar(aes(x = SHORT_NAME, y = number, fill = GROUP),
           stat = "identity",
           show.legend = FALSE) +
  scale_x_discrete(limits = rev(dfNumCsres$SHORT_NAME)) +
  scale_fill_manual(values = groupColors) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank()) +
  coord_flip()

#### base pairs of csres in each cell type
lengthCsres <- sum(width(grl))
dfLengthCsres <- ref %>%
  mutate(length = lengthCsres[EID])
plotLengthCsres <-
  ggplot(dfLengthCsres) +
  geom_bar(aes(x = SHORT_NAME, y = length / 1e6, fill = GROUP),
           stat = "identity",
           show.legend = FALSE) +
  scale_x_discrete(limits = rev(dfLengthCsres$SHORT_NAME)) +
  scale_fill_manual(values = groupColors) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() +
  ylab("length (Mb)")

#### boxplot of length of csres
dfBox <- reshape2::melt(as.list(width(grl))) %>%
  dplyr::rename(EID = L1, length = value) %>%
  mutate(EID = factor(EID, sortedEID),
         SHORT_NAME = ref[match(EID, ref$EID), "SHORT_NAME"],
         GROUP = ref[match(EID, ref$EID), "GROUP"])
plotBox <-
  ggplot(dfBox) +
  geom_boxplot(aes(x = SHORT_NAME, y = log10(length), fill = GROUP),
               show.legend = FALSE,
               fatten = 1,
               outlier.size = 0.2) +
  scale_x_discrete(limits = rev(ref$SHORT_NAME)) +
  scale_fill_manual(values = groupColors) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() +
  ylab(expression(log[10](length)))

#### number of neighbor genes
# txln <- transcriptLengths(txdb, with.cds_len = TRUE)
tx <- transcripts(txdb, filter = list(tx_chrom = nameChrs))
NM <- tx[grep("^NM_", tx$tx_name), ]
prNM <- promoters(NM, 2000, 2000)
NMWithPr <- punion(NM, prNM)
mcols(NMWithPr) <- mcols(NM)
numOfOverlappingGenes <- countOverlaps(grl, NMWithPr)
source("R/functions/getGenesForCsres.R")
# genesByCellType <- getOvGenes(grl, NMWithPr)
genesByCellType <- getNearestGenes(grl, NM)

dfNumOfOverlappingGenes <- ref %>%
  mutate(numGenes = numOfOverlappingGenes[EID])

plotNumGenes <-
  ggplot(dfNumOfOverlappingGenes) +
  geom_bar(aes(x = SHORT_NAME, y = numGenes, fill = GROUP),
           stat = "identity",
           show.legend = TRUE) +
  scale_x_discrete(limits = rev(dfNumOfOverlappingGenes$SHORT_NAME)) +
  scale_fill_manual(values = groupColors) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() +
  ylab("#neighboring genes")


library(egg)
ggTotal <-
  ggarrange(plotNumCsres, plotLengthCsres,
            plotBox, plotNumGenes,
            nrow = 1)
png(file.path(dirFigureCsre, "figS2_basicStats.png"),
    width = 7, height = 9, units = "in",
    res = 1200)
grid::grid.draw(ggTotal)
dev.off()
