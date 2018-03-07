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

#### figure 5C --------------------------------------------------------
grE003 <- subset(gr, gr$cellType == "E003")
oneE003 <- subsetByOverlaps(grE003, GRanges("chr22:28192342-28198164"))
zE003 <- oneE003$profile_zscore
grE009 <- subset(gr, gr$cellType == "E009")
oneE009 <- subsetByOverlaps(grE009, GRanges("chr22:28192342-28198164"))
zE009 <- oneE009$profile_zscore
grE082 <- subset(gr, gr$cellType == "E082")
oneE082 <- subsetByOverlaps(grE082, GRanges("chr22:28192342-28198164"))
zE082 <- oneE082$profile_zscore
grE070 <- subset(gr, gr$cellType == "E070")
oneE070 <- subsetByOverlaps(grE070, GRanges("chr22:28192342-28198164"))
zE070 <- oneE070$profile_zscore

df5 <- data.frame(EID = c("E003", "E009", "E082", "E070"),
           rbind(zE003, zE009, zE082, zE070)) %>%
  gather(mark, zscore, -EID) %>%
  mutate(mark = factor(mark, nameMarks),
         EID = factor(EID, c("E003", "E009", "E082", "E070")))
p5 <- ggplot(df5) +
  geom_bar(aes(x = mark, y = zscore, fill = mark),
           stat = "identity", show.legend = FALSE) +
  facet_wrap(~ EID, nrow = 1) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(size = rel(1.25),
                                   color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = rel(1.25),
                                   color = "black"),
        panel.grid.minor = element_blank()) +
  coord_flip() +
  ylab("z-score")
png(file.path(dirFigureCsre, "fig5_dynamic.png"),
    width = 5.5, height = 1.6, units = "in",
    res = 1200)
print(p5)
dev.off()

#### figure 5B --------------------------------------------------------
source("R/functions/getSpec.R")
rpkm <- readRDS("extData/roadmap/rpkm56.rds")
rpkm <- rpkm[rowMaxs(rpkm) >= 0.5, ]
logRpkm <- log2(rpkm + 1)
zLogRpkm <- getZscore(logRpkm)

symbolGene <- "MN1"
ensemblGene <- mapIds(org.Hs.eg.db, symbolGene,
                      "ENSEMBL", "SYMBOL")
expGene <- data.frame(EID = colnames(logRpkm),
                      expr = logRpkm[ensemblGene, ],
                      stringsAsFactors = FALSE) %>%
  left_join(ref, by = c("EID" = "EID"))
expGene$z <- zLogRpkm[ensemblGene, ]

#### figure 5B
ggExp <-
  ggplot(expGene, aes(x = 1, y = expr)) +
  geom_jitter(aes(color = GROUP), width = 0.25, show.legend = FALSE) +
  scale_color_manual(values = groupColors) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = rel(1.25))) +
  ylab(expression(paste(log[2], "(RPKM + 1)"))) +
  scale_x_continuous(limits = c(0.7, 1.3)) +
  scale_y_continuous(position = "right")

png(file.path(dirFigureCsre, "fig5_exp.png"),
    width = 4, height = 1.2, units = "in",
    res = 1200)
print(ggExp)
dev.off()


#### figure 6B --------------------------------------------------------

grE015 <- subset(gr, gr$cellType == "E015")
oneE015 <- subsetByOverlaps(grE015, GRanges("chr22:22890030-22890433"))
zE015 <- oneE015$profile_zscore
grE050 <- subset(gr, gr$cellType == "E050")
oneE050 <- subsetByOverlaps(grE050, GRanges("chr22:22890030-22890433"))
zE050 <- oneE050$profile_zscore
grE123 <- subset(gr, gr$cellType == "E123")
oneE123 <- subsetByOverlaps(grE123, GRanges("chr22:22890030-22890433"))
zE123 <- oneE123$profile_zscore

grNE123 <- subset(gr, gr$cellType != "E123")
grE123ExtremeSpec <- grE123[-queryHits(findOverlaps(grE123, grNE123))]

top <- grE123ExtremeSpec[order(grE123ExtremeSpec$score_pd_sh, decreasing = T)]

df6 <- data.frame(EID = c("E015", "E050", "E123"),
                  rbind(zE015, zE050, zE123)) %>%
  gather(mark, zscore, -EID) %>%
  mutate(mark = factor(mark, nameMarks),
         EID = factor(EID, c("E015", "E050", "E123")))
p6 <- ggplot(df6) +
  geom_bar(aes(x = mark, y = zscore, fill = mark),
           stat = "identity", show.legend = FALSE) +
  facet_wrap(~ EID, nrow = 1) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(size = rel(1.25),
                                   color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = rel(1.25),
                                   color = "black"),
        panel.grid.minor = element_blank()) +
  coord_flip() +
  ylab("z-score")
png(file.path(dirFigureCsre, "fig6_dynamic.png"),
    width = 4.5, height = 1.45, units = "in",
    res = 1200)
print(p6)
dev.off()

#### figure 6C --------------------------------------------------------

symbolGene <- "PRAME"
ensemblGene <- mapIds(org.Hs.eg.db, symbolGene,
                      "ENSEMBL", "SYMBOL")
expGene <- data.frame(EID = colnames(logRpkm),
                      expr = logRpkm[ensemblGene, ],
                      stringsAsFactors = FALSE) %>%
  left_join(ref, by = c("EID" = "EID"))
expGene$z <- zLogRpkm[ensemblGene, ]

#### figure 5B
ggExp <-
  ggplot(expGene, aes(x = 1, y = expr)) +
  geom_jitter(aes(color = GROUP), width = 0.25, show.legend = FALSE) +
  scale_color_manual(values = groupColors) +
  # coord_flip() +
  theme_classic(base_size = 8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25))) +
  ylab(expression(paste(log[2], "(RPKM + 1)"))) +
  scale_x_continuous(limits = c(0.7, 1.3)) +
  scale_y_continuous(position = "left")

png(file.path(dirFigureCsre, "fig6_exp.png"),
    width = 1.2, height = 1.8, units = "in",
    res = 1200)
print(ggExp)
dev.off()

