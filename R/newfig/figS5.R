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

allCate <- readRDS(file.path(dirFigureCsre, "allCate.rds"))
matCate <- allCate$matCate

dfCellTypeCate <- data.frame(cellType = decode(mcols(gr)[["cellType"]]),
                             matCate)

longCellTypeCate <- dfCellTypeCate %>%
  gather("category", "overlapped", -cellType, factor_key = TRUE) %>%
  mutate(category = factor(category, rev(levels(category)))) %>%
  within(levels(category) <-
           stringr::str_replace_all(levels(category),
                                    c("Pre$" = "",
                                      "five" = "5'",
                                      "three" = "3'"))) %>%
  mutate(cellType = as.character(cellType)) %>%
  left_join(ref, by = c("cellType" = "EID"))

longCellTypeCateToPlot <- longCellTypeCate

plotProp <-
  ggplot(longCellTypeCateToPlot) +
  geom_bar(aes(x = SHORT_NAME, y = ..count..,
               fill = category,
               weight = overlapped),
           position = "fill", color = "black",
           show.legend = TRUE) +
  scale_x_discrete(limits = rev(levels(longCellTypeCateToPlot$SHORT_NAME))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set3", direction = -1) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 2,
                             title.position = "top")) +
  labs(x = NULL, y = "Proportion",
       fill = "Genomic regions") +
  theme_bw(base_size = 5) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank(),
        legend.justification = "right",
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
  coord_flip()
png(file.path(dirFigureCsre, "figS5_proportion.png"),
    width = 3.5, height = 10, units = "in",
    res = 1200)
print(plotProp)
dev.off()


dfOv <- as.data.frame(allCate$matOv)
bgLenGenetic <- sum(as(width(geneticGRL), "NumericList"))
smOvLen <- cbind(cellType = mcols(gr)[["cellType"]],
                 totalLen = width(gr), dfOv) %>%
  gather("category", "length", -cellType, -totalLen,
         factor_key = TRUE) %>%
  group_by(cellType, category) %>%
  summarise(length = sum(length), totalLen = sum(totalLen)) %>%
  ungroup() %>%
  mutate(foldChange =
           length / totalLen /
           (bgLenGenetic[as.character(category)] / lenAllChrs)) %>%
  mutate(category =
           factor(category,
                  labels =
                    stringr::str_replace_all(levels(category),
                                             c("Bp$" = "",
                                               "five" = "5'",
                                               "three" = "3'")))) %>%
  mutate(cellType = as.character(cellType)) %>%
  left_join(ref, by = c("cellType" = "EID"))
smOvLenToPlot <- smOvLen
plotFc <- ggplot(smOvLenToPlot) +
  geom_boxplot(aes(x = GROUP, y = foldChange, fill = GROUP),
               fatten = 1,
               outlier.size = 0.2) +
  geom_hline(yintercept = 1, linetype = 2, col = "red") +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.37,
                                   color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.spacing.y = unit(0, units = "pt")) +
  scale_fill_manual(values = groupColors) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 4,
                             title.position = "top")) +
  ylab("Fold Change")
png(file.path(dirFigureCsre, "figS5_fc_gg.png"),
    width = 3.5, height = 10, units = "in",
    res = 1200)
print(plotFc)
dev.off()
