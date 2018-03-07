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

tx <- transcripts(txdb, filter = list(tx_chrom = nameChrs))
NM <- tx[grep("^NM_", tx$tx_name), ]
prNM <- promoters(NM, 2000, 2000)
NMWithPr <- punion(NM, prNM)
mcols(NMWithPr) <- mcols(NM)
source("R/functions/getOvEleAndLen.R")
# NMOverlappedByBody <- getOverlappedElements(gr, NM)
# NMAllCsre <- unlist(strsplit(NMOverlappedByBody, split = ","))
library(org.Hs.eg.db)
# refseq <- keys(org.Hs.eg.db, "REFSEQ")
# nameNM <- refseq[grep("^NM_", refseq)]
nameNM <- unique(mcols(NM)$tx_name)

source("R/functions/getGenesForCsres.R")

# genesByCellType <- getOvGenes(grl, NMWithPr)
genesByCellType <- getNearestGenes(grl, NM)

hkVersion <- 2013
if (hkVersion == 2011) {
  hk <- data.table::fread("data/hkGenes.txt", sep = "\t", header = TRUE,
                          data.table = FALSE)
  NMinHk <- grep("^NM_",
                 unique(unlist(strsplit(hk$RefSeq, split = " \\/\\/\\/ "))),
                 value = TRUE)
}
if (hkVersion == 2013) {
  hk <- data.table::fread("data/HK_genes_2013.txt", sep = "\t", header = FALSE,
                          data.table = FALSE)
  NMinHk <- grep("^NM_", unique(hk[, 2]), value = TRUE)
}

# NMinHk <- grep("^NM_",
#                unique(unlist(mapIds(org.Hs.eg.db, keys = hk$`Gene name`,
#                                     column = c("REFSEQ"),
#                                     keytype = "SYMBOL", multiVals = "list"))),
#                value = TRUE)
NMinHkAndHg <- NMinHk[NMinHk %in% nameNM]
NMinHg <- nameNM
genesByCellTypeOvHk <- lapply(genesByCellType,
                              function(x, y) intersect(x, y),
                              NMinHkAndHg)
enrichHk <- batchFisherTest(lengths(genesByCellTypeOvHk),
                            lengths(genesByCellType),
                            length(NMinHkAndHg),
                            length(NMinHg),
                            alternative = "less")
numCsres <- lengths(grl)
numOfOverlappingGenes <- lengths(genesByCellType)
longEnrichHk <- data.frame(enrichHk) %>%
  mutate(cellType = rownames(.), logp = -logp) %>%
  tidyr::gather(score, value, -cellType) %>%
  mutate(score = factor(score, levels = c("logp", "fc"))) %>%
  left_join(ref, by = c("cellType" = "EID")) %>%
  spread(score, value) %>%
  mutate(pv = 10^(logp),
         qv = p.adjust(pv, method = "fdr"),
         sig = ifelse(qv < 0.01, "*", "")) %>%
  mutate(numCsres = numCsres[match(cellType, names(numCsres))],
         numOfOverlappingGenes = numOfOverlappingGenes[match(cellType, names(numOfOverlappingGenes))]) %>%
  arrange(fc) %>%
  mutate(SHORT_NAME = factor(SHORT_NAME, unique(SHORT_NAME)))


plothk <- ggplot(longEnrichHk) +
  geom_bar(aes(x = SHORT_NAME, y = fc, fill = GROUP), stat = "identity") +
  geom_text(aes(x = SHORT_NAME, y = fc + 0.03, label = sig),
            nudge_x = -0.5) +
  scale_fill_manual(values = groupColors) +
  coord_flip() +
  geom_abline(intercept = 1, slope = 0, linetype = 2, color = "red") +
  theme_bw(base_size = 6) +
  theme(axis.text.x = element_text(color = "black",
                                   size = rel(1.25)),
        axis.text.y = element_text(color = "black",
                                   size = rel(1.25)),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank()) +
  ylab("Fold Change")

png(file.path(dirFigureCsre, "figS5_hk_new.png"),
    width = 6.6, height = 10, units = "in",
    res = 1200)
print(plothk)
dev.off()

# ggplot(longEnrichHk) +
#   geom_boxplot(aes(x = GROUP, y = fc, fill = GROUP)) +
#   scale_fill_manual(values = groupColors)
