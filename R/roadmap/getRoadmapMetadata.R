## should define fileCsre before it is sourced

library(tidyverse)
library(GenomicRanges)
source("R/functions/csreGetAndConvert.R")

## load signal
fileCsre <- "result/roadmap/csre/csre_nlp_c_z_ssz_2n_z_pd_sh_0.1_60.txt"
csreSignal <-
  getCsreGRFromFile(fileCsre,
                    "df", nameChrs, nameSamples)
csreSignalGR <-
  getCsreGRFromFile(fileCsre,
                    "gr", nameChrs, nameSamples, seqInfo)

## change name to standard name
csre <- csreSignal
csreGR <- csreSignalGR

## additional gr or grl may be used
# csreGRL <- splitByCellType(csreGR)
# csreBinGR <- bpToBin(csreGR, sizeBin, keepLastIncompleteBin)
# csreBinGRL <- splitByCellType(csreBinGR)

## metadata
sizeChrs <- seqlengths(seqInfo)[nameChrs]
lenAllChrs <- sum(as.double(sizeChrs))
# source("R/functions/getFeatures.R")
# geneticGRL <- getGeneticFeatures(txdb = txdb,
#                                  of = "NM")
geneticGRL <- readRDS("data/geneticGRLHg19Refseq.rds")

nameSamplesENCODE <- paste0("E", 114:129)
nameSamplesRoadmap <- setdiff(nameSamples, nameSamplesENCODE)
nameSamplesCsd <- paste0("E", c("003", "123", "116", "118", "122",
                                "120", "128", "127", "119"))
sortedEID <- readRDS("data/sortedEID.rds")
shortName <- readRDS("data/shortName.rds")
eid <- data.table::fread("data/EID_metadata.tab", data.table = FALSE)

factorByOrder <- function(f) factor(f, unique(f))
ref <- eid[match(sortedEID, eid$EID), ] %>%
  mutate(GROUP = factorByOrder(GROUP),
         MNEMONIC = factorByOrder(MNEMONIC),
         TYPE = factorByOrder(TYPE),
         STD_NAME = factorByOrder(STD_NAME),
         EDACC_NAME = factorByOrder(EDACC_NAME),
         FULL_NAME = paste(EID, STD_NAME),
         FULL_NAME = factorByOrder(FULL_NAME)) %>%
  left_join(shortName, by = "EID") %>%
  mutate(SHORT_NAME = factorByOrder(SHORT_NAME))
groupColors <- unique(eid[, c("COLOR", "GROUP")])
groupColors <- setNames(groupColors$COLOR, groupColors$GROUP)
typeColors <- c(CellLine = "#84B3D7",
                PrimaryCulture = "#90D8BF",
                ESCDerived = "#C0B5DE",
                PrimaryCell = "#F9F8A2",
                PrimaryTissue = "#F36467")
mnemonicColors <- unique(eid[, c("COLOR", "MNEMONIC")])
mnemonicColors <- setNames(mnemonicColors$COLOR, mnemonicColors$MNEMONIC)
eidColors <- unique(eid[, c("COLOR", "EID")])
eidColors <- setNames(eidColors$COLOR, eidColors$EID)
shortNameColors <- unique(ref[, c("COLOR", "SHORT_NAME")])
shortNameColors <- setNames(shortNameColors$COLOR, shortNameColors$SHORT_NAME)
