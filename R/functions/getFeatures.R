#' get a GRL of genetic features from a txdb
#'
#' @param txdb the txdb from which to derive genetic features
#' @param upstream upstream in promoters
#' @param downstream downstream in promoters
#' @param of get genetic features from all tx or NM or NR
#'
#' @return a GRL of features
#' @export
#'
#' @examples geneticGRL <- getGeneticFeatures(TxDb.Hsapiens.UCSC.hg18.refGene)
getGeneticFeatures <- function(txdb,
                               upstream = 2000L, downstream = 2000L,
                               of = c("allTx", "NM", "NR")) {
  # if (missing(txdb)) library(TxDb.Hsapiens.UCSC.hg18.refGene)
  of <- match.arg(of)
  library(GenomicFeatures)
  tx <- transcripts(txdb)
  if (of == "allTx") {
    id <- seq_along(tx)
  } else if (of == "NM") {
    id <- grep("^NM_", tx$tx_name)
  } else {
    id <- grep("^NR_", tx$tx_name)
  }
  # gn <- genes(txdb)
  # txbg <- transcriptsBy(txdb, "gene")
  #### promoters
  pm <- trim(promoters(tx, upstream = 2000, downstream = 2000))
  promoterBp <- GenomicRanges::reduce(pm[id], ignore.strand = TRUE)
  #### 5'UTR
  fiveUTRbt <- fiveUTRsByTranscript(txdb)
  ## fiveUTRbt belong to only NM
  fiveUTRBp <- GenomicRanges::reduce(unlist(fiveUTRbt[names(fiveUTRbt) %in% id]),
                                     ignore.strand = TRUE)
  ## 3'UTR
  threeUTRbt <- threeUTRsByTranscript(txdb)
  ## threeUTRbt belong to only NM
  threeUTRBp <-
    GenomicRanges::reduce(unlist(threeUTRbt[names(threeUTRbt) %in% id]),
                          ignore.strand = TRUE)
  #### exons
  # ex <- exons(txdb)
  exbt <- exonsBy(txdb, "tx")
  ## exbt belong to all tx
  # exonBp <- GenomicRanges::reduce(ex, ignore.strand = TRUE)
  ## the above and below are identical
  exonBp <-
    GenomicRanges::reduce(unlist(exbt[names(exbt) %in% id]),
                          ignore.strand = TRUE)
  #### introns
  intrbt <- intronsByTranscript(txdb)
  ## intrbt belong to all tx
  intronBp <-
    GenomicRanges::reduce(unlist(intrbt[names(intrbt) %in% id]),
                          ignore.strand = TRUE)
  #### intergenic
  inge <-
    gaps(Reduce(GenomicRanges::union,
                list(promoterBp, fiveUTRBp, threeUTRBp,
                     exonBp, intronBp)))
  intergenicBp <- inge[strand(inge) == "*"]
  return(GRangesList(promoterBp = promoterBp,
                     fiveUTRBp = fiveUTRBp,
                     threeUTRBp = threeUTRBp,
                     exonBp = exonBp,
                     intronBp = intronBp,
                     intergenicBp = intergenicBp))
}

#' get GC and CpG stats from a BSgenome
#'
#' @param bsg a BSgenome
#' @param nameChrs name of chromosomes
#'
#' @return a list of totalGC, totalCpG, totalLn, totalNoN, bgGCFreq and
#'   bgCpGFreq
#' @export
#'
#' @examples cgGlobal <- getGCCpGGlobal(BSgenome.Hsapiens.UCSC.hg18)
getGCCpGGlobal <- function(bsg,
                           nameChrs = paste0("chr", c(1:22, "X"))) {
  # if (missing(bsg)) library(BSgenome.Hsapiens.UCSC.hg18)
  totalGC <- 0 ## total length of G|Cs
  totalCpG <- 0 ## total number of CpGs
  totalLn <- 0 ## total length of chrs
  totalNoN <- 0 ## total length of non N bps
  for (chr in nameChrs) {
    cat(chr, "")
    seqChr <- bsg[[chr]]
    totalGC <- totalGC + letterFrequency(seqChr, letters = "GC")
    totalCpG <- totalCpG + countPattern("CG", seqChr)
    totalLn <- totalLn + length(seqChr)
    totalNoN <- totalNoN + length(seqChr) - letterFrequency(seqChr, letters = "N")
  }
  bgGCFreq <- totalGC / totalNoN
  bgCpGFreq <- unname(totalCpG * 2 / totalNoN)
  return(list(totalGC = unname(totalGC), totalCpG = totalCpG,
              totalLn = totalLn, totalNoN = unname(totalNoN),
              bgGCFreq = unname(bgGCFreq), bgCpGFreq = bgCpGFreq))
}

#' get tss from a txdb
#'
#' @param txdb a txdb
#' @param of get tss of all tx or NM or NR
#'
#' @return a GRanges of tss
#' @export
#'
#' @examples
getTss <- function(txdb, of = c("allTx", "NM", "NR")) {
  # if (missing(txdb)) library(TxDb.Hsapiens.UCSC.hg18.refGene)
  of <- match.arg(of)
  tx <- transcripts(txdb)
  if (of == "allTx") {
    tss <- resize(tx, 1)
  } else if (of == "NM") {
    NM <- tx[grep("^NM_", tx$tx_name), ]
    tss <- resize(NM, 1)
  } else {
    NR <- tx[grep("^NR_", tx$tx_name), ]
    tss <- resize(NR, 1)
  }
  return(tss)
}
