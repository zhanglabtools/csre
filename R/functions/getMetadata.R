#' get metadata of csd or roadmap
#'
#' @param whichData "csd" or "roadmap"
#'
#' @return generate variables of metadata
#' @export
#'
#' @examples
getMetadata <- function(whichData = c("csd", "roadmap",
                                      "roadmap_1e4",
                                      "roadmap_1e5",
                                      "roadmap_1e6",
                                      "roadmap_200")) {
  whichData <- match.arg(whichData)
  if (whichData == "csd") {
    nameChrs <<- paste0("chr", c(1:22, "X"))
    nameSamples <<- c("H1", "K562", "GM12878", "HepG2", "Huvec",
                      "HSMM", "NHLF", "NHEK", "HMEC")
    nameMarks <<- c("CTCF", "H3K27ac", "H3K27me3",
                    "H3K36me3",
                    "H3K4me1", "H3K4me2", "H3K4me3",
                    "H3K9ac",
                    "H4K20me1",
                    "WCE")
    # numChrs <<- length(nameChrs)
    # numSamples <<- length(nameSamples)
    # numMarks <<- length(nameMarks)
    sortedNameMarks <<- c("CTCF",
                          "H3K27me3",
                          "H3K27ac", "H3K4me1", "H3K4me2",
                          "H3K4me3", "H3K9ac",
                          "H3K36me3",
                          "H4K20me1",
                          "WCE")
    ### seqInfo
    ## the below line would also load some pkgs automatically
    load("data/seqInfoHg18.Rdata")
    seqInfo <<- seqInfo
    # seqInfo <<- Seqinfo(genome = "hg18")
    sizeBin <<- 200L
    # sizeChrs <<- seqlengths(seqInfo)[nameChrs]
    # roundFun <- floor
    # numBins <<- roundFun(sizeChrs / sizeBin)
    thresNlp <<- 4
    maxNlpOfMarks <<- NA
    dirDataHome <<- "result/csd"
    dirNlp <<- "result/csd/nlp"
    dirFigure <<- "result/csd/figure"
    keepLastIncompleteBin <<- FALSE
    library(TxDb.Hsapiens.UCSC.hg18.refGene)
    txdb <<- TxDb.Hsapiens.UCSC.hg18.refGene
    library(BSgenome.Hsapiens.UCSC.hg18)
    bsg <<- BSgenome.Hsapiens.UCSC.hg18
  }
  if (whichData == "roadmap") {
    nameChrs <<- paste0("chr", c(1:22, "X"))
    nameSamples <<- paste0("E", sprintf("%03d", seq_len(129)[-c(60, 64)]))
    nameMarks <<- c("H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
    # numChrs <<- length(nameChrs)
    # numSamples <<- length(nameSamples)
    # numMarks <<- length(nameMarks)
    ### seqInfo
    ## the below line would also load some pkgs automatically
    load("data/seqInfoHg19.Rdata")
    seqInfo <<- seqInfo
    # seqInfo <<- Seqinfo(genome = "hg19")
    sizeBin <<- 25L
    # sizeChrs <<- seqlengths(seqInfo)[nameChrs]
    # roundFun <- ceiling
    # numBins <<- roundFun(sizeChrs / sizeBin)
    thresNlp <<- 2
    ## no truncation
    maxNlpOfMarks <<- NA
    dirDataHome <<- "result/roadmap"
    dirNlp <<- "result/roadmap/nlp"
    dirFigure <<- "result/roadmap/figure"
    keepLastIncompleteBin <<- TRUE
    library(TxDb.Hsapiens.UCSC.hg19.refGene)
    txdb <<- TxDb.Hsapiens.UCSC.hg19.refGene
    library(BSgenome.Hsapiens.UCSC.hg19)
    bsg <<- BSgenome.Hsapiens.UCSC.hg19
  }
  if (whichData == "roadmap_1e4") {
    nameChrs <<- paste0("chr", c(1:22, "X"))
    nameSamples <<- paste0("E", sprintf("%03d", seq_len(129)[-c(60, 64)]))
    nameMarks <<- c("H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
    # numChrs <<- length(nameChrs)
    # numSamples <<- length(nameSamples)
    # numMarks <<- length(nameMarks)
    ### seqInfo
    ## the below line would also load some pkgs automatically
    load("data/seqInfoHg19.Rdata")
    seqInfo <<- seqInfo
    # seqInfo <<- Seqinfo(genome = "hg19")
    sizeBin <<- 25L
    # sizeChrs <<- seqlengths(seqInfo)[nameChrs]
    # roundFun <- ceiling
    # numBins <<- roundFun(sizeChrs / sizeBin)
    thresNlp <<- 2
    ## 99.99%
    maxNlpOfMarks <<- c(23.93463, 160.14935, 19.04578, 17.83909, 13.52525)
    dirDataHome <<- "result/roadmap_1e4"
    dirNlp <<- "result/roadmap/nlp"
    dirFigure <<- "result/roadmap_1e4/figure"
    keepLastIncompleteBin <<- TRUE
    library(TxDb.Hsapiens.UCSC.hg19.refGene)
    txdb <<- TxDb.Hsapiens.UCSC.hg19.refGene
    library(BSgenome.Hsapiens.UCSC.hg19)
    bsg <<- BSgenome.Hsapiens.UCSC.hg19
  }
  if (whichData == "roadmap_1e5") {
    nameChrs <<- paste0("chr", c(1:22, "X"))
    nameSamples <<- paste0("E", sprintf("%03d", seq_len(129)[-c(60, 64)]))
    nameMarks <<- c("H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
    # numChrs <<- length(nameChrs)
    # numSamples <<- length(nameSamples)
    # numMarks <<- length(nameMarks)
    ### seqInfo
    ## the below line would also load some pkgs automatically
    load("data/seqInfoHg19.Rdata")
    seqInfo <<- seqInfo
    # seqInfo <<- Seqinfo(genome = "hg19")
    sizeBin <<- 25L
    # sizeChrs <<- seqlengths(seqInfo)[nameChrs]
    # roundFun <- ceiling
    # numBins <<- roundFun(sizeChrs / sizeBin)
    thresNlp <<- 2
    ## 99.999%
    maxNlpOfMarks <<- c(33.20958, 198.78615, 27.41964, 26.56604, 24.18820)
    dirDataHome <<- "result/roadmap_1e5"
    dirNlp <<- "result/roadmap/nlp"
    dirFigure <<- "result/roadmap_1e5/figure"
    keepLastIncompleteBin <<- TRUE
    library(TxDb.Hsapiens.UCSC.hg19.refGene)
    txdb <<- TxDb.Hsapiens.UCSC.hg19.refGene
    library(BSgenome.Hsapiens.UCSC.hg19)
    bsg <<- BSgenome.Hsapiens.UCSC.hg19
  }
  if (whichData == "roadmap_1e6") {
    nameChrs <<- paste0("chr", c(1:22, "X"))
    nameSamples <<- paste0("E", sprintf("%03d", seq_len(129)[-c(60, 64)]))
    nameMarks <<- c("H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
    # numChrs <<- length(nameChrs)
    # numSamples <<- length(nameSamples)
    # numMarks <<- length(nameMarks)
    ### seqInfo
    ## the below line would also load some pkgs automatically
    load("data/seqInfoHg19.Rdata")
    seqInfo <<- seqInfo
    # seqInfo <<- Seqinfo(genome = "hg19")
    sizeBin <<- 25L
    # sizeChrs <<- seqlengths(seqInfo)[nameChrs]
    # roundFun <- ceiling
    # numBins <<- roundFun(sizeChrs / sizeBin)
    thresNlp <<- 2
    ## no truncation
    # maxNlpOfMarks <<- NA
    ## 99.999%
    # maxNlpOfMarks <<- c(33.20958, 198.78615, 27.41964, 26.56604, 24.18820)
    ## 99.9999%
    maxNlpOfMarks <<- c(43.57203, 232.11694, 36.88951, 37.17780, 41.60989)
    dirDataHome <<- "result/roadmap_1e6"
    dirNlp <<- "result/roadmap/nlp"
    dirFigure <<- "result/roadmap_1e6/figure"
    keepLastIncompleteBin <<- TRUE
    library(TxDb.Hsapiens.UCSC.hg19.refGene)
    txdb <<- TxDb.Hsapiens.UCSC.hg19.refGene
    library(BSgenome.Hsapiens.UCSC.hg19)
    bsg <<- BSgenome.Hsapiens.UCSC.hg19
  }
  if (whichData == "roadmap_200") {
    nameChrs <<- paste0("chr", c(1:22, "X"))
    nameSamples <<- paste0("E", sprintf("%03d", seq_len(129)[-c(60, 64)]))
    nameMarks <<- c("H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
    # numChrs <<- length(nameChrs)
    # numSamples <<- length(nameSamples)
    # numMarks <<- length(nameMarks)
    ### seqInfo
    ## the below line would also load some pkgs automatically
    load("data/seqInfoHg19.Rdata")
    seqInfo <<- seqInfo
    # seqInfo <<- Seqinfo(genome = "hg19")
    sizeBin <<- 200L
    # sizeChrs <<- seqlengths(seqInfo)[nameChrs]
    # roundFun <- ceiling
    # numBins <<- roundFun(sizeChrs / sizeBin)
    thresNlp <<- 2
    ## no truncation
    maxNlpOfMarks <<- NA
    dirDataHome <<- "result/roadmap_200"
    dirNlp <<- "result/roadmap_200/nlp"
    dirFigure <<- "result/roadmap_200/figure"
    keepLastIncompleteBin <<- TRUE
    library(TxDb.Hsapiens.UCSC.hg19.refGene)
    txdb <<- TxDb.Hsapiens.UCSC.hg19.refGene
    library(BSgenome.Hsapiens.UCSC.hg19)
    bsg <<- BSgenome.Hsapiens.UCSC.hg19
  }
  invisible(NULL)
}


