library(GenomicRanges)
library(rtracklayer)
library(rhdf5)
# seqInfo <- Seqinfo(genome = "hg19")
seqInfo <- readRDS("data/seqInfoHg19.rds")

nameMarks <- c("H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
nameSamples <- paste0("E", sprintf("%03d", seq_len(129)[-c(60, 64)]))
nameChrs <- paste0("chr", c(1:22, "X"))

numMarks <- length(nameMarks)
numSamples <- length(nameSamples)
numChrs <- length(nameChrs)

# 25L is used in our paper. If you just want to test the pipeline, I recommend you use 200L or some integer larger.
sizeBin <- 25L
sizeChrs <- seqlengths(seqInfo)[nameChrs]
numBins <- ceiling(sizeChrs / sizeBin)

dirBw <- "E:/Users/cwang/roadmap/bigwig" # change to a dir containing bigwig files of -log10(p-value)
# dirBw <- "/share_bio/nas5/amsszhangsh_group/wangcan/data/roadmap/bigwig"
dirNlp <- "result/roadmap/nlp"
if (!file.exists(dirNlp)) dir.create(dirNlp, recursive = TRUE)

chunk <- c(10 ^ 5, 1, 1)
level <- 0

args <- commandArgs(trailingOnly = TRUE)
whichChrs <- as.integer(args)
stopifnot(length(setdiff(whichChrs, seq_len(numChrs))) == 0)

#' get gr of binned genome
#'
#' @param chr the chromosome to be tiled (binned)
#' @param sizeBin size of bin
#' @param seqInfo a Seqinfo
#'
#' @return a gr of binned genome
#' @export
#'
#' @examples
getTile <- function(chr, sizeBin, seqInfo) {
  ### 20x faster than tileGenome() and 8x faster than slidingWindows() +++++++++
  library(GenomicRanges)
  chrLen <- seqlengths(seqInfo)[[chr]]
  starts <- seq(1, chrLen, by = sizeBin)
  ends <- starts + sizeBin - 1L
  ends[length(ends)] <- min(ends[length(ends)], chrLen)
  chrTile <- GRanges(chr, IRanges(starts, ends))
}

tms <- Sys.time()
for (iChr in whichChrs) {
  chr <- nameChrs[iChr]
  cat(chr, "\n", sep = "")
  chrTile <- getTile(chr, sizeBin, seqInfo)
  ### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  chrFile <- file.path(dirNlp, paste0("nlp_", chr, ".h5"))
  numBinsInChr <- unname(numBins[chr])
  h5createFile(chrFile)
  h5createDataset(file = chrFile, dataset = chr,
                  dims = c(numBinsInChr, numSamples, numMarks),
                  storage.mode = "double",
                  chunk = chunk, level = level)
  for (iMark in seq_len(numMarks)) {
    mark <- nameMarks[iMark]
    cat(mark, "\n", sep = "")
    for (iSample in seq_len(numSamples)) {
      sample <- nameSamples[iSample]
      cat(sample, " ", sep = "")
      basenameBw <- paste0(sample, "-", mark, ".pval.signal.bigwig")
      fileName <- file.path(dirBw, basenameBw)
      bwf <- BigWigFile(fileName)
      stopifnot(chr %in% seqnames(seqinfo(bwf)))
      rlBw <- import.bw(fileName, as = "Rle", which = chrTile)[chr]
      ba <- binnedAverage(chrTile, rlBw, "mean")$mean
      stopifnot(length(ba) == numBinsInChr)
      h5write(obj = ba, file = chrFile, name = chr,
              index = list(NULL, iSample, iMark))
      # tmp <- h5read(chrFile, chr,
      #               index = list(NULL, iSample, iMark))
      # dim(tmp) <- NULL
      # stopifnot(identical(tmp, ba))
    }
    cat("\n")
  }
  cat("metadata")
  h5write(obj = c("bin", "sample", "mark"), file = chrFile, name = "dims",
          level = 0)
  h5write(obj = nameSamples, file = chrFile, name = "nameSamples",
          level = 0)
  h5write(obj = nameMarks, file = chrFile, name = "nameMarks",
          level = 0)
  cat("\n")
}
tme <- Sys.time()
cat("Total time:", as.numeric(tme - tms, units = "mins"), "mins\n")


