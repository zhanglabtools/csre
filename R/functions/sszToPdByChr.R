#' Get pd from ssz by normalize, zscore and prod by chr
#'
#' @param chr chromosome
#' @param fSszInChr ssz file in chr
#' @param nameSamples name of samples
#' @param nameMarks name of marks
#' @param factorTwoNorm the normalization factor to be multiplied to ssz of each
#'   sample
#' @param sizePartForChr size of the part of chr
#' @param dirPd dir of prod
#' @param chunk chunk in h5
#' @param level level in h5
#' @param dir2n dir of 2n
#' @param dirZ dir of zscore
#'
#' @return write out pd and intermediate result
#' @export
#'
#' @examples
sszToPdByChr <- function(chr, fSszInChr,
                         nameSamples, nameMarks,
                         factorTwoNorm, sizePartForChr,
                         dirPd = NULL,
                         chunk, level,
                         dir2n = NULL, dirZ = NULL) {
  library(rhdf5)
  source("R/functions/checkH5.R")
  source("R/functions/pipelineFuns.R")
  source("R/functions/getSpec.R")
  outPd <- !is.null(dirPd)
  out2n <- !is.null(dir2n)
  outZ <- !is.null(dirZ)
  numSamples <- length(nameSamples)
  cat(chr, "\n", sep = "")
  #### check input
  cat("check input h5\n")
  checkMatOneChrH5(fSszInChr, chr, nameSamples, nameMarks)
  numBinsInChr <- readNumBinsInMatOneChrH5(fSszInChr)
  #### create output
  cat("create output h5\n")
  source("R/functions/h5Utils.R")
  if (outPd) {
    fPd <- createDirAndGiveFileName(dirPd, chr, "pd")
    stopifnot(h5createFile(fPd))
    h5createDataset(file = fPd, dataset = chr,
                    dims = c(numBinsInChr, numSamples),
                    storage.mode = "double",
                    chunk = chunk, level = level)
  }
  if (out2n) {
    f2n <- createDirAndGiveFileName(dir2n, chr, "2n")
    stopifnot(h5createFile(f2n))
    h5createDataset(file = f2n, dataset = chr,
                    dims = c(numBinsInChr, numSamples),
                    storage.mode = "double",
                    chunk = chunk, level = level)
  }
  if (outZ) {
    fZ <- createDirAndGiveFileName(dirZ, chr, "z")
    stopifnot(h5createFile(fZ))
    h5createDataset(file = fZ, dataset = chr,
                    dims = c(numBinsInChr, numSamples),
                    storage.mode = "double",
                    chunk = chunk, level = level)
  }
  #### split chr to parts
  cat("split chr to parts\n")
  partsOfChr <- splitIdx(numBinsInChr, sizePartForChr)
  numParts <- length(partsOfChr)
  for (iPart in seq_len(numParts)) {
    cat("Part ", iPart, "\n", sep = "")
    idxBinsInPart <- partsOfChr[[iPart]]
    numBinsInPart <- length(idxBinsInPart)
    #### get ssz in part +++++++++++++++++++++++++++++++++++++++++
    cat("get ssz in part ", iPart, "\n", sep = "")
    matSszInPart <- h5read(file = fSszInChr, name = chr,
                           index = list(idxBinsInPart, NULL))
    #### normalize ssz +++++++++++++++++++++++++++++++++++++++++++
    cat("normalize ssz\n")
    matTwoNormSszInPart <-
      matSszInPart * rep(factorTwoNorm, each = numBinsInPart)
    if (out2n) {
      cat("write h5 of 2n in part", iPart, "\n", sep = "")
      h5write(obj = matTwoNormSszInPart, file = f2n, name = chr,
              index = list(idxBinsInPart, NULL))
    }
    #### get zscore by sample ++++++++++++++++++++++++++++++++++++
    cat("get zscore by sample\n")
    matZscoreInPart <- getZscore(matTwoNormSszInPart)
    if (outZ) {
      cat("write h5 of zscore in part", iPart, "\n", sep = "")
      h5write(obj = matZscoreInPart, file = fZ, name = chr,
              index = list(idxBinsInPart, NULL))
    }
    #### multiply zscore by normalized ssz +++++++++++++++++++++++
    cat("multiply zscore by normalized ssz\n")
    matPdInPart <- matZscoreInPart * matTwoNormSszInPart
    if (outPd) {
      #### write h5 in part
      cat("write h5 of pd in part", iPart, "\n", sep = "")
      h5write(obj = matPdInPart, file = fPd, name = chr,
              index = list(idxBinsInPart, NULL))
    }
  }
  if (out2n) {
    cat("write h5 metadata of 2n\n")
    writeMatH5Metadata(f2n, nameSamples, nameMarks)
  }
  if (outZ) {
    cat("write h5 metadata of zscore\n")
    writeMatH5Metadata(fZ, nameSamples, nameMarks)
  }
  if (outPd) {
    #### write h5 metadata
    cat("write h5 metadata of pd\n")
    writeMatH5Metadata(fPd, nameSamples, nameMarks)
  }
  invisible(NULL)
}
