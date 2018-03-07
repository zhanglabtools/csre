#' Get ssz from input, correction, zscore and ssz by chr
#'
#' @param chr chromosome
#' @param fInputInChr input file in chr
#' @param nameSamples name of samples
#' @param nameMarks name of marks
#' @param thresInput threshold of input
#' @param maxInputOfMarks max input of marks
#' @param sizePartForChr size of the part of chr
#' @param dirSsz dir of ssz
#' @param dirStats dir of stats
#' @param chunkMat chunk of mat h5
#' @param levelMat level of mat h5
#' @param nameInput name of input, e.g. "nlp"
#' @param dirCtd dir of corrected score
#' @param dirZscore dir of zscore
#' @param chunkTsr chunk of tensor
#' @param levelTsr level of  tensor
#'
#' @return write out ssz and stats and intermediate result
#' @export
#'
#' @examples
inputToSszByChr <- function(chr, fInputInChr,
                            nameSamples, nameMarks,
                            thresInput, maxInputOfMarks,
                            sizePartForChr,
                            dirSsz = NULL,
                            dirStats = NULL,
                            chunkMat, levelMat,
                            nameInput,
                            dirCtd = NULL,
                            dirZscore = NULL,
                            chunkTsr, levelTsr) {
  library(rhdf5)
  source("R/functions/checkH5.R")
  source("R/functions/pipelineFuns.R")
  numSamples <- length(nameSamples)
  numMarks <- length(nameMarks)
  outSsz <- !is.null(dirSsz)
  outStats <- !is.null(dirStats)
  outCtd <- !is.null(dirCtd)
  outZscore <- !is.null(dirZscore)
  cat(chr, "\n", sep = "")
  #### check input
  cat("check input h5\n")
  checkTsrOneChrH5(fInputInChr, chr, nameSamples, nameMarks)
  numBinsInChr <- readNumBinsInTsrOneChrH5(fInputInChr)
  #### create output
  cat("create output h5\n")
  source("R/functions/h5Utils.R")
  if (outStats) {
    if (!file.exists(dirStats)) dir.create(dirStats, recursive = TRUE)
  }
  if (outSsz) {
    fSsz <- createDirAndGiveFileName(dirSsz, chr, "ssz")
    stopifnot(h5createFile(fSsz))
    h5createDataset(file = fSsz, dataset = chr,
                    dims = c(numBinsInChr, numSamples),
                    storage.mode = "double",
                    chunk = chunkMat, level = levelMat)
  }
  if (outCtd) {
    fCtd <- createDirAndGiveFileName(dirCtd, chr, "c")
    stopifnot(h5createFile(fCtd))
    h5createDataset(file = fCtd, dataset = chr,
                    dims = c(numBinsInChr, numSamples, numMarks),
                    storage.mode = "double",
                    chunk = chunkTsr, level = levelTsr)
  }
  if (outZscore) {
    fZscore <- createDirAndGiveFileName(dirZscore, chr, "z")
    stopifnot(h5createFile(fZscore))
    h5createDataset(file = fZscore, dataset = chr,
                    dims = c(numBinsInChr, numSamples, numMarks),
                    storage.mode = "double",
                    chunk = chunkTsr, level = levelTsr)
  }
  #### create stats
  if (outStats) {
    cat("create stats\n")
    numSigBinsOfSampleMarkInChr <-
      array(0L, dim = c(numSamples, numMarks),
            dimnames = list(nameSamples, nameMarks))
    numSigBinsOfMarkInChr <- setNames(integer(numMarks), nameMarks)
    sumSquaredSszOfSampleInChr <- setNames(double(numSamples), nameSamples)
  }
  #### split chr to parts
  cat("split chr to parts\n")
  partsOfChr <- splitIdx(numBinsInChr, sizePartForChr)
  numParts <- length(partsOfChr)
  #### main
  for (iPart in seq_len(numParts)) {
    cat("Part ", iPart, "\n", sep = "")
    idxBinsInPart <- partsOfChr[[iPart]]
    numBinsInPart <- length(idxBinsInPart)
    #### get input in part +++++++++++++++++++++++++++++++++++++++++
    cat("get input in part ", iPart, "\n", sep = "")
    tsrInputInPart <- h5read(file = fInputInChr, name = chr,
                             index = list(idxBinsInPart, NULL, NULL))
    #### correction by mark +++++++++++++++++++++++++++++++++++++++++
    cat("correct by mark\n")
    resCorrectTsrInput <-
      correctTsrInput(tsrInputInPart, thresInput, maxInputOfMarks,
                      numBinsInPart, nameSamples, nameMarks)
    tsrCtdInputInPart <- resCorrectTsrInput$tsrCtdInput
    if (outStats) {
      numSigBinsOfSampleMarkInChr <-
        numSigBinsOfSampleMarkInChr + resCorrectTsrInput$numSigBinsOfSampleMark
      numSigBinsOfMarkInChr <-
        numSigBinsOfMarkInChr + resCorrectTsrInput$numSigBinsOfMark
    }
    rm(resCorrectTsrInput, tsrInputInPart)
    if (outCtd) {
      cat("write h5 of corrected score in part", iPart, "\n", sep = "")
      h5write(obj = tsrCtdInputInPart, file = fCtd, name = chr,
              index = list(idxBinsInPart, NULL, NULL))
    }
    #### get zscore by mark +++++++++++++++++++++++++++++++++++++++++
    cat("get zscore by mark\n")
    tsrZscoreInPart <- getTsrZscore(tsrCtdInputInPart, nameMarks)
    rm(tsrCtdInputInPart)
    if (outZscore) {
      cat("write h5 of zscore in part", iPart, "\n", sep = "")
      h5write(obj = tsrZscoreInPart, file = fZscore, name = chr,
              index = list(idxBinsInPart, NULL, NULL))
    }
    #### get sum squared zscore by sample +++++++++++++++++++++++++++
    cat("get sum squared zscore by sample\n")
    resGetMatSsz <- getMatSsz(tsrZscoreInPart, nameSamples)
    matSszInPart <- resGetMatSsz$matSsz
    if (outStats) {
      sumSquaredSszOfSampleInChr <-
        sumSquaredSszOfSampleInChr + resGetMatSsz$sumSquaredSsz
    }
    rm(resGetMatSsz, tsrZscoreInPart)
    if (outSsz) {
      cat("write h5 of ssz in part", iPart, "\n", sep = "")
      h5write(obj = matSszInPart, file = fSsz, name = chr,
              index = list(idxBinsInPart, NULL))
    }
  }
  if (outCtd) {
    cat("write h5 metadata of corrected score\n")
    writeTsrH5Metadata(fCtd, nameSamples, nameMarks)
  }
  if (outZscore) {
    cat("write h5 metadata of zscore\n")
    writeTsrH5Metadata(fZscore, nameSamples, nameMarks)
  }
  if (outSsz) {
    cat("write h5 metadata of ssz\n")
    writeMatH5Metadata(fSsz, nameSamples, nameMarks)
  }
  if (outStats) {
    fStats <- file.path(dirStats, paste0(nameInput, "ToSsz_", chr, ".rds"))
    if (!file.exists(fStats)) {
      #### write stats
      cat("write stats\n")
      saveRDS(list(chr = chr,
                   numSigBinsOfSampleMarkInChr = numSigBinsOfSampleMarkInChr,
                   numSigBinsOfMarkInChr = numSigBinsOfMarkInChr,
                   sumSquaredSszOfSampleInChr = sumSquaredSszOfSampleInChr),
              fStats)
    }
  }
  invisible(NULL)
}
