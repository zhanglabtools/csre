#' Get ssz from input, correction, zscore and ssz by chr
#'
#' @param chr chromosome
#' @param fInputInChr input file in chr
#' @param nameSamples name of samples
#' @param nameMarks name of marks
#' @param thresInput threshold of input
#' @param sizePartForChr size of the part of chr
#' @param dirSsz dir of ssz
#' @param dirStats dir of stats
#' @param chunk chunk in h5
#' @param level level in h5
#' @param nameInput name of input, e.g. "nlp"
#'
#' @return write out ssz and stats
#' @export
#'
#' @examples
inputToSszByChrOld <- function(chr, fInputInChr,
                               nameSamples, nameMarks,
                               thresInput, sizePartForChr,
                               dirSsz, dirStats,
                               chunk, level, nameInput) {
  library(rhdf5)
  source("R/functions/checkH5.R")
  source("R/functions/pipelineFuns.R")
  numSamples <- length(nameSamples)
  numMarks <- length(nameMarks)
  cat(chr, "\n", sep = "")
  #### check input
  cat("check input h5\n")
  checkTsrOneChrH5(fInputInChr, chr, nameSamples, nameMarks)
  numBinsInChr <- readNumBinsInTsrOneChrH5(fInputInChr)
  #### create output
  cat("create output h5\n")
  if (!file.exists(dirSsz)) dir.create(dirSsz, recursive = TRUE)
  if (!file.exists(dirStats)) dir.create(dirStats, recursive = TRUE)
  fSsz <- file.path(dirSsz, paste0("ssz_", chr, ".h5"))
  stopifnot(h5createFile(fSsz))
  h5createDataset(file = fSsz, dataset = chr,
                  dims = c(numBinsInChr, numSamples),
                  storage.mode = "double",
                  chunk = chunk, level = level)
  #### create stats
  cat("create stats\n")
  numSigBinsOfSampleMarkInChr <- array(0L, dim = c(numSamples, numMarks),
                                       dimnames = list(nameSamples, nameMarks))
  numSigBinsOfMarkInChr <- setNames(integer(numMarks), nameMarks)
  sumSquaredSszOfSampleInChr <- setNames(double(numSamples), nameSamples)
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
    cat("correction by mark\n")
    resCorrectTsrInput <-
      correctTsrInput(tsrInputInPart, thresInput,
                      numBinsInPart, nameSamples, nameMarks)
    tsrCtdInputInPart <- resCorrectTsrInput$tsrCtdInput
    numSigBinsOfSampleMarkInChr <-
      numSigBinsOfSampleMarkInChr + resCorrectTsrInput$numSigBinsOfSampleMark
    numSigBinsOfMarkInChr <-
      numSigBinsOfMarkInChr + resCorrectTsrInput$numSigBinsOfMark
    rm(resCorrectTsrInput, tsrInputInPart)
    #### get zscore by mark +++++++++++++++++++++++++++++++++++++++++
    cat("get zscore by mark\n")
    tsrZscoreInPart <- getTsrZscore(tsrCtdInputInPart, nameMarks)
    rm(tsrCtdInputInPart)
    #### get sum squared zscore by sample +++++++++++++++++++++++++++
    cat("get sum squared zscore by sample\n")
    resGetMatSsz <- getMatSsz(tsrZscoreInPart, nameSamples)
    matSszInPart <- resGetMatSsz$matSsz
    sumSquaredSszOfSampleInChr <-
      sumSquaredSszOfSampleInChr + resGetMatSsz$sumSquaredSsz
    rm(resGetMatSsz, tsrZscoreInPart)
    #### write h5 in part
    cat("write h5 in part", iPart, "\n", sep = "")
    h5write(obj = matSszInPart, file = fSsz, name = chr,
            index = list(idxBinsInPart, NULL))
    cat("\n")
  }
  #### write h5 metadata
  cat("write h5 metadata\n")
  h5write(obj = c("bin", "sample"), file = fSsz, name = "dims", level = 0)
  h5write(obj = nameSamples, file = fSsz, name = "nameSamples", level = 0)
  h5write(obj = nameMarks, file = fSsz, name = "summedMarks", level = 0)
  #### write stats
  cat("write stats\n")
  saveRDS(list(chr = chr,
               numSigBinsOfSampleMarkInChr = numSigBinsOfSampleMarkInChr,
               numSigBinsOfMarkInChr = numSigBinsOfMarkInChr,
               sumSquaredSszOfSampleInChr = sumSquaredSszOfSampleInChr),
          file.path(dirStats, paste0(nameInput, "ToSsz_", chr, ".rds")))
  invisible(NULL)
}
