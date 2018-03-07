#' Get split indices
#'
#' @param sizeTotal total length to be split
#' @param sizePart length of part
#'
#' @return a list of indices of each part
#' @export
#'
#' @examples
splitIdx <- function(sizeTotal, sizePart) {
  parts <- split(seq_len(sizeTotal),
                 findInterval(seq_len(sizeTotal),
                              seq(1, sizeTotal, by = sizePart)))
  return(parts)
}

#' Correct input tensor by mark
#'
#' @param tsrInput tensor of input
#' @param thresInput threshold of input
#' @param maxInputOfMarks max input of marks
#' @param nameSamples names of smaples
#' @param nameMarks names of marks
#' @param numBinsInPart number of bins in part
#'
#' @return list of tsrCtdInput. numSigBinsOfSampleMark and numSigBinsOfMark
#' @export
#'
#' @examples
correctTsrInput <- function(tsrInput, thresInput, maxInputOfMarks,
                            numBinsInPart, nameSamples, nameMarks) {
  numSamples <- length(nameSamples)
  numMarks <- length(nameMarks)
  tsrCtdInput <- array(0, dim = dim(tsrInput))
  numSigBinsOfSampleMark <- array(0L, dim = c(numSamples, numMarks))
  numSigBinsOfMark <- integer(numMarks)
  for (iMark in seq_len(numMarks)) {
    mark <- nameMarks[iMark]
    cat("", mark)
    matInputInMark <- tsrInput[, , iMark]
    matBinaryInMark <- (matInputInMark >= thresInput) * 1L
    numSigBinsOfSampleMark[, iMark] <- colSums(matBinaryInMark)
    isSigAtBinInMark <- (rowSums(matBinaryInMark) > 0) * 1L
    numSigBinsOfMark[iMark] <- sum(isSigAtBinInMark)
    matKeptAtBinSampleInMark <- rep(isSigAtBinInMark, times = numSamples)
    dim(matKeptAtBinSampleInMark) <- c(numBinsInPart, numSamples)
    matCtdInputInMark <- matInputInMark * matKeptAtBinSampleInMark
    if (!any(is.na(maxInputOfMarks))) {
      matCtdInputInMark <- pmin(matCtdInputInMark, maxInputOfMarks[iMark])
    }
    tsrCtdInput[, , iMark] <- matCtdInputInMark
    # .Internal(inspect(tsrCtdInput)) # no change
  }
  cat("\n")
  return(list(tsrCtdInput = tsrCtdInput,
              numSigBinsOfSampleMark = numSigBinsOfSampleMark,
              numSigBinsOfMark = numSigBinsOfMark))
}

#' Get zscore by mark
#'
#' @param tsrScore tensor of score
#' @param nameMarks names of marks
#'
#' @return tensor of zscore
#' @export
#'
#' @examples
getTsrZscore <- function(tsrScore, nameMarks) {
  numMarks <- length(nameMarks)
  source("R/functions/getSpec.R")
  tsrZscore <- array(0, dim = dim(tsrScore))
  for (iMark in seq_len(numMarks)) {
    mark <- nameMarks[iMark]
    cat("", mark)
    matScoreInMark <- tsrScore[, , iMark]
    matZscoreInMark <- getZscore(matScoreInMark)
    tsrZscore[, , iMark] <- matZscoreInMark
  }
  cat("\n")
  return(tsrZscore)
}

#' Get sum squared zscore by sample
#'
#' @param tsrZscore tensor of zscore
#' @param nameSamples names of samples
#'
#' @return list of matSsz and sumSquaredSsz
#' @export
#'
#' @examples
getMatSsz <- function(tsrZscore, nameSamples) {
  numSamples <- length(nameSamples)
  matSsz <- array(0, dim = dim(tsrZscore)[1:2])
  for (iSample in seq_len(numSamples)) {
    sample <- nameSamples[iSample]
    cat("", sample)
    matZscoreInSample <- tsrZscore[, iSample, ]
    sszInSample <- rowSums(matZscoreInSample ^ 2)
    matSsz[, iSample] <- sszInSample
  }
  sumSquaredSsz <- colSums(matSsz ^ 2)
  cat("\n")
  return(list(matSsz = matSsz, sumSquaredSsz = sumSquaredSsz))
}

#' Binarize input tensor by mark
#'
#' @param tsrInput tensor of input
#' @param thresInput threshold of input
#' @param nameSamples names of smaples
#' @param nameMarks names of marks
#' @param numBinsInPart number of bins in part
#'
#' @return list of tsrBinaryInput. numSigBinsOfSampleMark and numSigBinsOfMark
#' @export
#'
#' @examples
binarizeTsrInput <- function(tsrInput, thresInput,
                            numBinsInPart, nameSamples, nameMarks) {
  numSamples <- length(nameSamples)
  numMarks <- length(nameMarks)
  tsrBinaryInput <- array(0L, dim = dim(tsrInput))
  numSigBinsOfSampleMark <- array(0L, dim = c(numSamples, numMarks))
  numSigBinsOfMark <- integer(numMarks)
  for (iMark in seq_len(numMarks)) {
    mark <- nameMarks[iMark]
    cat("", mark)
    matInputInMark <- tsrInput[, , iMark]
    matBinaryInMark <- (matInputInMark >= thresInput) * 1L
    numSigBinsOfSampleMark[, iMark] <- colSums(matBinaryInMark)
    isSigAtBinInMark <- (rowSums(matBinaryInMark) > 0) * 1L
    numSigBinsOfMark[iMark] <- sum(isSigAtBinInMark)
    tsrBinaryInput[, , iMark] <- matBinaryInMark
    # .Internal(inspect(tsrBinaryInput)) # no change
  }
  cat("\n")
  return(list(tsrBinaryInput = tsrBinaryInput,
              numSigBinsOfSampleMark = numSigBinsOfSampleMark,
              numSigBinsOfMark = numSigBinsOfMark))
}

#' Get hamming by mark
#'
#' @param tsrScore tensor of score
#' @param nameMarks names of marks
#'
#' @return tensor of hamming
#' @export
#'
#' @examples
getTsrHamming <- function(tsrScore, nameMarks) {
  numMarks <- length(nameMarks)
  source("R/functions/getSpec.R")
  tsrHamming <- array(0L, dim = dim(tsrScore))
  for (iMark in seq_len(numMarks)) {
    mark <- nameMarks[iMark]
    cat("", mark)
    matScoreInMark <- tsrScore[, , iMark]
    matHammingInMark <- getHamming(matScoreInMark)
    tsrHamming[, , iMark] <- matHammingInMark
  }
  cat("\n")
  return(tsrHamming)
}

#' Get sum absolute hamming by sample
#'
#' @param tsrHamming tensor of hamming
#' @param nameSamples names of samples
#'
#' @return list of matSah and sumSquaredSah
#' @export
#'
#' @examples
getMatSah <- function(tsrHamming, nameSamples) {
  numSamples <- length(nameSamples)
  matSah <- array(0, dim = dim(tsrHamming)[1:2])
  for (iSample in seq_len(numSamples)) {
    sample <- nameSamples[iSample]
    cat("", sample)
    matHammingInSample <- tsrHamming[, iSample, ]
    sahInSample <- rowSums(abs(matHammingInSample))
    matSah[, iSample] <- sahInSample
  }
  sumSquaredSah <- colSums(matSah ^ 2)
  cat("\n")
  return(list(matSah = matSah, sumSquaredSah = sumSquaredSah))
}

