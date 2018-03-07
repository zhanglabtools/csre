#' Check metadata of h5 file of the track tensor of one chromosome
#'
#' @param file h5 file
#' @param chr the chromosome of that h5 file
#' @param nameSamples names of samples to be checked
#' @param nameMarks names of marks to be checked
#' @param numBinsInChr number of bins in the chr to be checked
#' @param dims names of dims to be checked
#'
#' @return if the check were passed then return invisible NULL, otherwise stop
#' @export
#'
#' @examples
checkTsrOneChrH5 <- function(file, chr = NULL,
                             nameSamples = NULL,
                             nameMarks = NULL,
                             numBinsInChr = NULL,
                             dims = c("bin", "sample", "mark")) {
  library(rhdf5)
  info <- h5ls(file)
  isChr <- grepl("^chr", info$name)
  if (!is.null(chr))
    stopifnot(identical(chr, info$name[isChr]))
  nameMetadata <- info$name[!isChr]
  metadata <- list()
  for (meta in nameMetadata)
    metadata[[meta]] <- h5read(file = file, name = meta)
  ## I want dims to be a vector of length 3, so the below should not run
  # dims <- match.arg(dims)
  stopifnot(identical(dims, as.vector(metadata[["dims"]])))
  lenDims <- as.double(unlist(strsplit(info$dim[isChr], " x ")))
  if (!is.null(numBinsInChr))
    stopifnot(numBinsInChr == lenDims[1])
  if (!is.null(nameSamples))
    stopifnot(identical(nameSamples, as.vector(metadata[["nameSamples"]])),
              length(nameSamples) == lenDims[2])
  if (!is.null(nameMarks))
    stopifnot(identical(nameMarks, as.vector(metadata[["nameMarks"]])),
              length(nameMarks) == lenDims[3])
  invisible(NULL)
}

#' Get number of bins of h5 file of the track tensor of one chromosome after
#' checking metadata of it
#'
#' @param file the h5 file to get number of bins in its chromosome
#'
#' @return number of bins in the chromosome of that h5 file
#' @export
#'
#' @examples
readNumBinsInTsrOneChrH5 <- function(file) {
  library(rhdf5)
  info <- h5ls(file)
  isChr <- grepl("^chr", info$name)
  numBinsInChr <- as.double(unlist(strsplit(info$dim[isChr], " x ")))[1]
  return(numBinsInChr)
}

#' Check metadata of h5 file of the track matrix of one chromosome
#'
#' @param file h5 file
#' @param chr the chromosome of that h5 file
#' @param nameSamples names of samples to be checked
#' @param summedMarks names of marks to be checked
#' @param numBinsInChr number of bins in the chr to be checked
#' @param dims names of dims to be checked
#'
#' @return if the check were passed then return invisible NULL, otherwise stop
#' @export
#'
#' @examples
checkMatOneChrH5 <- function(file, chr = NULL,
                             nameSamples = NULL,
                             summedMarks = NULL,
                             numBinsInChr = NULL,
                             dims = c("bin", "sample")) {
  library(rhdf5)
  info <- h5ls(file)
  isChr <- grepl("^chr", info$name)
  if (!is.null(chr))
    stopifnot(identical(chr, info$name[isChr]))
  nameMetadata <- info$name[!isChr]
  metadata <- list()
  for (meta in nameMetadata)
    metadata[[meta]] <- h5read(file = file, name = meta)
  ## I want dims to be a vector of length 3, so the below should not run
  # dims <- match.arg(dims)
  stopifnot(identical(dims, as.vector(metadata[["dims"]])))
  lenDims <- as.double(unlist(strsplit(info$dim[isChr], " x ")))
  if (!is.null(numBinsInChr))
    stopifnot(numBinsInChr == lenDims[1])
  if (!is.null(nameSamples))
    stopifnot(identical(nameSamples, as.vector(metadata[["nameSamples"]])),
              length(nameSamples) == lenDims[2])
  if (!is.null(summedMarks))
    stopifnot(identical(summedMarks, as.vector(metadata[["summedMarks"]])))
  invisible(NULL)
}

#' Get number of bins of h5 file of the track matrix of one chromosome after
#' checking metadata of it
#'
#' @param file the h5 file to get number of bins in its chromosome
#'
#' @return number of bins in the chromosome of that h5 file
#' @export
#'
#' @examples
readNumBinsInMatOneChrH5 <- function(file) {
  library(rhdf5)
  info <- h5ls(file)
  isChr <- grepl("^chr", info$name)
  numBinsInChr <- as.double(unlist(strsplit(info$dim[isChr], " x ")))[1]
  return(numBinsInChr)
}

#' Check metadata of h5 file of the track matrix of all chromosomes
#'
#' @param file h5 file
#' @param nameChrs the chromosomes of that h5 file
#' @param nameSamples names of samples to be checked
#' @param summedMarks names of marks to be checked
#' @param numBins number of bins in all chrs to be checked
#' @param dims names of dims to be checked
#'
#' @return if the check were passed then return invisible NULL, otherwise stop
#' @export
#'
#' @examples
checkMatAllChrsH5 <- function(file, nameChrs = NULL,
                              nameSamples = NULL,
                              summedMarks = NULL,
                              numBins = NULL,
                              dims = c("bin", "sample")) {
  library(rhdf5)
  info <- h5ls(file)
  isChr <- grepl("^chr", info$name)
  if (!is.null(nameChrs))
    stopifnot(identical(sort(nameChrs), sort(info$name[isChr])))
  nameMetadata <- info$name[!isChr]
  metadata <- list()
  for (meta in nameMetadata)
    metadata[[meta]] <- h5read(file = file, name = meta)
  stopifnot(identical(dims, as.vector(metadata[["dims"]])))
  lenDims <- as.double(unlist(strsplit(info$dim[isChr], " x ")))
  if (!is.null(nameSamples))
    stopifnot(unique(lenDims[seq(2, length(lenDims), by = 2)]) ==
                length(nameSamples),
              identical(nameSamples, as.vector(metadata[["nameSamples"]])))
  if (!is.null(summedMarks))
    stopifnot(identical(summedMarks, as.vector(metadata[["summedMarks"]])))
  if (!is.null(numBins)) {
    numBinsInFile <- lenDims[seq(1, length(lenDims), by = 2)]
    names(numBinsInFile) <- info$name[isChr]
    stopifnot(identical(sort(numBinsInFile), sort(numBins)))
  }
  invisible(NULL)
}

#' Get number of bins of h5 file of the track matrix of all chromosomes after
#' checking
#'
#' @param file the h5 file to get number of bins in all chromosomes
#' @param nameChrs user supplied chrs to check chrs in H5 and reorder the result
#'
#' @return number of bins in all chromosomes of that h5 file in the order of
#'   nameChrs
#' @export
#'
#' @examples
readNumBinsInMatAllChrsH5 <- function(file, nameChrs) {
  library(rhdf5)
  info <- h5ls(file)
  isChr <- grepl("^chr", info$name)
  lenDims <- as.double(unlist(strsplit(info$dim[isChr], " x ")))
  numBinsInFile <- lenDims[seq(1, length(lenDims), by = 2)]
  names(numBinsInFile) <- info$name[isChr]
  stopifnot(identical(sort(names(numBinsInFile)), sort(nameChrs)))
  numBins <- numBinsInFile[nameChrs]
  return(numBins)
}
