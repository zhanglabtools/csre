#' Write metadata to matrix h5 file
#'
#' @param file h5 file
#' @param nameSamples names of samples
#' @param summedMarks names of marks
#' @param level level of h5 about metadata
#'
#' @return invisible()
#' @export
#'
#' @examples
writeMatH5Metadata <- function(file, nameSamples, summedMarks, level = 0) {
  h5write(obj = c("bin", "sample"), file = file, name = "dims", level = level)
  h5write(obj = nameSamples, file = file, name = "nameSamples", level = level)
  h5write(obj = nameMarks, file = file, name = "summedMarks", level = level)
  invisible()
}

#' Write metadata to tensor h5 file
#'
#' @param file h5 file
#' @param nameSamples names of samples
#' @param nameMarks names of marks
#' @param level level of h5 about metadata
#'
#' @return invisible()
#' @export
#'
#' @examples
writeTsrH5Metadata <- function(file, nameSamples, nameMarks, level = 0) {
  h5write(obj = c("bin", "sample", "mark"), file = file, name = "dims",
          level = level)
  h5write(obj = nameSamples, file = file, name = "nameSamples",
          level = level)
  h5write(obj = nameMarks, file = file, name = "nameMarks",
          level = level)
  invisible()
}

#' Create dir if it does not exist and give the file name of the chr to write
#'
#' @param theDir the dir to create
#' @param chr the chromosome
#' @param prefix the file name would be prefix_chr.h5
#'
#' @return the file name corresponds to the chromosome
#' @export
#'
#' @examples
createDirAndGiveFileName <- function(theDir, chr, prefix) {
  if (!file.exists(theDir)) dir.create(theDir, recursive = TRUE)
  f <- file.path(theDir, paste0(prefix, "_", chr, ".h5"))
  return(f)
}

#' Get h5 file names of nameChrs in dir with prefix
#'
#' @param dir dir of h5 files
#' @param prefix prefix of h5 files
#' @param nameChrs name of chromosome correspond to h5 files
#'
#' @return file names of h5 files of chr in dir with names of chr
#' @export
#'
#' @examples
getH5FilesWithChrName <- function(dir, prefix, nameChrs) {
  basename <- paste0(prefix, "_", nameChrs, ".h5")
  stopifnot(basename %in% list.files(dir))
  f <- file.path(dir, basename)
  fWithChrName <- setNames(f, nameChrs)
  return(fWithChrName)
}
