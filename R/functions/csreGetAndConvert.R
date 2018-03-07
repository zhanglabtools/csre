#' Get GR of CSRE from data.frame
#'
#' @param csre data.frame with chr, start, end, cellType
#' @param seqInfo a Seqinfo object to give the final seqinfo
#'
#' @return a GR of CSRE
#' @export
#'
#' @examples
getCsreGR <- function(csre, seqInfo) {
  library(GenomicRanges)
  csreGR <- with(csre,
                 GRanges(chr,
                         IRanges(start = start + 1L,
                                 end = end),
                         cellType = as(factor(cellType, unique(cellType)),
                                       "Rle")))
  seqinfo(csreGR) <- GenomicRanges::intersect(seqInfo, seqinfo(csreGR))
  # seqinfo(csreGR) <- seqInfo
  return(csreGR)
}

#' Get BinGR of CSRE from data.frame
#'
#' @param csreBin data.frame with chr, startBin, endBin, cellType
#' @param seqInfo a Seqinfo object to be binned to the final seqinfo
#' @param sizeBin how many base pairs in the non-overlapping bins
#' @param keepLastIncompleteBin should the last bin be kept when it is not
#'   complete
#'
#' @return
#' @export
#'
#' @examples
getCsreBinGR <- function(csreBin, seqInfo, sizeBin,
                         keepLastIncompleteBin = TRUE) {
  library(GenomicRanges)
  csreBinGR <- with(csreBin,
                    GRanges(chr,
                            IRanges(start = startBin,
                                    end = endBin),
                            cellType = as(factor(cellType, unique(cellType)),
                                          "Rle")))
  seqinfo(csreBinGR) <- GenomicRanges::intersect(seqInfo, seqinfo(csreBinGR))
  # seqinfo(csreBinGR) <- seqInfo
  genome(csreBinGR) <- paste0(unique(genome(csreBinGR)), "Bin", sizeBin)
  if (keepLastIncompleteBin) {
    seqlengths(csreBinGR) <- (seqlengths(csreBinGR) - 1L) %/% sizeBin + 1L
  } else {
    seqlengths(csreBinGR) <- seqlengths(csreBinGR) %/% sizeBin
    csreBinGR <- trim(csreBinGR)
  }
  return(csreBinGR)
}

#' Convert GRanges from bin based to base pair based
#'
#' @param gr a GRanges object
#' @param sizeBin how many base pairs in the non-overlapping bins
#' @param seqInfo a Seqinfo object to give the final seqinfo but only seqlevels
#'   in gr are keeped
#'
#' @return the converted GRanges object
#' @export
#'
#' @examples
binToBp <- function(gr, sizeBin, seqInfo) {
  library(GenomicRanges)
  seqLevels <- seqlevels(gr)
  seqinfo(gr) <- seqInfo
  seqlevels(gr) <- seqLevels
  end(gr) <- end(gr) * sizeBin
  ## the below is identical to the above, but more consistent with treat of
  ## start
  # end(gr) <- (end(gr) - 1L) * sizeBin + sizeBin
  start(gr) <- (start(gr) - 1L) * sizeBin + 1L
  gr <- trim(gr)
  return(gr)
}

#' Convert GRanges from base pair based to bin based
#'
#' @param gr a GRanges object
#' @param sizeBin how many base pairs in the non-overlapping bins
#' @param keepLastIncompleteBin should the last bin be kept when it is not
#'   complete
#'
#' @return the converted GRanges object
#' @export
#'
#' @examples
bpToBin <- function(gr, sizeBin, keepLastIncompleteBin = TRUE) {
  library(GenomicRanges)
  start(gr) <- (start(gr) - 1L) %/% sizeBin + 1L
  end(gr) <- (end(gr) - 1L) %/% sizeBin + 1L
  genome(gr) <- paste0(genome(gr), "Bin", sizeBin)
  if (keepLastIncompleteBin) {
    seqlengths(gr) <- (seqlengths(gr) - 1L) %/% sizeBin + 1L
  } else {
    seqlengths(gr) <- seqlengths(gr) %/% sizeBin
    gr <- trim(gr)
  }
  return(gr)
}

#' Split gr to grl by cellType in mcols
#'
#' @param gr a gr with cellType in mcols
#'
#' @return a grl split by cellType
#' @export
#'
#' @examples
splitByCellType <- function(gr) {
  split(gr, with(mcols(gr), factor(decode(cellType),
                                   levels = unique(cellType))))
}

#' Convert gr to bed
#'
#' @param gr the gr
#'
#' @return a data.frame of bed
#' @export
#'
#' @examples
grToBed <- function(gr) {
  return(data.frame(chr = seqnames(gr), start = start(gr) - 1L, end = end(gr)))
}

#' Get csre GR from file
#'
#' @param file file path
#' @param type one of df, gr, grl, binGr, binGrl
#' @param nameChrs names of chromosome
#' @param nameSamples names of samples
#' @param seqInfo an Seqinfo
#' @param sizeBin size of bin
#' @param keepLastIncompleteBin keepLastIncompleteBin in bpToBin()
#'
#' @return a gr of csre
#' @export
#'
#' @examples
getCsreGRFromFile <- function(file,
                              type = c("df", "gr", "grl", "binGr", "binGrl"),
                              nameChrs, nameSamples,
                              seqInfo, sizeBin, keepLastIncompleteBin) {
  type <- match.arg(type)
  input <- data.table::fread(file, data.table = FALSE)
  library(magrittr)
  csre <- input %>%
    dplyr::transmute(cellType = factor(Cell_Type, nameSamples),
                     chr = factor(Chromosome, nameChrs),
                     start = Start_ZeroBased, end = End_ZeroBased,
                     len = Length)
  if (type == "df") return(csre)
  csreGR <- getCsreGR(csre, seqInfo)
  if (type == "gr") return(csreGR)
  if (type == "grl") {
    csreGRL <- splitByCellType(csreGR)
    return(csreGRL)
  }
  csreBinGR <- bpToBin(csreGR, sizeBin, keepLastIncompleteBin)
  if (type == "binGr") return(csreBinGR)
  csreBinGRL <- splitByCellType(csreBinGR)
  return(csreBinGRL)
}

#' Get csre GR from file of old Chen
#'
#' @param file file path
#' @param type one of df, gr, grl, binGr, binGrl
#' @param nameChrs names of chromosome
#' @param nameSamples names of samples
#' @param seqInfo an Seqinfo
#' @param sizeBin size of bin
#' @param keepLastIncompleteBin keepLastIncompleteBin in bpToBin()
#'
#' @return a gr of csre
#' @export
#'
#' @examples
getCsreGRFromFileOldChen <-
  function(file = "data/CSREsOfDCMA.txt",
           type = c("df", "gr", "grl", "binGr", "binGrl"),
           nameChrs, nameSamples,
           seqInfo, sizeBin, keepLastIncompleteBin) {
  type <- match.arg(type)
  input <- data.table::fread(file, data.table = FALSE)
  library(magrittr)
  csre <- input %>%
    dplyr::transmute(cellType = factor(Cell_Type, nameSamples),
                     chr = factor(Chromosome, nameChrs),
                     start = Start_Position, end = End_Position + 1L,
                     len = Length)
  if (type == "df") return(csre)
  csreGR <- getCsreGR(csre, seqInfo)
  if (type == "gr") return(csreGR)
  if (type == "grl") {
    csreGRL <- splitByCellType(csreGR)
    return(csreGRL)
  }
  csreBinGR <- bpToBin(csreGR, sizeBin, keepLastIncompleteBin)
  if (type == "binGr") return(csreBinGR)
  csreBinGRL <- splitByCellType(csreBinGR)
  return(csreBinGRL)
}

