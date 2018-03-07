## getRandomIntervals

## We want to sample each valid permutation (set of intervals) with same
## probability. To sample a valid permutation is hard, so we use the following
## methods. We sample each interval uniformly and independently along the chr.
## In this way each valid set of intervals will be get in same probability.
## Although this method can give invalid set of intervals (that is to say the
## sum of probabilities of each valid permutation is less than 1), we can sample
## again and again until a valid set is sampled. We can prove that this process
## can get each valid permutation with same probability. The key is that in each
## sample turn all the valid permutation have the same probability to be
## sampled.

# bed <- csre[csre[, 1] == "H1", c(2:4,1)]
# bed <- bed[order(bed[, 1]), ]

#' replace each interval randomly
#'
#' @param bed a data.frame contains the intervals to be replaced, only first
#'   columns used
#' @param sizechrs named vector with size of chrs
#' @param method to what extent intervals on the same chr can be overlapped
#' @param numtries max number of tries to sample intervals on each chr
#' @param seed set seed for sample.int
#' @param verbose show messages
#'
#' @return a data.frame with one to one replaced intervals
#' @export
#'
#' @examples
getRandomIntervals <- function(bed, sizechrs,
                               method = c("arbitrary",
                                          "nonoverlapped",
                                          "nonconsectutive"),
                               numtries = 1e6,
                               seed = NA,
                               verbose = 0) {
  method <- match.arg(method)
  if (!is.na(seed)) set.seed(seed)
  ## listbed <- split(bed, factor(bed[[1]], levels = unique(bed[[1]]))) split
  ## bed as the order that each chr first come out, not the order of levels of
  ## the bed[[1]] or as.factor(bed[[1]]), because the levels of bed[[1]] may be
  ## created by read.table(), which is sorted(unique(bed[[1]])), not the
  ## unique(bed[[1]]) as expected by users
  stopifnot(identical(setdiff(bed[[1]], names(sizeChrs)), character(0)))
  if (is.character(bed[[1]])) {
    isFactorChr <- FALSE
    listbed <- split(bed, bed[[1]])
  } else if (is.factor(bed[[1]])) {
    isFactorChr <- TRUE
    listbed <- split(bed, droplevels(bed[[1]]))
  } else {
    stop("chr is not chatacter or factor")
  }

  listran <- list()
  # for (curchr in unique(bed[, 1])) {
  for (curchr in names(listbed)) {
    curbed <- listbed[[curchr]]
    if (verbose > 0) cat(curchr, "\n")
    num <- nrow(curbed)
    if (num > 0) {
      cursizechr <- sizechrs[curchr]
      curlens <- curbed[, 3] - curbed[, 2]
      k = 1
      done = FALSE
      while (!done & (k <= numtries)) {
        if (method == "arbitrary") {
          ranstart <- sample.int(cursizechr, num, replace = TRUE) - 1L
        } else {
          ## Even though replace = FALSE makes the intervals not entirely
          ## indenpent like what replace = TRUE does, it can also gives all
          ## valid permutations the same probability. And replace = FALSE can
          ## have higher probability to sample a valid set of nonoverlapped
          ## intervals than replace = TRUE.
          ranstart <- sample.int(cursizechr, num, replace = FALSE) - 1L
        }
        ranintervals <- cbind(ranstart, ranstart + curlens)
        isinbound <- all(ranintervals[, 2] <= cursizechr)
        if (method == "arbitrary") {
          done <- isinbound
        } else {
          ## Use ordered to test whether the random intervals are overlapped or
          ## consectutive. It is compatibible for one row ranintervals because
          ## all(logical(0)) is TRUE
          ordered <- ranintervals[order(ranintervals[, 1]), , drop = FALSE]
          if (method == "nonconsectutive") {
            isnonconsectutive <- all(ordered[, 2][-num] - ordered[, 1][-1] < 0)
            done <- isinbound & isnonconsectutive
          } else {
            isnonoverlapped <- all(ordered[, 2][-num] - ordered[, 1][-1] <= 0)
            done <- isinbound & isnonoverlapped
          }
        }
        if (verbose > 1) cat(k, "")
        k = k + 1
      }
    } else {
      ranintervals <- array(0, dim = c(0, 2))
      done <- TRUE
    }
    stopifnot(done)
    if (nrow(ranintervals) == 0) {
      listran[[curchr]] <- data.frame(chr = character(0),
                                      start = double(0),
                                      end = double(0))
    } else {
      listran[[curchr]] <- data.frame(curchr, ranintervals,
                                      stringsAsFactors = FALSE)
    }
    # stringsAsFactors = FALSE to let unsplit do right
    names(listran[[curchr]]) <- c("chr", "start", "end")
    rownames(listran[[curchr]]) <- rownames(listbed[[curchr]])
    if (verbose > 0) cat("\n")
  }
  # res <- Reduce(rbind, listran, data.frame())
  ## or we can ues res <- data.frame(data.table::rbindlist(listran))
  ## keep the original format of bed
  # oldlevel <- levels(bed[[1]])
  # if (is.null(oldlevel)) {
  #   res[[1]] <- as.character(res[[1]])
  # } else {
  #   res[[1]] <- factor(res[[1]], levels = oldlevel)
  # }
  if (!isFactorChr) {
    res <- unsplit(listran, bed[[1]])
  } else {
    res <- unsplit(listran, droplevels(bed[[1]]))
  }

  oldlevel <- levels(bed[[1]])
  if (!is.null(oldlevel))
    res[[1]] <- factor(res[[1]], levels = oldlevel)
  stopifnot(identical(bed[[1]], res[[1]]))
  # attributes(res)$row.names <- attributes(bed)$row.names
  ## the above still can not let the two row.names be the same
  ## the below can really let the two row.names be the same
  attr(res, "row.names") <- .row_names_info(bed, type = 0)
  ## the above and below have same effects
  # attributes(res)$row.names <- .row_names_info(bed, type = 0)
  stopifnot(identical(.row_names_info(bed, type = 0),
                      .row_names_info(res, type = 0)))
  return(res)
}

#' replace each interval randomly
#'
#' @param gr a gr contains the intervals to be replaced, only seqnames, start
#'   and end are used
#' @param sizeChrs named vector with size of chrs
#' @param withMcols whether keep mcols of gr
#' @param ... other args of getRandomIntervals
#'
#' @return a gr with one to one replaced intervals
#' @export
#'
#' @examples
getRandomGR <- function(gr, sizeChrs, withMcols = TRUE, ...) {
  stopifnot(is.logical(withMcols))
  source("R/functions/csreGetAndConvert.R")
  bed <- grToBed(gr)
  randomBed <- getRandomIntervals(bed, sizeChrs, ...)
  randomGR <- with(randomBed,
                   GRanges(chr,
                           IRanges(start = start + 1L,
                                   end = end),
                           seqinfo = seqinfo(gr)))
  if (withMcols)
    mcols(randomGR) <- mcols(gr)
  stopifnot(identical(seqnames(randomGR), seqnames(gr)),
            identical(width(randomGR), width(gr)),
            identical(seqinfo(randomGR), seqinfo(gr)))
  return(randomGR)
}

#' replace each interval randomly in supplied region, intervals may overlap
#'
#' @param gr a gr contains the intervals to be replaced
#' @param sizeChrs named vector with size of chrs
#' @param regionGR supplied region that random intervals should be within
#' @param withMcols whether keep mcols of gr
#' @param verbose verbose
#' @param ... other args of getRandomGR
#'
#' @return a gr
#' @export
#'
#' @examples
getRandomGRInRegion <- function(gr, sizeChrs, regionGR,
                                withMcols = TRUE, verbose = 0, ...) {
  # set.seed(1)
  tmp <- gr
  grNow <- gr
  if (!withMcols) {
    mcols(tmp) <- NULL
    mcols(grNow) <- NULL
  }
  idxNow <- seq_along(gr)
  complementGR <- subset(gaps(regionGR), strand == "*")
  i <- 0
  while (1) {
    i <- i + 1
    if (verbose > 0)
      cat("(", i, "," ,length(idxNow), ") ")
    rdmGR <- getRandomGR(grNow, sizeChrs, ...)
    subIdxReSample <- unique(queryHits(findOverlaps(rdmGR, complementGR)))
    subIdxOK <- setdiff(seq_along(idxNow), subIdxReSample)
    tmp[idxNow[subIdxOK]] <- rdmGR[subIdxOK]
    idxNow <- idxNow[subIdxReSample]
    if (length(idxNow) == 0) break
    grNow <- grNow[subIdxReSample]
  }
  return(tmp)
}

#' rep gr in grl
#'
#' This version can succeed in a larger probability but is 2 times slower
#' compared to repGRInGRL_2
#'
#' @param csreGRL grl
#' @param sizeChrs size of chromosome
#' @param numReps number of reps
#' @param verbose verbose
#' @param ... other args of getRandomGR
#'
#' @return a list of cell type's grl of rep random gr
#' @export
#'
#' @examples
repGRInGRL <- function(csreGRL, sizeChrs, numReps = 1e2, verbose = 0, ...) {
  library(GenomicRanges)
  source("R/enrichment/getRandomIntervals.R")
  listCtGRLOfRepRdmGR <- list()
  for (i in seq_along(csreGRL)) {
    if (verbose > 0) cat(i, "")
    gr <- csreGRL[[i]]
    listCtGRLOfRepRdmGR[[i]] <-
      GRangesList(replicate(numReps, getRandomGR(gr, sizeChrs, ...)))
  }
  if (verbose > 0) cat("\n")
  names(listCtGRLOfRepRdmGR) <- names(csreGRL)
  numCsres <- lengths(csreGRL)
  return(listCtGRLOfRepRdmGR)
}

#' rep gr in grl, different intervals may overlap
#'
#' @param csreGRL grl
#' @param sizeChrs size of chromosome
#' @param numReps number of reps
#' @param verbose verbose
#'
#' @return a list of cell type's grl of rep random gr
#' @export
#'
#' @examples
repGRInGRL_2 <- function(csreGRL, sizeChrs, numReps = 1e2, verbose = 0) {
  library(GenomicRanges)
  source("R/enrichment/getRandomIntervals.R")
  source("R/functions/csreGetAndConvert.R")
  grCt <- unlist(csreGRL)
  if (verbose > 0) cat("Rep grl now\n")
  listRepRdmGRL <-
    replicate(numReps, splitByCellType(getRandomGR(grCt, sizeChrs)))
  listCtGRLOfRepRdmGR <- list()
  for (name in names(csreGRL)) {
    if (verbose > 0) cat(name, "")
    listCtGRLOfRepRdmGR[[name]] <-
      GRangesList(lapply(listRepRdmGRL, `[[`, name))
  }
  if (verbose > 0) cat("\n")
  return(listCtGRLOfRepRdmGR)
}
