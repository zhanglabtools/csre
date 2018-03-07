getSpecByEntropy <- function(track, byrow = TRUE)
{
  if (!byrow)
  {
    track <- t(track)
  }
  track <- pmax(track, array(0, dim = dim(track)))
  dims <- dim(track)[2]
  track <- track + (rowSums(track) == 0)
  track <- track / rowSums(track)
  track[track == 0] <- 1
  track <- -track * log(track, base = dims)
  track[is.nan(track)] <- 0
  entropy <- rowSums(track)
  spec <- 1 - entropy
  return(spec)
}

getSpecByCV <- function(track, byrow = TRUE, biasedsd = FALSE, symmetry = FALSE)
{
  if (!byrow)
  {
    track <- t(track)
  }
  track <- pmax(track, array(0, dim = dim(track)))
  track <- track + (rowSums(track) == 0)
  rowmeans <- rowMeans(track)
  dims <- dim(track)[2]
  sdtrack <- sqrt(rowSums((track - rowmeans) ^ 2) / (dims - 1))
  cv <- sdtrack / rowmeans
  if (biasedsd)
  {
    cv <- cv * sqrt((dims - 1) / dims)
  }
  if (symmetry)
  {
    rowidx <- seq_len(dim(track)[1])
    rowmaxs <- track[cbind(rowidx, max.col(m = track))]
    rowmins <- track[cbind(rowidx, max.col(m = -track))]
    symmeans <- rowmaxs + rowmins - rowmeans
    cv <- cv * sqrt(rowmeans / symmeans)
  }
  return(cv)
}

getSpecTrackByPairwiseMinus <- function(
  track, pairbyrow = TRUE, above = 0, absinpair = FALSE, scale = TRUE,
  divisor = c("mean", "midrange", "geomean", "minmean", "sd", "mad",
              "ave"))
{
  divisor <- match.arg(divisor)
  if (!pairbyrow)
  {
    track <- t(track)
  }
  if (!is.na(above))
  {
    track <- pmax(track, array(above, dim = dim(track)))
  }
  track <- track + (rowSums(track) == 0)
  rowmeans <- rowMeans(track)
  spectrack <- array(0, dim = dim(track))
  colnames(spectrack) <- colnames(track)
  dims <- dim(track)[2]
  if (absinpair)
  {
    for (isam in seq_len(dims))
    {
      spectrack[, isam] <- rowSums(abs(track[, isam] - track)) / dims
    }
    if (scale)
    {
      if (divisor == "ave")
      {
        ave <- rowMeans(spectrack)
        ave <- ave + (ave == 0)
        spectrack <- spectrack / ave
      }
    }
  } else
  {
    spectrack <- track - rowmeans
    if (scale)
    {
      if (divisor == "mean")
      {
        spectrack <- spectrack / rowmeans
      } else if (divisor %in% c("midrange", "geomean", "minmean"))
      {
        rowidx <- seq_len(dim(track)[1])
        rowmaxs <- track[cbind(rowidx, max.col(m = track))]
        rowmins <- track[cbind(rowidx, max.col(m = -track))]
        symmeans <- rowmaxs + rowmins - rowmeans
        if (divisor == "midrange")
          spectrack <- spectrack / ((rowmeans + symmeans) / 2)
        if (divisor == "geomean")
          spectrack <- spectrack / sqrt(rowmeans * symmeans)
        if (divisor == "minmean")
          spectrack <- spectrack / pmin(rowmeans, symmeans)
      } else
      {
        if (divisor == "sd")
        {
          sdtrack <- sqrt(rowSums(spectrack ^ 2) / (dims - 1))
          sdzero <- sdtrack == 0
          spectrack <- spectrack / (sdtrack + sdzero)
        }
        if (divisor == "mad")
        {
          mad <- rowMeans(abs(spectrack))
          madzero <- mad == 0
          spectrack <- spectrack / (mad + madzero)
        }
      }
    }
  }
  return(spectrack)
}

#' transform each row (or col) of a matrix to zscore
#'
#' @param track the matrix, each col is a track, each row is a bin
#' @param pairbyrow tranform row to zscore if TRUE, otherwise col
#' @param above the lower bound of all numbers if the matrix
#'
#' @return the tansformed matrix
#' @export
#'
#' @examples
getZscore <- function(track, pairbyrow = TRUE, above = NA)
{
  if (!pairbyrow)
  {
    track <- t(track)
  }
  if (!is.na(above))
  {
    track <- pmax(track, array(above, dim = dim(track)))
  }
  rowmeans <- rowMeans(track)
  dims <- dim(track)[2]
  trackminusmean <- track - rowmeans
  sdtrack <- sqrt(rowSums(trackminusmean ^ 2) / (dims - 1))
  sdzero <- sdtrack == 0
  track <- trackminusmean / (sdtrack + sdzero)
  return(track)
}

#' transform each row (or col) of a binary matrix to hamming
#'
#' @param track the binary matrix, each col is a track, each row is a bin
#' @param pairbyrow tranform row to hamming if TRUE, otherwise col
#' @param abs whether the difference is signed
#'
#' @return the tansformed matrix
#' @export
#'
#' @examples
getHamming <- function(track, pairbyrow = TRUE, abs = FALSE) {
  if (ncol(track) < 2)
    stop("number of columns should be no less than 2")
  if (!pairbyrow) {
    if (nrow(track) < 2)
      stop("number of rows should be no less than 2")
    track <- t(track)
  }
  nc <- ncol(track)
  if (abs) {
    # ham <- array(0L, dim = dim(track), dimnames = dimnames(track))
    # for (i in 1:(nc - 1)) {
    #   for (j in (i + 1):nc) {
    #     tmp <- track[, i] != track[, j]
    #     ham[, i] <- ham[, i] + tmp
    #     ham[, j] <- ham[, j] + tmp
    #   }
    # }
    ## the above is quadratic to ncol VS the below is linear ncol, when ncol is
    ## big, the below is significantly faster than the above (I have tested
    ## that)
    rs <- rowSums(track)
    storage.mode(rs) <- "integer"
    ## for ones, abs hamming is number of zeros in its row, for zeros, abs
    ## hamming is number of ones in its row
    ham <- track * (nc - rs) + (1L - track) * rs
  } else {
    # ham <- array(0L, dim = dim(track), dimnames = dimnames(track))
    # for (i in 1:(nc - 1)) {
    #   for (j in (i + 1):nc) {
    #     tmp <- track[, i] - track[, j]
    #     ham[, i] <- ham[, i] + tmp
    #     ham[, j] <- ham[, j] - tmp
    #   }
    # }
    ## the above is quadratic to ncol VS the below is linear ncol, when ncol is
    ## big, the below is significantly faster than the above (I have tested
    ## that)
    rs <- rowSums(track)
    storage.mode(rs) <- "integer"
    ## for ones, hamming is number of zeros in its row, for zeros, hamming is
    ## negative number of ones in its row
    ham <- track * (nc - rs) + (track - 1L) * rs
  }
  return(ham)
}
