#' Get fold change of two gr
#'
#' @param gr1 gr1
#' @param gr2 gr2
#' @param lenAllChrs length of all chromosome
#'
#' @return a scalar
#' @export
#'
#' @examples
fcGrGr <- function(gr1, gr2, lenAllChrs) {
  library(GenomicRanges)
  gr1 <- reduce(gr1)
  gr2 <- reduce(gr2)
  ov <- intersect(gr1, gr2)
  ovLen <- sum(as.numeric(width(ov)))
  gr1Len <- sum(as.numeric(width(gr1)))
  gr2Len <- sum(as.numeric(width(gr2)))
  fc <- ovLen * lenAllChrs / (gr1Len * gr2Len)
  ## the below is right but would cause asymmetry of result of gr1 and gr2
  # fc <- ovLen * lenAllChrs / gr1Len / gr2Len
  stopifnot(all(fc >= 0))
  return(fc)
}


#' Get fold change between a gr and each gr in grl
#'
#' @param gr a gr
#' @param grl a grl
#' @param lenAllChrs length of all chromosome
#'
#' @return a numeric vector
#' @export
#'
#' @examples
fcGrGrl <- function(gr, grl, lenAllChrs) {
  library(tidyverse)
  library(GenomicRanges)
  stopifnot(is(gr, "GRanges"))
  stopifnot(is(grl, "GRangesList"))
  if (is.null(names(grl)) | "" %in% names(grl))
    names(grl) <- seq_along(grl)
  ## duplication of names would cause big error on the result
  if (any(duplicated(names(grl)))) {
    warning("Names of grl are not unique. Here apply make.unique to them.")
    names(grl) <- make.unique(names(grl), sep = "_")
  }
  ## reduce is a must!
  gr <- reduce(gr)
  grl <- reduce(grl)
  gr1 <- unlist(grl)
  gr2 <- gr
  ov <- findOverlaps(gr1, gr2)
  qh <- queryHits(ov)
  sh <- subjectHits(ov)
  inters <- pintersect(gr1[qh], gr2[sh])
  wd <- as.numeric(width(inters))
  df <- data.frame(name1 = names(gr1)[qh],
                   ovLen = wd)
  un1 <- unique(names(gr1))
  len1 <- sum(as(width(grl), "NumericList"))
  len2 <- sum(as.numeric(width(gr)))
  dfCom <- df %>%
    mutate(name1 = factor(name1, un1)) %>%
    group_by(name1) %>%
    summarise(ovLen = sum(ovLen)) %>%
    ungroup() %>%
    complete(name1, fill = list(ovLen = 0)) %>%
    mutate(len1 = len1[name1],
           len2 = len2,
           fc = lenAllChrs * ovLen / (len1 * len2))
  fc <- setNames(dfCom$fc, dfCom$name1)
  stopifnot(all(fc >= 0))
  return(fc)
}


#' Get fold change between each gr in grl and a gr
#'
#' @param grl a grl
#' @param gr a gr
#' @param lenAllChrs length of all chromosome
#'
#' @return a numeric vector
#' @export
#'
#' @examples
fcGrlGr <- function(grl, gr, lenAllChrs) {
  library(tidyverse)
  library(GenomicRanges)
  stopifnot(is(gr, "GRanges"))
  stopifnot(is(grl, "GRangesList"))
  if (is.null(names(grl)) | "" %in% names(grl))
    names(grl) <- seq_along(grl)
  if (any(duplicated(names(grl)))) {
    warning("Names of grl are not unique. Here apply make.unique to them.")
    names(grl) <- make.unique(names(grl), sep = "_")
  }
  gr <- reduce(gr)
  grl <- reduce(grl)
  gr1 <- unlist(grl)
  gr2 <- gr
  ov <- findOverlaps(gr1, gr2)
  qh <- queryHits(ov)
  sh <- subjectHits(ov)
  inters <- pintersect(gr1[qh], gr2[sh])
  wd <- as.numeric(width(inters))
  df <- data.frame(name1 = names(gr1)[qh],
                   ovLen = wd)
  un1 <- unique(names(gr1))
  len1 <- sum(as(width(grl), "NumericList"))
  len2 <- sum(as.numeric(width(gr)))
  dfCom <- df %>%
    mutate(name1 = factor(name1, un1)) %>%
    group_by(name1) %>%
    summarise(ovLen = sum(ovLen)) %>%
    ungroup() %>%
    complete(name1, fill = list(ovLen = 0)) %>%
    mutate(len1 = len1[name1],
           len2 = len2,
           fc = lenAllChrs * ovLen / (len1 * len2))
  fc <- setNames(dfCom$fc, dfCom$name1)
  stopifnot(all(fc >= 0))
  return(fc)
}


#' Get fold change between each gr in grl1 and each gr in grl2
#'
#' @param grl1 grl1
#' @param grl2 grl2
#' @param lenAllChrs length of all chromosome
#'
#' @return a numeric matrix, rows for grl1, cols for grl2
#' @export
#'
#' @examples
fcGrlGrl <- function(grl1, grl2, lenAllChrs) {
  library(tidyverse)
  library(GenomicRanges)
  stopifnot(is(grl1, "GRangesList"))
  stopifnot(is(grl2, "GRangesList"))
  if (is.null(names(grl1)) | "" %in% names(grl1))
    names(grl1) <- seq_along(grl1)
  if (is.null(names(grl2)) | "" %in% names(grl2))
    names(grl2) <- seq_along(grl2)
  if (any(duplicated(names(grl1)))) {
    warning("Names of grl1 are not unique. Here apply make.unique to them.")
    names(grl1) <- make.unique(names(grl1), sep = "_")
  }
  if (any(duplicated(names(grl2)))) {
    warning("Names of grl2 are not unique. Here apply make.unique to them.")
    names(grl2) <- make.unique(names(grl2), sep = "_")
  }
  grl1 <- reduce(grl1)
  grl2 <- reduce(grl2)
  gr1 <- unlist(grl1)
  gr2 <- unlist(grl2)
  ov <- findOverlaps(gr1, gr2)
  qh <- queryHits(ov)
  sh <- subjectHits(ov)
  inters <- pintersect(gr1[qh], gr2[sh])
  wd <- as.numeric(width(inters))
  df <- data.frame(name1 = names(gr1)[qh],
                   name2 = names(gr2)[sh],
                   ovLen = wd)
  un1 <- unique(names(gr1))
  un2 <- unique(names(gr2))
  len1 <- sum(as(width(grl1), "NumericList"))
  len2 <- sum(as(width(grl2), "NumericList"))
  dfCom <- df %>%
    mutate(name1 = factor(name1, un1),
           name2 = factor(name2, un2)) %>%
    group_by(name1, name2) %>%
    summarise(ovLen = sum(ovLen)) %>%
    ungroup() %>%
    complete(name1, name2, fill = list(ovLen = 0)) %>%
    mutate(len1 = len1[name1],
           len2 = len2[name2],
           fc = lenAllChrs * ovLen / (len1 * len2))
  ar <- reshape2::acast(dfCom, name1 ~ name2, value.var = "fc")
  stopifnot(all(ar >= 0))
  return(ar)
}


#' Get fold change between each gr in grl1 and each gr in grl2. Use it when
#' fcGrlGrl cannot be used for grl1 or grl2 is too large
#'
#' @param grl1 grl1
#' @param grl2 grl2
#' @param lenAllChrs length of all chromosome
#' @param para whether parallel or not. If TRUE, para on grl2.
#'
#' @return a numeric matrix, rows for grl1, cols for grl2
#' @export
#'
#' @examples
fcGrlGrl_2 <- function(grl1, grl2, lenAllChrs, para = TRUE,
                       paraBy = c("second", "first")) {
  if (para) {
    library(parallel)
    paraBy <- match.arg(paraBy)
    # print(paraBy)
    source("R/functions/fc.R")
    if (paraBy == "second") {
      numCores <- min(detectCores(), length(grl2))
      cl <- makeCluster(numCores)
      on.exit(stopCluster(cl))
      fc <- clusterApplyLB(cl, grl2, fcGrGrl, grl1, lenAllChrs)
      names(fc) <- names(grl2)
      fc <- do.call(cbind, fc)
    } else {
      numCores <- min(detectCores(), length(grl1))
      cl <- makeCluster(numCores)
      on.exit(stopCluster(cl))
      fc <- clusterApplyLB(cl, grl1, fcGrGrl, grl2, lenAllChrs)
      names(fc) <- names(grl1)
      fc <- do.call(rbind, fc)
    }
  } else {
    fc <- vapply(grl2, fcGrGrl, double(length(grl1)), grl1, lenAllChrs)
  }
  stopifnot(all(fc >= 0))
  return(fc)
}





# fcGrGrl_2 <- function(gr, grl, lenAllChrs) {
#   fc <- vapply(grl, fcGrGr, double(1), gr, lenAllChrs)
#   return(fc)
# }

# fcGrlGr_2 <- function(grl, gr, lenAllChrs) {
#   fc <- vapply(grl, fcGrGr, double(1), gr, lenAllChrs)
#   return(fc)
# }


# library(tidyverse)
# library(GenomicRanges)
# source("R/functions/getMetadata.R")
# getMetadata("roadmap")
# sizeChrs <- seqlengths(seqInfo)[nameChrs]
# lenAllChrs <- sum(as.numeric(sizeChrs))

# gg <- fcGrGr(csreGRL[[1]], geneticGRL[[1]], lenAllChrs)
# gl <- fcGrGrl(csreGRL[[1]], csreGRL[1:3], lenAllChrs)
# gl2 <- fcGrlGr(csreGRL, csreGRL[[1]], lenAllChrs)
