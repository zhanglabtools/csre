#' get overlapped elements for each csre
#'
#' @param csreGR GR of csre
#' @param elementGR GR of element
#' @param sort sort overlapped element of each csre
#' @param uniq uniq overlapped element of each csre
#' @param resType return chr or gr with mcols
#'
#' @return overlapped elements for each csre
#' @export
#'
#' @examples
getOverlappedElements <- function(csreGR, elementGR,
                                  sort = TRUE, uniq = TRUE,
                                  resType = c("chr", "gr")) {
  stopifnot(is.logical(uniq), is.logical(uniq))
  library(GenomicRanges)
  resType <- match.arg(resType)
  ov <- findOverlaps(csreGR, elementGR)
  qh <- queryHits(ov)
  sh <- subjectHits(ov)
  nameSh <- elementGR$tx_name[sh]
  f <- factor(qh, levels = seq_len(queryLength(ov)))
  nameByCsre <- unname(split(nameSh, f))
  if (sort)
    nameByCsre <- lapply(nameByCsre, sort)
  if (uniq)
    nameByCsre <- lapply(nameByCsre, unique)
  nameByCsre <- unlist(lapply(nameByCsre, paste, collapse = ","))
  if (resType == "chr")
    return(nameByCsre)
  mcols(csreGR)$elements <- nameByCsre
  return(csreGR)
}
# getOverlappedLengthWrong <- function(csreGR, elementGR, resType = c("int", "gr")) {
#   resType <- match.arg(resType)
#   if (!isDisjoint(elementGR))
#     elementGR <- reduce(elementGR, ignore.strand = TRUE)
#   ovp <- findOverlapPairs(csreGR, elementGR)
#   inters <- pintersect(ovp)
#   wd <- width(inters)
#   ## if some range is duplicated n times, then the answer would be n times right
#   ## length, beacause then belong to the same level of f, so I use findOverlaps
#   ## which can distinguish them by id
#   f <- as.factor(first(ovp))
#   wd <- split(wd, f)
#   len <- vapply(wd, sum, integer(1))
#   gr <- GRanges(names(len))
#   mcols(gr)$len <- len
#   mt <- match(csreGR, gr)
#   resLen <- len[mt]
#   resLen[is.na(resLen)] <- 0L
#   resLen <- unname(resLen)
#   if (resType == "int") return(resLen)
#   mcols(csreGR)$len <- as(resLen, "Rle")
#   return(csreGR)
# }

#' get overlapped lengths for each csre
#'
#' @param csreGR GR of csre
#' @param elementGR GR of element
#' @param resType return int or gr with mcols
#'
#' @return overlapped lengths for each csre
#' @export
#'
#' @examples
getOverlappedLength <- function(csreGR, elementGR,
                                resType = c("int", "gr")) {
  library(GenomicRanges)
  resType <- match.arg(resType)
  if (!isDisjoint(elementGR))
    elementGR <- reduce(elementGR, ignore.strand = TRUE)
  ov <- findOverlaps(csreGR, elementGR)
  qh <- queryHits(ov)
  sh <- subjectHits(ov)
  inters <- pintersect(csreGR[qh], elementGR[sh])
  wd <- width(inters)
  f <- factor(qh, levels = seq_len(queryLength(ov)))
  wd <- split(wd, f)
  len <- unname(vapply(wd, sum, integer(1)))
  if (resType == "int") return(len)
  mcols(csreGR)$len <- as(len, "Rle")
  return(csreGR)
}

#' get overlapped lengths of each genetic feature for each csre
#'
#' @param csreGR GR of csre
#' @param geneticGRL GRL of genetic features
#' @param resType return data.frame or matrix
#' @param para whether parallel or not
#'
#' @return data.frame of overlapped lengths of each genetic feature for each
#'   csre
#' @export
#'
#' @examples
getOverlappedGeneticLength <- function(csreGR, geneticGRL,
                                       resType = c("df", "mat"),
                                       para = TRUE) {
  if (para == TRUE) {
    library(parallel)
    f <- function(x) getOverlappedLength(csreGR, x)
    numCores <- min(detectCores(), length(geneticGRL))
    cl <- makeCluster(numCores)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, source("R/functions/getOvEleAndLen.R"))
    #### the below is to avoid copying csreGR to each node multiple times
    ## to avoid copying variables in current env
    environment(f) <- .GlobalEnv
    ## library pkg first before exporting csreGR, to avoid risk of non-S4 funs
    clusterEvalQ(cl, library("GenomicRanges"))
    ## envir should be environment(), or global env will be used
    clusterExport(cl, "csreGR", envir = environment())
    ## use balanced clusterApplyLB
    l <- clusterApplyLB(cl, geneticGRL, f)
    ## assign names to make it like lapply
    names(l) <- names(geneticGRL)
  } else {
    l <- lapply(geneticGRL, function(x) getOverlappedLength(csreGR, x))
  }
  resType <- match.arg(resType)
  if (resType == "df") {
    res <- data.frame(l)
  } else {
    res <- do.call(cbind, l)
  }
  return(res)
}

#' get genetic category of each csre
#'
#' @param csreGR GR of csre
#' @param geneticGRL GRL of genetic features
#' @param precedence if TRUE, the first to last element (feature) in geneticGRL
#'   has highest to lowest priority and a csre only belongs to no more than one
#'   feature; if FALSE, a csre belongs to each feature it overlaps
#' @param resTypeOfPre type of result when precedence is TRUE
#'
#' @return a 0-1 matrix where 1 indicate the csre in that row belongs to feature
#'   in that column, or a chraracter, factor or factor-Rle when precedence is
#'   TRUE
#' @export
#'
#' @examples
getCategory <- function(csreGR, geneticGRL, precedence = TRUE,
                        resTypeOfPre = c("mat", "chr", "fct", "Rle", "all")) {
  stopifnot(is.logical(precedence))
  resTypeOfPre <- match.arg(resTypeOfPre)
  csreGeneticRes <- getOverlappedGeneticLength(csreGR, geneticGRL,
                                               resType = "mat")
  cn <- sub("Bp", "", colnames(csreGeneticRes))
  category <- (csreGeneticRes > 0L) * 1L
  colnames(category) <- cn
  if (!precedence)
    return(category)
  library(matrixStats)
  cumCategory <- rowCummaxs(category)
  categoryWithPrecedence <- cbind(category[, 1], rowDiffs(cumCategory))
  colnames(categoryWithPrecedence) <- paste0(cn, "Pre")
  if (resTypeOfPre == "mat")
    return(categoryWithPrecedence)
  idx <- which(categoryWithPrecedence == 1L, arr.ind = TRUE)
  idx <- idx[order(idx[, "row"]), "col"]
  chrCate <- cn[idx]
  if (resTypeOfPre == "chr")
    return(chrCate)
  fctCate <- factor(chrCate, levels = cn)
  if (resTypeOfPre == "fct")
    return(fctCate)
  RleCate <- Rle(fctCate)
  if (resTypeOfPre == "Rle")
    return(RleCate)
  if (resTypeOfPre == "all")
    return(list(matOv = csreGeneticRes,
                matCate = categoryWithPrecedence,
                chrCate = chrCate,
                fctCate = fctCate,
                RleCate = RleCate))
}

#' find intervals in each gr of grl that overlaps another gr
#'
#' countOverlaps can only do the reverse thing so I write this function
#' @param grl a grl
#' @param gr another gr
#'
#' @return a integer vector indicate how many intervals in each gr of grl that
#'   overlaps another gr
#' @export
#'
#' @examples
revCountOverlapsGRLGR <- function(grl, gr) {
  library(GenomicRanges)
  ov <- findOverlaps(unlist(grl), gr)
  qh <- queryHits(ov)
  res <- unname(table(findInterval(unique(qh), start(grl@partitioning))))
  return(res)
}
