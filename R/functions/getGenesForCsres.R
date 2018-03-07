#' Get overlapped genes
#'
#' @param grl grl of CSREs
#' @param NMWithPr gr of NMs with promoters
#'
#' @return list of overlapped genes for each gr of CSREs
#' @export
#'
#' @examples
getOvGenes <- function(grl, NMWithPr) {
  library(GenomicRanges)
  genesByOv <- lapply(grl, function(gr) {
    sg <- subsetByOverlaps(NMWithPr, gr)
    ge <- sort(unique(sg$tx_name))
    return(ge) } )
  return(genesByOv)
}

#' Get nearest genes
#'
#' @param grl grl of CSREs
#' @param NM gr of NMs
#'
#' @return list of overlapped genes for each gr of CSREs
#' @export
#'
#' @examples
getNearestGenesOld <- function(grl, NM) {
  library(GenomicRanges)
  tss <- promoters(NM, 0, 1)
  genesByNearest <- lapply(grl, function(gr) {
    sg <- tss[nearest(gr, tss)]
    ge <- sort(unique(sg$tx_name))
    return(ge) } )
  return(genesByNearest)
}

#' Get nearest genes
#'
#' @param grl grl of CSREs
#' @param NM gr of NMs
#'
#' @return list of overlapped genes for each gr of CSREs
#' @export
#'
#' @examples
getNearestGenes <- function(grl, NM) {
  library(GenomicRanges)
  tss <- promoters(NM, 0, 1)
  genesByNearest <- lapply(grl, function(gr) {
    sg <- tss[unique(subjectHits(nearest(gr, tss,
                                         ignore.strand = TRUE,
                                         select = "all")))]
    ge <- sort(unique(sg$tx_name))
    return(ge) } )
  return(genesByNearest)
}

#' Pairwise overlap of list og genes
#'
#' @param listGenes list of genes
#'
#' @return matrix of pairwise overlap of genes
#' @export
#'
#' @examples
pOvGenes <- function(listGenes) {
  num <- length(listGenes)
  matOv <- array(0, dim = rep(num, 2))
  for (i in 1:num) {
    for (j in i:num) {
      matOv[i, j] <- matOv[j, i] <-
        length(intersect(listGenes[[i]],
                         listGenes[[j]]))
    }
  }
  return(matOv)
}
