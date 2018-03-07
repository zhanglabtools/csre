

#' Convert list of genes (from REFSEQ to ENTREZID)
#'
#' @param lg list of genes
#' @param orgDb org.db
#' @param from an id type
#' @param to an id type
#'
#' @return converted list of genes
#' @export
#'
#' @examples
convertListGeneName <- function(lg, orgDb = org.Hs.eg.db,
                                from = "REFSEQ", to = "ENTREZID") {
  library(org.Hs.eg.db)
  ag <- unique(unlist(lg))
  tb <- select(orgDb, ag, to, from)
  entrezByCellType <- lapply(lg, function(x)
    unique(
      na.omit(
        # subset(tb, REFSEQ %in% x, ENTREZID, drop = TRUE)
        tb[[to]][tb[[from]] %in% x]
        )
      )
    )
  return(entrezByCellType)
}

#' Get table about enrichment of GO of list of genes
#'
#' @param lg list of genes
#' @param cellType names of list of genes
#' @param orgDb org.db
#' @param GODb GO.db
#'
#' @return a data.frame
#' @export
#'
#' @examples
GOEnrichListGeneSet <- function(lg, cellType = names(lg),
                                orgDb = org.Hs.egGO2ALLEGS,
                                GODb = GO.db) {
  # orgDb <- org.Hs.egGO2ALLEGS
  # orgDb <- org.Hs.egGO2EG
  if (missing(orgDb)) library(org.Hs.eg.db)
  if (missing(GODb)) library(GO.db)
  source("R/functions/batchFisherTest.R")
  bgGenesByGO <- as.list(orgDb)
  # saveRDS(bgGenesByGO, "data/bgGenesByGO.rds")
  # bgGenesByGO <- readRDS("data/bgGenesByGO.rds")
  lenUniNoNa <- function(x) length(unique(na.omit(x)))
  sizeBgGenesByGO <- vapply(bgGenesByGO, lenUniNoNa, integer(1))
  bgGenes <- mappedLkeys(orgDb)
  sizeBgGenes <- length(bgGenes)
  GOId <- mappedRkeys(orgDb)
  GOontology <- Ontology(GOId)
  tb <- select(GO.db, GOId, c("ONTOLOGY" ,"TERM"), "GOID")
  ## all pairs are considered
  GOEnrichPerGeneSet <- function(geneSet, curCellType) {
    subOrgDb <- subset(orgDb, Lkeys = geneSet)
    geneSetByGO <- as.list(subOrgDb)
    sizeGeneSetByGO <- vapply(geneSetByGO, lenUniNoNa, integer(1))
    sizeGeneSet <- sum(geneSet %in% bgGenes)
    # sizeGeneSet <- length(mappedLkeys(subOrgDb))
    resFs <- data.frame(batchFisherTest(sizeGeneSetByGO,
                                        sizeBgGenesByGO,
                                        sizeGeneSet,
                                        sizeBgGenes,
                                        alternative = "greater"))
    resFs <- cbind(cellType = curCellType,
                   resFs,
                   ontology = GOontology,
                   sizeOv = sizeGeneSetByGO,
                   sizeTerm = sizeBgGenesByGO,
                   sizeQueryGenes = sizeGeneSet,
                   sizeBgGenes = sizeBgGenes,
                   stringsAsFactors = FALSE)
    # tmp <- select(GODb, rownames(resFs),
    #               c("ONTOLOGY" ,"TERM"), "GOID")
    tmp <- subset(tb, GOID == GOId)
    resFs <- merge(tmp,
                   resFs,
                   by.x = "GOID", by.y = 0)
    resFs <- resFs[order(resFs$logp, decreasing = TRUE), ]
    return(resFs)
  }
  if (is.null(cellType) | "" %in% cellType)
    cellType <- seq_along(lg)
  res <- mapply(GOEnrichPerGeneSet, lg, cellType, SIMPLIFY = FALSE)
  res <- do.call(rbind, res)
  res$cellType <- factor(res$cellType, cellType)
  attr(res, "row.names") <- .set_row_names(nrow(res))
  library(dplyr)
  res <- res %>%
    mutate(logqByPairs = -log(p.adjust(10 ^ -logp, method = "fdr",
                                       n = n()),
                              base = 10)) %>%
    group_by(cellType) %>%
    mutate(logqByCellType = -log(p.adjust(10 ^ -logp, method = "fdr",
                                          n = n()),
                                 base = 10),
           fcMultiplyNlp = fc * logqByCellType) %>%
    as.data.frame()
  return(res)
}

#' Get table about enrichment of GWAS of csreGR
#'
#' @param csreGR a gr
#' @param gwasGR a gr of GWAS
#' @param gwasIdName id of GWAS
#'
#' @return a data.frame
#' @export
#'
#' @examples
gwasEnrich <- function(csreGR, gwasGR, gwasIdName) {
  library(dplyr)
  library(GenomicRanges)
  source("R/functions/batchFisherTest.R")
  gwasGRL <- split(gwasGR, mcols(gwasGR)[, gwasIdName], drop = TRUE)
  ## not an elegent way to get factor pairs
  study <- unique(mcols(unlist(gwasGRL)[, gwasIdName]))
  numStudy <- length(gwasGRL)
  gwasGRLNoOv <- unique(gwasGRL)
  # gwasGRLNoOv <- endoapply(gwasGRL, unique)
  csreGRL <- split(csreGR, csreGR$cellType)
  numSamples <- length(csreGRL)
  ## unlike GOEnrichListGeneSet(), only overlapped pairs are used to calculate
  ## p-value
  ov <- findOverlaps(csreGRL, gwasGRLNoOv)
  widthCsreGRL <- sum(width(csreGRL))
  numSNPsInCellType <- countOverlaps(csreGRL, gwasGR)
  numSNPsInGWAS <- sum(!duplicated(gwasGR))
  widthGwasGRLNoOv <- sum(width(gwasGRLNoOv))
  qh <- queryHits(ov)
  sh <- subjectHits(ov)
  pOv <- intersect(csreGRL[qh], gwasGRLNoOv[sh])
  mcols(ov) <- DataFrame(
    lengthQuery = widthCsreGRL[qh],
    numSNPsInQuery = numSNPsInCellType[qh],
    lengthSubject = widthGwasGRLNoOv[sh],
    lengthOv = sum(width(pOv)),
    ovGR = pOv)
  lenAllChrs <- sum(as.numeric(seqlengths(csreGR)))
  enrichSNP <- batchFisherTest(mcols(ov)$lengthOv,
                               mcols(ov)$numSNPsInQuery,
                               mcols(ov)$lengthSubject,
                               numSNPsInGWAS,
                               alternative = "greater")
  dfEnrichSNP <- cbind(cellType = factor(names(csreGRL)[qh]),
                       as.data.frame(study[sh, , drop = FALSE],
                                     optional = TRUE),
                       as.data.frame(enrichSNP),
                       as.data.frame(mcols(ov), optional = TRUE))
  res <- dfEnrichSNP %>%
    mutate(logqByPairs = -log(p.adjust(10 ^ -logp, method = "fdr",
                                       n = numStudy * numSamples),
                              base = 10)) %>%
    dplyr::rename(lenCsreOfCellType = lengthQuery,
                  numUniqueSNPsOfStudy = lengthSubject,
                  numHits = lengthOv) %>%
    mutate(recall = numHits / numUniqueSNPsOfStudy) %>%
    group_by(cellType) %>%
    mutate(logqCellType = -log(p.adjust(10 ^ -logp, method = "fdr",
                                        n = numStudy),
                               base = 10)) %>%
    arrange(cellType, desc(logqByPairs), desc(fc)) %>%
    dplyr::select(cellType, logqByPairs, logqCellType, logp, fc, recall,
                  numHits, numUniqueSNPsOfStudy, lenCsreOfCellType,
                  everything()) %>%
    select_(.dots = c(gwasIdName, names(.))) %>%
    as.data.frame()
  return(res)
}
