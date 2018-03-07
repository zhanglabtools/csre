#' Batch fisher.test
#'
#' @param overlappedLen length of overlapped region
#' @param csreLen length of csre
#' @param elementLen length of element
#' @param totalLen length of total background
#' @param log.p whether get negtive log p value
#' @param alternative auto or greater or less
#'
#' @return fold change and p value
#' @export
#'
#' @examples
batchFisherTest <- function(overlappedLen, csreLen, elementLen, totalLen,
                            log.p = TRUE,
                            alternative = c("auto", "greater", "less")) {
  alternative <- match.arg(alternative, several.ok = FALSE)
  nonElementLen <- totalLen - elementLen
  fc <- overlappedLen * totalLen / elementLen / csreLen
  if (log.p) {
    if (alternative %in% c("auto", "greater"))
      logpg <- phyper(overlappedLen - 1, elementLen, nonElementLen, csreLen,
                      lower.tail = FALSE, log.p = TRUE) / -log(10)
    if (alternative %in% c("auto", "less"))
      logpl <- phyper(overlappedLen, elementLen, nonElementLen, csreLen,
                      lower.tail = TRUE, log.p = TRUE) / -log(10)
    if (alternative == "auto") {
      logpm <- pmax(logpg, logpl)
      return(list(fc = fc, logp = logpm))
    }
    if (alternative == "greater")
      return(list(fc = fc, logp = logpg))
    if (alternative == "less")
      return(list(fc = fc, logp = logpl))
  } else {
    if (alternative %in% c("auto", "greater"))
      pg <- phyper(overlappedLen - 1, elementLen, nonElementLen, csreLen,
                   lower.tail = FALSE, log.p = FALSE)
    if (alternative %in% c("auto", "less"))
      pl <- phyper(overlappedLen, elementLen, nonElementLen, csreLen,
                   lower.tail = TRUE, log.p = FALSE)
    if (alternative == "auto") {
      pm <- pmin(pg, pl)
      return(list(fc = fc, pvalue = pm))
    }
    if (alternative == "greater")
      return(list(fc = fc, pvalue = pg))
    if (alternative == "less")
      return(list(fc = fc, pvalue = pl))
  }
}
