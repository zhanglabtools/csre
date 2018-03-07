#' smooth input by chr
#'
#' @param chr chromosome
#' @param fInputInChr input file in chr
#' @param nameSamples name of samples
#' @param nameMarks name of marks
#' @param dirSh dir of smooth
#' @param chunk chunk in h5
#' @param level level in h5
#' @param numMaxCores max number of cores given by user
#' @param isBsub called by a R script by bsub or not
#'
#' @return write out sh
#' @export
#'
#' @examples
shByChr <- function(chr, fInputInChr,
                    nameSamples, nameMarks,
                    dirSh,
                    chunk, level,
                    numMaxCores, isBsub) {
  ts <- Sys.time()
  library(rhdf5)
  library(parallel)
  source("R/functions/checkH5.R")
  getMra <- function(x) {
    library(wavelets)
    mra <- mra(x, n.levels = 7, fast = TRUE, method = "modwt")
    res <- mra@S[[3]][, 1]
    rm(mra)
    gc() # rm and gc together can have effect
    ## the above two lines are equal to the below in this context
    # on.exit({rm(mra);gc()})
    ## res can not removed in on.exit(), I don't know why
    return(res)
  }
  numSamples <- length(nameSamples)
  cat(chr, "\n", sep = "")
  #### check input
  cat("check input h5\n")
  checkMatOneChrH5(fInputInChr, chr, nameSamples, nameMarks)
  numBinsInChr <- readNumBinsInMatOneChrH5(fInputInChr)
  #### read data
  cat("read data of ", numBinsInChr, " and ", numSamples, " samples\n",
      sep = "")
  matInputInChr <- h5read(file = fInputInChr, name = chr, index = NULL)
  listInput <- list()
  for (iSample in seq_len(numSamples))
    listInput[[iSample]] <- matInputInChr[, iSample]
  rm(matInputInChr)
  gc()
  #### mra
  if (!isBsub) {
    numAvailClusters <- detectCores()
    numCores <- min(numSamples, numAvailClusters, numMaxCores)
    cat("make cluster of ", numCores, " nodes on localhost\n", sep = "")
    cl <- makeCluster(numCores)
  } else {
    hosts <- strsplit(Sys.getenv("LSB_HOSTS"), split = " ")[[1]]
    usedHosts <- hosts[-1]
    if (length(usedHosts) > numSamples)
      usedHosts <- usedHosts[seq_len(numSamples)]
    numCores <- length(usedHosts)
    cat("make cluster of ", numCores, " nodes on hosts:\n", sep = "")
    cat(usedHosts, "\n", sep = " ")
    cl <- makeCluster(usedHosts)
  }
  ## Set enclosing env to .GlovalEnv, otherwise all objects (very large and
  ## unnecessary in my scene) of its enclosing env will be sent to each cluster
  ## node along with getMra it self by clusterApplyLB and other such functions.
  ## I noticed the phenomenon that each cluster copy all objects in shByChr, but
  ## I didn't know why and how to prevent it, until I read page 25 in this book
  ## (http://detritus.fundacioace.com/pub/books/Oreilly.Parallel.R.Oct.2011.pdf)
  environment(getMra) <- .GlobalEnv
  cat("mra on cluster of ", numCores, " nodes\n", sep = "")
  results <- clusterApplyLB(cl, listInput, getMra)
  # Sys.sleep(5)
  cat("stop cluster of ", numCores, " nodes\n", sep = "")
  stopCluster(cl)
  rm(listInput)
  gc()
  mra <- do.call("cbind", results)
  rm(results)
  gc()
  #### write results
  cat("write results of ", numBinsInChr, " and ", numSamples, " samples\n",
      sep = "")
  if (!file.exists(dirSh)) dir.create(dirSh, recursive = TRUE)
  fSh <- file.path(dirSh, paste0("sh_", chr, ".h5"))
  stopifnot(h5createFile(fSh))
  h5createDataset(file = fSh, dataset = chr,
                  dims = c(numBinsInChr, numSamples),
                  storage.mode = "double",
                  chunk = chunk, level = level)
  h5write(obj = mra, file = fSh, name = chr)
  rm(mra)
  gc()
  #### write h5 metadata
  cat("write h5 metadata\n")
  h5write(obj = c("bin", "sample"), file = fSh, name = "dims", level = 0)
  h5write(obj = nameSamples, file = fSh, name = "nameSamples", level = 0)
  h5write(obj = nameMarks, file = fSh, name = "summedMarks", level = 0)
  te <- Sys.time()
  cat("Time cost:", as.numeric(te - ts, units = "mins"), "mins\n")
}
