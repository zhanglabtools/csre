tms <- Sys.time()
source("R/functions/getMetadata.R")
byHand <- FALSE
if (byHand) {
  #### get args by hand
  whichData <- "roadmap"
} else {
  #### get args from shell
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1L)
  whichData <- args
}
getMetadata(whichData)

dirStats <- file.path(dirDataHome, "stats")
prefixStats <- "nlpToSsz"
numChrs <- length(nameChrs)
numSamples <- length(nameSamples)
numMarks <- length(nameMarks)
numSigBinsOfChrSampleMark <-
  array(0, dim = c(numChrs, numSamples, numMarks),
        dimnames = list(nameChrs, nameSamples, nameMarks))
numSigBinsOfChrMark <- array(0, dim = c(numChrs, numMarks),
                             dimnames = list(nameChrs, nameMarks))
sumSquaredSszOfChrSample <- array(0, dim = c(numChrs, numSamples),
                                  dimnames = list(nameChrs, nameSamples))
for (iChr in seq_along(nameChrs)) {
  chr <- nameChrs[iChr]
  fStatsInChr <- file.path(dirStats, paste0(prefixStats, "_", chr, ".rds"))
  statsInChr <- readRDS(fStatsInChr)
  stopifnot(identical(dimnames(numSigBinsOfChrSampleMark)[-1],
                      dimnames(statsInChr$numSigBinsOfSampleMarkInChr)),
            identical(colnames(numSigBinsOfChrMark),
                      names(statsInChr$numSigBinsOfMarkInChr)),
            identical(colnames(sumSquaredSszOfChrSample),
                      names(statsInChr$sumSquaredSszOfSampleInChr)))

  numSigBinsOfChrSampleMark[chr, , ] <- statsInChr$numSigBinsOfSampleMarkInChr
  numSigBinsOfChrMark[chr, ] <- statsInChr$numSigBinsOfMarkInChr
  sumSquaredSszOfChrSample[chr, ] <- statsInChr$sumSquaredSszOfSampleInChr
}
numSigBinsOfSampleMark <- colSums(numSigBinsOfChrSampleMark)
numSigBinsOfMark <- colSums(numSigBinsOfChrMark)
sumSquaredSszOfSample <- colSums(sumSquaredSszOfChrSample)

saveRDS(list(numSigBinsOfChrSampleMark = numSigBinsOfChrSampleMark,
             numSigBinsOfChrMark = numSigBinsOfChrMark,
             sumSquaredSszOfChrSample = sumSquaredSszOfChrSample,
             numSigBinsOfSampleMark = numSigBinsOfSampleMark,
             numSigBinsOfMark = numSigBinsOfMark,
             sumSquaredSszOfSample = sumSquaredSszOfSample),
        file.path(dirStats, paste0(prefixStats, "_", "whole", ".rds")))
tme <- Sys.time()
cat("Total time:", as.numeric(tme - tms, units = "mins"), "mins\n")
