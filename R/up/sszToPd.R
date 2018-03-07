tms <- Sys.time()
source("R/functions/getMetadata.R")
byHand <- FALSE
if (byHand) {
  #### get args by hand
  whichData <- "roadmap"
  getMetadata(whichData)
  whichChrs <- seq_along(nameChrs)
  # # whichChrs <- 21
} else {
  #### get args from shell
  args <- commandArgs(trailingOnly = TRUE)
  whichData <- args[1]
  getMetadata(whichData)
  whichChrs <- as.integer(args[-1])
  stopifnot(length(setdiff(whichChrs, seq_along(nameChrs))) == 0L)
}

dirSsz <- file.path(dirDataHome, paste0("nlp", "_c_z_ssz"))
prefixSsz <- "ssz"
source("R/functions/h5Utils.R")
fSsz <- getH5FilesWithChrName(dirSsz, prefixSsz, nameChrs)

dirStats <- file.path(dirDataHome, "stats")
prefixStats <- "nlpToSsz"
stats <- readRDS(file.path(dirStats, paste0(prefixStats, "_", "whole", ".rds")))
sumSquaredSszOfSample <- stats$sumSquaredSszOfSample

#### get normalization factor
twoNormSszOfSample <- sqrt(sumSquaredSszOfSample)
constTwoNorm <- 1000
factorTwoNorm <- constTwoNorm / twoNormSszOfSample

## chunk of pd h5
chunk <- c(1e5, 1)
## if level == 0, then chunk will be ignored by h5createDataset()
level <- 7
sizePartForChr <- 1e6
dir2n <- paste0(dirSsz, "_2n")
# dir2n <- NULL
dirZ <- paste0(dirSsz, "_2n_z")
# dirZ <- NULL
dirPd <- paste0(dirSsz, "_2n_z_pd")
source("R/functions/sszToPdByChr.R")

for (iChr in whichChrs) {
  chr <- nameChrs[iChr]
  fSszInChr <- fSsz[iChr]
  sszToPdByChr(chr, fSszInChr,
               nameSamples, nameMarks,
               factorTwoNorm, sizePartForChr,
               dirPd,
               chunk, level,
               dir2n, dirZ)
}
tme <- Sys.time()
cat("Total time:", as.numeric(tme - tms, units = "mins"), "mins\n")
