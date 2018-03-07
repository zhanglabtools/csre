tms <- Sys.time()
source("R/functions/getMetadata.R")
byHand <- FALSE
if (byHand) {
  #### get args by hand
  whichData <- "roadmap_1e5"
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
dirInput <- file.path(dirDataHome, "nlp_c_z_ssz_2n_z_pd")
prefixInput <- "pd"
source("R/functions/h5Utils.R")
fInput <- getH5FilesWithChrName(dirInput, prefixInput, nameChrs)

dirSh <- paste0(dirInput, "_sh")
## chunk of sh h5
chunk <- c(1e5, 1)
## if level == 0, then chunk will be ignored by h5createDataset()
level <- 7

for (iChr in whichChrs) {
  chr <- nameChrs[iChr]
  fInputInChr <- fInput[iChr]
  source("R/functions/shByChr.R")
  shByChr(chr, fInputInChr,
          nameSamples, nameMarks,
          dirSh,
          chunk, level,
          Inf, TRUE)
}
tme <- Sys.time()
cat("Total time:", as.numeric(tme - tms, units = "mins"), "mins\n")
