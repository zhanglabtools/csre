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
## input can be count, nlp or any track else
dirInput <- dirNlp
prefixInput <- "nlp"
source("R/functions/h5Utils.R")
fInput <- getH5FilesWithChrName(dirInput, prefixInput, nameChrs)

sizePartForChr <- 1e6
dirCtd <- file.path(dirDataHome, paste0(prefixInput, "_c"))
# dirCtd <- NULL
dirZscore <- file.path(dirDataHome, paste0(prefixInput, "_c_z"))
# dirZscore <- NULL
dirSsz <- file.path(dirDataHome, paste0(prefixInput, "_c_z_ssz"))
dirStats <- file.path(dirDataHome, "stats")
thresInput <- thresNlp
maxInputOfMarks <- maxNlpOfMarks
## chunk of ssz h5
chunkMat <- c(1e5, 1)
## if level == 0, then chunk will be ignored by h5createDataset()
levelMat <- 7
chunkTsr <- c(1e5, 1, 1)
levelTsr <- 7
source("R/functions/inputToSszByChr.R")

for (iChr in whichChrs) {
  chr <- nameChrs[iChr]
  fInputInChr <- fInput[iChr]
  inputToSszByChr(chr, fInputInChr,
                  nameSamples, nameMarks,
                  thresInput, maxInputOfMarks,
                  sizePartForChr,
                  dirSsz = dirSsz, dirStats = dirStats,
                  chunkMat, levelMat, prefixInput,
                  dirCtd = dirCtd, dirZscore = dirZscore,
                  chunkTsr, levelTsr)
}
tme <- Sys.time()
cat("Total time:", as.numeric(tme - tms, units = "mins"), "mins\n")
