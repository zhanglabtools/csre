tms <- Sys.time()
source("R/functions/getMetadata.R")
byHand <- FALSE
isBsub <- TRUE
if (byHand) {
  #### get args by hand
  whichData <- "roadmap"
  # regionHeight <- 1.4907
  regionHeight <- 0.1
  regionWidth <- 60
} else {
  #### get args from shell
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 3L)
  whichData <- args[1]
  regionHeight <- as.double(args[2])
  regionWidth <- as.double(args[3])
  stopifnot(!is.na(regionHeight), !is.na(regionWidth))
}
getMetadata(whichData)

suffix <- "nlp_c_z_ssz_2n_z_pd_sh"
# suffix <- "nlp_b_h_sah_2n_z_pd_sh"
dirInput <- file.path(dirDataHome, suffix)
prefixInput <- "sh"
source("R/functions/h5Utils.R")
fInput <- getH5FilesWithChrName(dirInput, prefixInput, nameChrs)

dirOutput <- file.path(dirDataHome, "csre")
dirBed <- file.path(dirDataHome, "bed")

library(parallel)
numMaxCores <- Inf
numChrs <- length(nameChrs)
if (!isBsub) {
  numAvailClusters <- detectCores()
  numCores <- min(numChrs, numAvailClusters, numMaxCores)
  cat("make cluster of ", numCores, " nodes on localhost\n", sep = "")
  cl <- makeCluster(numCores)
} else {
  hosts <- strsplit(Sys.getenv("LSB_HOSTS"), split = " ")[[1]]
  usedHosts <- hosts[-1]
  if (length(usedHosts) > numChrs)
    usedHosts <- usedHosts[seq_len(numChrs)]
  numCores <- length(usedHosts)
  cat("make cluster of ", numCores, " nodes on hosts:\n", sep = "")
  cat(usedHosts, "\n", sep = " ")
  cl <- makeCluster(usedHosts)
}
getCsreByChr <- function(iChr, nameChrs, fInput,
                         nameSamples, nameMarks,
                         regionHeight, regionWidth,
                         dirWd = getwd()) {
  ### Once setwd() is called, it would have effect until another setwd() was
  ### called in this sub session. When dirWd is missing, dirWd would be getwd()
  ### in the *current* sub session, so the wd would change in that case.
  ### However, my aim to keep dirWd as an arg is to change wd of the sub session
  ### to a proper path to make source find the checkH5.R.
  setwd(dirWd)
  cat("The current working directory of this sub session is:\n", dirWd, "\n",
      sep = "")
  library(rhdf5)
  ### If cl is made by host names that are not localhost, the default wd of sub
  ### sessions are home dir (~). So the source() would not find checkH5.R, and
  ### an error occurs. So I pass dirWd from the master session to set the sub
  ### session wd to dmge home dir. Thus the source() would succeed!
  source("R/functions/checkH5.R")
  chr <- nameChrs[iChr]
  fInputInChr <- fInput[iChr]
  numSamples <- length(nameSamples)
  cat(chr, "\n", sep = "")
  #### check input
  cat("check input h5\n")
  checkMatOneChrH5(fInputInChr, chr, nameSamples, nameMarks)
  numBinsInChr <- readNumBinsInMatOneChrH5(fInputInChr)
  #### read data
  cat("read data of ", numBinsInChr, " and ", numSamples, " samples\n")
  matInputInChr <- h5read(file = fInputInChr, name = chr, index = NULL)

  csreList <- list()
  #### WigToBedLike
  cat("find CSREs in samples\n")
  for (iSample in seq_len(numSamples)) {
    v <- matInputInChr[, iSample] > regionHeight
    regionStart <- which(diff(c(0, v)) == 1)
    regionEnd <- which(diff(c(v, 0)) == -1)
    ## both start and end are inclusive
    regionLength <- regionEnd - regionStart + 1L
    idx <- regionLength >= regionWidth
    regionStart <- regionStart[idx]
    regionEnd <- regionEnd[idx]
    regionLength <- regionLength[idx]
    if (length(regionStart) > 0) {
      cat(nameSamples[iSample], "\n", sep = "")
      csreList[[length(csreList) + 1L]] <-
        data.frame(Cell_Type = nameSamples[iSample],
                   Chromosome = chr,
                   Start_Bin = regionStart,
                   End_Bin = regionEnd,
                   Length_Bin = regionLength,
                   stringsAsFactors = FALSE)
    }
  }
  rm(matInputInChr)
  gc()
  cat("\n")
  csreByChr <- do.call(rbind, csreList)
  rm(csreList)
  gc()
  return(csreByChr)
}
cat("get CSRE on cluster of ", numCores, " nodes\n", sep = "")
### The following getwd() in clusterApplyLB() would be called before sent to
### each sub session, you can see that by source code of clusterApplyLB(). So it
### would be the wd of master session. That's what I want to make getCsreByChr()
### source checkH5.R in the right place.
wholeCsreList <- clusterApplyLB(cl, seq_len(numChrs), getCsreByChr,
                                nameChrs, fInput,
                                nameSamples, nameMarks,
                                regionHeight, regionWidth,
                                getwd())
stopCluster(cl)
csre <- do.call(rbind, wholeCsreList)
rm(wholeCsreList)
## reorder the factor levels
cat("reorder the factor levels\n")
csre$Cell_Type <- factor(csre$Cell_Type, levels = nameSamples)
csre$Chromosome <- factor(csre$Chromosome, levels = nameChrs)
## arrange the table
cat("arrange the table\n")
csre <- csre[order(csre$Cell_Type, csre$Chromosome,
                   csre$Start_Bin, csre$Length), ]
## plus zerobased intervals
cat("plus zerobased intervals\n")
csre$Start_ZeroBased <- (csre$Start_Bin - 1L) * sizeBin
csre$End_ZeroBased <- csre$End_Bin * sizeBin
csre$Length <- csre$Length_Bin * sizeBin
## keep End_ZeroBased no more than end bp of each chromosome
## and keep Start_ZeroBased strictly less than end bp of each chromosome
cat("keep intervals bounded in chr\n")
sizeChrs <- GenomeInfoDb::seqlengths(seqInfo)[nameChrs]
upperbound <- sizeChrs[match(csre$Chromosome, names(sizeChrs))]
if (any(csre$End_ZeroBased > upperbound)) {
  cat("Number of intervals to be bounded:",
      sum(csre$End_ZeroBased > upperbound),
      "\n")
  ## csre$Start_ZeroBased can not go out of bound, so the below is not necessary
  # csre$Start_ZeroBased <- pmin(csre$Start_ZeroBased, upperbound - 1)
  csre$End_ZeroBased <- pmin(csre$End_ZeroBased, upperbound)
  csre$Length <- csre$End_ZeroBased - csre$Start_ZeroBased
}

## write table
cat("write table of CSREs of all chromosomes and samples\n")
if (!file.exists(dirOutput)) dir.create(dirOutput, recursive = TRUE)
fCsre <- file.path(dirOutput,
                   paste0("csre", "_", suffix, "_",
                          regionHeight, "_", regionWidth, ".txt"))
write.table(csre, fCsre, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
## write bed
cat("write bed of CSREs of all chromosomes and samples\n")
if (!file.exists(dirBed)) dir.create(dirBed, recursive = TRUE)
library(magrittr)
csreForBed <- csre %>%
  dplyr::transmute(cellType = factor(Cell_Type, nameSamples),
                   chr = factor(Chromosome, nameChrs),
                   start = Start_ZeroBased, end = End_ZeroBased,
                   len = Length)
source("R/functions/csreGetAndConvert.R")
csreGR <- getCsreGR(csreForBed, seqInfo)
names(mcols(csreGR)) <- "name"
csreGR <- sort(csreGR)
if (grepl("roadmap", whichData)) {
  eid <- data.table::fread("data/EID_metadata.tab", data.table = FALSE)
  eidColors <- unique(eid[, c("COLOR", "EID")])
  eidColors <- setNames(eidColors$COLOR, eidColors$EID)
  csreGR$itemRgb <- eidColors[as.character(csreGR$name)]
}
fBed <- file.path(dirBed,
                  paste0("csre", "_", suffix, "_",
                         regionHeight, "_", regionWidth, ".bed"))
export.bed(csreGR, fBed, index = TRUE)

tme <- Sys.time()
cat("Total time:", as.numeric(tme - tms, units = "mins"), "mins\n")
