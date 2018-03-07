source("R/functions/getMetadata.R")
whichData <- "roadmap"
getMetadata(whichData)
sizeChrs <- seqlengths(seqInfo)[nameChrs]


# get bw of one score -----------------------------------------------------

suffix <- "nlp_c_z_ssz"
dirInput <- file.path(dirDataHome, suffix)
prefixInput <- "ssz"
source("R/functions/h5Utils.R")
fInput <- getH5FilesWithChrName(dirInput, prefixInput, nameChrs)
dirOutput <- file.path(dirDataHome, "bigwig", prefixInput)
if (!file.exists(dirOutput)) dir.create(dirOutput, recursive = TRUE)
source("R/functions/io.R")
nameSamplesENCODE <- paste0("E", 114:129)
nameSamplesRoadmap <- setdiff(nameSamples, nameSamplesENCODE)
nameSamplesCsd <- paste0("E", c("003", "123", "116", "118", "122",
                                "120", "128", "127", "119"))
whichSamples <- nameSamplesCsd
for (sample in whichSamples) {
  message(sample, "_", prefixInput, ".bw")
  rleListSample <- getRleList(fInput, sample, nameSamples, nameMarks, nameChrs,
                              sizeBin, sizeChrs)
  fBwSample <- file.path(dirOutput, paste0(sample, "_", prefixInput, ".bw"))
  export.bw(object = rleListSample, con = fBwSample)
}


# get bw of scores of marks -----------------------------------------------

suffix <- "nlp_c_z"
dirInput <- file.path(dirDataHome, suffix)
prefixInput <- "z"
source("R/functions/h5Utils.R")
fInput <- getH5FilesWithChrName(dirInput, prefixInput, nameChrs)
dirOutput <- file.path(dirDataHome, "bigwig", prefixInput)
if (!file.exists(dirOutput)) dir.create(dirOutput, recursive = TRUE)

whichSamples <- nameSamplesCsd
whichMarks <- nameMarks
for (sample in whichSamples) {
  for (mark in whichMarks) {
    message(sample, "_", mark, "_", prefixInput, ".bw")
    rleListSampleMark <-
      getRleListOneMark(fInput, sample, mark,
                        nameSamples, nameMarks, nameChrs,
                        sizeBin, sizeChrs)
    fBwSampleMark <- file.path(dirOutput, paste0(sample, "_", mark, "_", prefixInput, ".bw"))
    export.bw(object = rleListSampleMark, con = fBwSampleMark)
  }
}


# get bed -----------------------------------------------------------------

dirInput <- file.path(dirDataHome, "csre")
dirOutput <- file.path(dirDataHome, "bed")
if (!file.exists(dirOutput)) dir.create(dirOutput, recursive = TRUE)
nameAllCsres <- dir(dirInput)
eid <- data.table::fread("data/EID_metadata.tab", data.table = FALSE)
eidColors <- unique(eid[, c("COLOR", "EID")])
eidColors <- setNames(eidColors$COLOR, eidColors$EID)
for (nameCsre in nameAllCsres) {
  # nameCsre <- "csre_nlp_c_z_ssz_2n_z_pd_sh_0.1_60.txt"
  fileCsre <- file.path(dirInput, nameCsre)
  source("R/functions/csreGetAndConvert.R")
  csreSignalGR <-
    getCsreGRFromFile(fileCsre,
                      "gr", nameChrs, nameSamples, seqInfo)
  names(mcols(csreSignalGR)) <- "name"
  csreSignalGR$itemRgb <- eidColors[as.character(csreSignalGR$name)]
  # csreSignalGR$mnemonic <- eid[decode(match(csreSignalGR$name, eid$EID)),
  #                              "MNEMONIC"]
  # csreSignalGR$group <- eid[decode(match(csreSignalGR$name, eid$EID)),
  #                              "GROUP"]
  csreSignalGR <- sort(csreSignalGR)
  fBed <- file.path(dirOutput, sub(".txt", ".bed", nameCsre))
  # fTabix <- file.path(dirOutput, sub(".txt", "", nameCsre))
  message(fBed)
  export.bed(csreSignalGR, fBed, index = FALSE)
  export.bed(csreSignalGR, fBed, index = TRUE)
  # exportToTabix(csreSignalGR, fBed)
}

