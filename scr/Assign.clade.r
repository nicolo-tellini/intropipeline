# Mon Jul 13 14:24:20 2020 

# Title: Assign clade to results
# Author: Nicolò T.
# Status: Draft

# Comments:

# Options ----

  rm(list = ls())
  options(warn = 1)
  options(stringsAsFactors = F)
  gc()
  gcinfo(FALSE)
  options(scipen=999)
  
# Variables ----

# baseDir <- '/home/nico/pipeline'
# allDir <- list.dirs()
# allFiles <- list.files()
  
argsVal <- commandArgs(trailingOnly = T)
baseDir <- argsVal[1]
  
# Libraries ----

#  library()

# body ----

# load file names,ploidy and clade and assign the clade

tab_info <- read.table(paste0(baseDir,"/seq/ploidy"))
RData <- list.files(paste0(baseDir,"/LOH/AllSegments"),pattern = "RData")

for (indRData in RData) {
  load(paste0(baseDir,"/LOH/AllSegments/",indRData))
  allEv$clade <- rep("",times=nrow(allEv))
  allEv$hybrid <- as.character(allEv$hybrid)
  for (indH in unique(allEv$hybrid)) {
    clade <- tab_info[tab_info[,1] == indH,3]
    if (length(allEv[allEv[,15] == indH,18]) > 0) {
    allEv[allEv[,15] == indH,18] <- clade
    }
  }
  indC <- sapply(strsplit(indRData,split = "\\."),"[[",2)
  save(allEv, file = file.path(paste0(baseDir,"/LOH/AllSegments"), paste("Events.clades", indC, "RData", sep = ".")))
}