# Fri Feb 11 18:17:02 2022 

# Title:
# Status: Complete
# Author: Tattini Lorenzo
# Modified: Nicolò Tellini

# Comments:

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

# Variables ----
# workDir <- "/home/ntellini/GITHUB"
argsVal <- commandArgs(trailingOnly = T)
workDir <- argsVal[1]

# Libraries ----

library(vcfR,verbose = F)

# body ----

outDir <- file.path(workDir, "mrk")

dir.create(outDir, showWarnings = F)

mrkdir <- file.path(workDir, "mrk")

vcffiles <- list.files(mrkdir, pattern = "\\.vcf\\.gz$")

outDirRes <- file.path(outDir,basename(mrkdir))

for (ind in vcffiles) {
  
  myVcf <- read.vcfR(file.path(mrkdir, ind), 
                     nrows = -1,convertNA = T, verbose = T)
  
  refPath <- unlist(strsplit(grep("reference", myVcf@meta, value = T), 
                             split = "##reference=file://", fixed = T))[2]
  
  refName <- unlist(strsplit(basename(refPath), split = "\\."))[[1]]
  
  myVcf@fix[, 1] <- gsub(pattern = "_.*$", "", myVcf@fix[, 1])
  
  sample <- unlist(strsplit(ind, split = "\\."))[1]
  
  outDirSamp <- file.path(outDir, sample)
  
  dir.create(outDirSamp, recursive = T, showWarnings = F)
  
  if ( length(grep("/",sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",1))) > 0 ) {
    
  } else {
    GP <- paste0(sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",1),"/",sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",1))
    PL <- sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",2)
    DP <- sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",3)
    SP <- sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",4)
    ADF <- sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",5)
    ADR <- sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",6)
    AD <- sapply(strsplit(as.character(as.data.frame(myVcf@gt)[,2]),split = ":"),"[[",7)
    myVcf@gt[,2] <- paste(GP,PL,DP,SP,ADF,ADR,ADF,AD,sep = ":")
  }
  
  # save RData
  outRData <- file.path(outDirSamp, paste(sample, refName, "vcf", "RData", sep = "."))
  vcfData <- list(myVcf@meta, myVcf@fix, myVcf@gt)
  names(vcfData) <- c("meta", "fix", "format")
  save(list = c("vcfData"), file = outRData)
}