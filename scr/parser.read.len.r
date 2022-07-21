# Thu Jul  9 15:30:41 2020 

# Title:
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

 # BaseDir <- '/home/nico'
 # allDir <- list.dirs()
 # allFiles <- list.files()
  
 argsVal <- commandArgs(trailingOnly = T)
 BaseDir <- argsVal[1]

# Libraries ----

#  library()

# body ----
  
  tab <- read.table(paste0(BaseDir,"/seq/read.len"),header = F)
  
  temptab1 <- tab[grep("-",tab[,2],invert = T),] 
  temptab <- tab[grep("-",tab[,2]),]
  
  temptab[,2] <- sapply(strsplit(temptab[,2],split = "-"),"[[",2)

  finaltab <- rbind(temptab,temptab1)
  
  write.table(finaltab,append = F,file = paste0(BaseDir,"/seq/read.len"),quote = FALSE,sep = "\t",row.names = F,col.names = F)
