# Thu Jul  9 15:30:41 2020 

# Author: Nicolò T.
# Status: Complete

# Comments:

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

# Variables ----

# BaseDir <- '/home/ntellini/GITHUB'
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

if ( nrow(temptab) > 0 ) {
  temptab[,2] <- sapply(strsplit(temptab[,2],split = "-"),"[[",2)
}

finaltab <- rbind(temptab1,temptab)

rm(temptab,temptab1,tab)

write.table(finaltab,append = F,file = paste0(BaseDir,"/seq/read.len"),quote = FALSE,sep = "\t",row.names = F,col.names = F)