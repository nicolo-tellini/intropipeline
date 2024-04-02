# Fri Feb 11 18:17:02 2022 

# Title:
# Status: Complete
# Author: Nicolò Tellini

# Comments:

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

# Variables ----

# baseDir <- "/home/nico/intropipeline2"
argsVal <- commandArgs(trailingOnly = T)
baseDir <- argsVal[1]

# Libraries ----

library(data.table)
library(R.utils)

# Functions ----

str_split_custom <- function(x,y,z) sapply(strsplit(x,split = z),"[[",y) # frequent strsplit

# body ----

mrkdir <- file.path(baseDir, "mrk")

mrkfiles <- list.files(path =mrkdir,pattern = "^Sites.")

refmrk <- str_split_custom(mrkfiles,2,"\\.")

mrksites <- lapply(mrkfiles, function(x) fread(paste0(mrkdir,"/",x),sep = "\t",data.table = T))

mrksites <- lapply(mrksites,function(x) { x[,"rank"] <- 1:nrow(x); return(x)})

mrksites <- lapply(mrksites, function(x) x[,V1:=.(gsub("_.*","",x[[1]]))])

refs <- str_split_custom(mrkfiles,2,"\\.")

mrksites <- mapply(function(x,y){ x[,lab:=paste0(x[[1]],x[[2]],y)]; return(x)},mrksites,refs,SIMPLIFY = F)

mrksites <- lapply(mrksites, function(x) {cols = names(x)[c(1:3)]; x[,(cols):=NULL]; return(x)} )
             
vcfs_names <- list.files(mrkdir, pattern = "\\.vcf\\.gz$")

names(mrksites) <- refs

for (i in vcfs_names) {
  
vcfs <- fread(paste0(mrkdir,"/",i),sep = "\t",data.table = T,drop = c(3,6,7,8,9))

vcfs_refs <- sapply(strsplit(i,split = "\\."),"[[",2)

colnames(vcfs)[1] <- "chr"

vcfs[,chr:=.(gsub("_.*","",vcfs[[1]])),]

all_samples <- str_split_custom(i,1,"\\.")

colnames(vcfs)[5] <- all_samples

ref <- str_split_custom(i,2,"\\.")

vcfs[,lab:=.(paste0(vcfs[[1]],vcfs[[2]],ref)),]

cols <- colnames(vcfs)[5]

vcfs[ , (cols) := lapply(.SD, function(x) str_split_custom(x,1,"\\:") ), .SDcols = cols]

vcfs <- vcfs[vcfs[[5]] != "./."]

vcfs[, rank:=.(match(lab,mrksites[[ref]][["lab"]])),]

vcfs[,multi:=nchar(ALT)]

vcfs <- vcfs[vcfs[["multi"]] == 1]

cols = names(vcfs)[c(6,8)]

vcfs[, (cols):=NULL]

out <- paste0(mrkdir,"/",all_samples,".",vcfs_refs,".tab.gz")

fwrite(x = vcfs,file = out,append = F,quote = F,sep = "\t",row.names = F,col.names = T,nThread = 2)

}
