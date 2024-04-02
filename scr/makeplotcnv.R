options(scipen=999)
args <- commandArgs()

nChrom <- as.numeric(args[6])
lastArg <- 7 + nChrom - 1

dataTable <- read.table(args[5], header=TRUE);

ratio <- data.frame(dataTable)
ploidy <- type.convert(args[4])

# png(filename = paste(args[5],".log2.png",sep = ""), width = 5900, height = 5900,
#     units = "px", pointsize = 20, bg = "white", res = NA)
# plot(1:10)
# op <- par(mfrow = c(5,5))

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
  if (ratio$Ratio[i]>maxLevelToPlot) {
    ratio$Ratio[i]=maxLevelToPlot;
  }
}

allChr <- args[7:lastArg]
# allChr <- c("I_DBVPG6765", "II_DBVPG6765", "III_DBVPG6765", "IV_DBVPG6765", "V_DBVPG6765", "VI_DBVPG6765", "VII_DBVPG6765", "VIII_DBVPG6765", 
#             "IX_DBVPG6765", "X_DBVPG6765", "XI_DBVPG6765", "XII_DBVPG6765", "XIII_DBVPG6765", "XIV_DBVPG6765", "XV_DBVPG6765", "XVI_DBVPG6765")

for (i in allChr) {
  tt <- which(ratio$Chromosome==i)
  pdf(file = paste0(args[5],".chr",i,".log2.pdf"))
  if (length(tt)>0) {
    plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste("Position [bp] \n", paste0("chr",i)),ylab = "RC",pch = ".",col = colors()[88])
    tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
    points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136])
    
    tt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)	
    points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex=4)
    
    tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
    points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[461])
    tt <- which(ratio$Chromosome==i)
    
    #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
    #points(ratio$Start[tt],ratio$CopyNumber[tt], pch = ".", col = colors()[24],cex=4)
    
  }
  tt <- which(ratio$Chromosome==i)
  
  #UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:
  points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = ".", col = colors()[463],cex=4)
  dev.off()
}



