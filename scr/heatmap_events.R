# Thu Jun 18 09:47:02 2020 

# Title: Heatmap AllSegments
# Author: Nicolò T.
# Status: complete

# Comments:

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F,scipen=999)

# Variables ----

argsVal <- commandArgs(trailingOnly = T)

BaseDir <- argsVal[1]
refBonLabel <- argsVal[2]
refBon <- argsVal[2]
setwd(BaseDir)

# comment the lines below ----
# 
# BaseDir <- "/home/ntellini/GITHUB"
# refBonLabel <- "Scc"
# refBon <- "Scc"
# setwd(BaseDir)

# Libraries ----

library(ggplot2)
library(data.table)

# body ----


allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

allRData <- list.files(paste0(BaseDir,"/int/AllSegments"),pattern = ".RData")
Events <- data.frame()
for (ind in allRData) {
  print(ind)
  load(paste0(BaseDir,"/int/AllSegments/",ind))
  allEv <- allEv[allEv$color != "blue" & allEv$strain == "Scc" & allEv$len > 1,]
  allEv <- allEv[,-c(8:13)]
  allEv <- allEv[,-c(11)]
  Events <- rbind(Events,allEv)
  rm(allEv)
}

allEv0 <- data.frame()
tabplot <- c()
samples <- c()
for (indC in allRData) {
  
  print(indC)
  # load dataframe ----
  load(paste0(BaseDir,"/int/AllSegments/",indC))
  
  # select events != blue, CBS != che è l'hybrid (fa rumore qui) ----
  allEv <- allEv[allEv$color != "blue" & allEv$strain == "Scc" & allEv$len > 1,]
  
  allEv$clade <- "" 
  
  if (nrow(allEv) > 0) {
    
    samples <- c(samples,as.vector(allEv$hybrid))
    
    # tieni solo chr, fisrt , last e Scc ----
    allEv <- allEv[,c(1,3,4,14)]
    
    allEv$Evlength <- allEv$last - allEv$first
    allEv$chr <- as.factor(allEv$chr)
    allEv <- allEv[order(allEv$Evlength,decreasing = T),]
    
    # frammenta e conta----
    counts <- rep(0,times=nrow(allEv))
    allEv$counts <- counts
    rownames(allEv) <- 1:nrow(allEv)
    
    allEv0temp <- allEv[allEv[,"Evlength"] == 0,]
    allEv <- allEv[!(allEv[,"Evlength"] == 0),]
    if ( nrow(allEv) > 0 ) {
      pointsall <- c()
      for (indE in 1:nrow(allEv)) {
        pointstemp <-  seq(from = allEv[indE,"first"],to = allEv[indE,"last"],by = 1)
        pointsall <- c(pointsall,pointstemp)
      }
      pointsall <- c(pointsall)
      
      allEv0 <- rbind(allEv0temp,allEv0)
      
      tab <- as.data.frame(table(pointsall))
      chrs <- rep(sapply(strsplit(indC,split = "\\."),"[[",2) ,times=nrow(tab))
      tab <- as.data.frame(cbind(chrs,tab))
      tab$pointsall <- as.numeric(as.vector(tab$pointsall))
      tabplot <- rbind(tabplot,tab)
      rm(tab,pointsall)
    }
  }
}

tabplot$chrs <- as.factor(tabplot$chrs)

centromeric <- read.delim(paste0(BaseDir,"/rep/Ann/Scc.centromere.txt"), header=FALSE)
chrlen <- read.delim(paste0(BaseDir,"/rep/Ann/Scc.chrs.txt"), header=FALSE)

if (nrow(chrlen) == 17) {
  chrlen <- chrlen[-nrow(chrlen),]
}    

xmin <- rep(0,times=nrow(chrlen))
chrlen$xmin <- xmin
colnames(chrlen)[1:2] <- c("chr","xmax")

yminimi <- seq(from=0, to=7.5, by=0.5)
ymassimi <- seq(from=0.25, to=7.75, by=0.5)

chrlen$ymin <- yminimi
chrlen$ymax <- ymassimi

centromeric$ymin <- yminimi
centromeric$ymax <- ymassimi

yminimo <- rep("",times=nrow(tabplot))
ymassimo <- rep("",times=nrow(tabplot))

tabplot$ymin <- yminimo
tabplot$ymax <- ymassimo

for (ind in 1:length(allChr)) {
  tabplot[tabplot[,1] == chrlen[ind,"chr"],"ymin"] <- chrlen[chrlen[,1] == chrlen[ind,"chr"],"ymin"]
  tabplot[tabplot[,1] == chrlen[ind,"chr"],"ymax"] <- chrlen[chrlen[,1] == chrlen[ind,"chr"],"ymax"]
}


tabplot$chrs <- factor(tabplot$chrs, levels=allChr)
tabplot <- tabplot[order(match(tabplot[[1]], allChr)), ]
tabplot$ymax <- as.numeric(tabplot$ymax)
tabplot$ymin <- as.numeric(tabplot$ymin)


tabplot2 <- data.frame()
for (indC in unique(tabplot$chrs) ) {
  temptab <- tabplot[tabplot[,1] == indC,]
  conte <- rle(temptab$Freq)[["lengths"]]
  print(indC)
  start <- 0
  end <- 0
  taBB <- data.frame()
  for ( ind in 1:length(conte) ) {
    print(ind)
    if ( ind == 1 ) { 
      inizio <- start + 1
      fine <- end + conte[ind]
    }
    
    tab1 <- temptab[inizio:fine,]
    chr <- unique(tab1$chrs)
    freq <- unique(tab1$Freq)
    
    if ( length(unique( diff(tab1$pointsall) ) ) != 1 & length(diff(tab1$pointsall)) != 0 ) {
      
      tab2 <- unname(tapply(tab1$pointsall, cumsum(c(1, diff(tab1$pointsall)) != 1), range))
      
      tab3 <- data.frame()
      for (ind2 in 1:length(tab2)) {
        tab4 <-  data.frame(tab2[[ind2]][1],tab2[[ind2]][2])
        tab3 <- rbind(tab3,tab4)
      }
      
      tab3$chr <- chr
      tab3$freq <- freq
      colnames(tab3)[1:2] <- c("start","end")
      
      if (nrow(taBB) != 0) {
        colnames(taBB)<- c("start","end","chr","freq")
      }
      
      taBB <- rbind(taBB,tab3)
    } else {
      riga <- c(min(tab1$pointsall),max(tab1$pointsall),as.character(tab1[1,1]),unique(tab1$Freq))
      taBB <- rbind(taBB,riga)
    }
    
    if (ind != length(conte)) {
      inizio <- fine + 1
      fine <- fine + (conte[ind+1]) 
    }
  }
  colnames(taBB)<- c("start","end","chr","freq")
  tabplot2 <- rbind(tabplot2,taBB)
}

tabplot2[,1] <- as.numeric(tabplot2[,1])
tabplot2[,2] <- as.numeric(tabplot2[,2])
tabplot2[,4] <- as.numeric(tabplot2[,4])

tabplot2$length <- tabplot2$end -  tabplot2$start
tabplot2 <- tabplot2[order(tabplot2$length,decreasing = T),]

tabplot2$ymin <- rep("",times=nrow(tabplot2))
tabplot2$ymax <- rep("",times=nrow(tabplot2))

for (indC in allChr) {
  tabplot2[tabplot2[,3] == indC,"ymin"] <- chrlen[chrlen[,1] == indC,"ymin"]
  tabplot2[tabplot2[,3] == indC,"ymax"] <- chrlen[chrlen[,1] == indC,"ymax"]
}

# Plot part ----

# Rescale ----

tabplot2[,1] <- tabplot2[,1] / 1000
tabplot2[,2] <- tabplot2[,2] / 1000
centromeric[,c(1,2)] <- centromeric[,c(1,2)] / 1000
chrlen[,c(2,3)] <- chrlen[,c(2,3)] /1000

tabplot2[,6] <- as.numeric(tabplot2[,6])
tabplot2[,7] <- as.numeric(tabplot2[,7])

pHeatRainbow <- ggplot(tabplot2) + 
  geom_rect(tabplot2,mapping = aes(xmin=start, xmax=end, ymin=ymin,ymax=ymax,fill=freq)) +
  geom_rect(chrlen, mapping = aes(xmin=0, xmax=xmax, ymin=ymin,ymax=ymax),fill="grey99",color="black", size=.1,alpha=0.000001) +
  scale_fill_distiller( type = "seq", palette = 8,direction = 1) +
  geom_point(centromeric,mapping = aes(x =((as.numeric(V1) + as.numeric(V2))/2), y=((as.numeric(ymin) + as.numeric(ymax))/2)),color="black",fill="black") +
  annotate(geom="text", x=as.numeric(-0.05), y=((as.numeric(yminimi) + as.numeric(ymassimi))/2), label=allChr ,color="black",hjust=1.5,size=3) +
  scale_x_continuous(name="Genomic Position (Kb)") +
  theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        ) + 
  labs(fill = "samples fraction",title = "Heatmap")

pPath <- file.path(paste0(BaseDir,"/int/Allchr.events.heatmap.pdf"))
pdf(file = pPath, width = 16, height = 10)
print(pHeatRainbow)
dev.off()

fwrite(tabplot2,paste0(BaseDir,"/int/Allchr.events.heatmap.txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)