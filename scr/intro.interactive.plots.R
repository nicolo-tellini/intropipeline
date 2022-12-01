# # # # # # # # # # # # # # # # # # # # # # 
# Plottare tutti i cromosomi per campione #
# # # # # # # # # # # # # # # # # # # # # # 

# Status: complete

# Opening ---- 

rm(list = ls())
gc()
gcinfo(FALSE)
options(warn = 1)
options(stringsAsFactors = F)
options(scipen=999)

# Dirs ----

argsVal <- commandArgs(trailingOnly = T)
BaseDir <- argsVal[1]
ref1 <- "Scc"

BaseDir <- "/home/ntellini/GITHUB"

# Libraries ----

library(ggplot2)

# Table ----

ListFiles <- list.files(path =  paste0(BaseDir,"/int/AllSegments"), pattern = ".RData" )

LOHTab <- data.frame()
for (ind in ListFiles) {
  load(file = paste0(BaseDir,"/int/AllSegments/",ind))
  LOHTab <- as.data.frame(rbind(LOHTab,allEv))
  rm(allEv)
}

centromeric <- read.delim(paste0(BaseDir,"/rep/Ann/Scc.centromere.txt"), header=FALSE)
chrlen <- read.delim(paste0(BaseDir,"/rep/Ann/Scc.chrs.txt"), header=FALSE)

if (nrow(chrlen) == 17) {
  chrlen <- chrlen[-nrow(chrlen),]
}    

samples <- as.character(unique(LOHTab$hybrid))

yminimi <- seq(from=0, to=7.5, by=0.5)
ymassimi <- seq(from=0.25, to=7.75, by=0.5)

emptytab <- as.data.frame(matrix(ncol=2,nrow = nrow(LOHTab)))
LOHTab <- cbind(LOHTab,emptytab)
colnames(LOHTab)[18:19] <- c("ymin","ymax")

z <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII","chrIX","chrX","chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

LOHTab$chr <- factor(LOHTab$chr, levels=z)

LOHTab <- LOHTab[order(match(LOHTab[[1]], z)), ]

ymin <- 0
ymax <- 0.25
for (indc in z) {
  
  LOHTab[LOHTab[,1] == indc,"ymin"] <- ymin
  LOHTab[LOHTab[,1] == indc,"ymax"] <- ymax
  
  ymin <- ymin + 0.5
  ymax <- ymax + 0.5
}
rm(ymin,ymax,emptytab)

LOHTab$color <- as.character(LOHTab$color)
LOHTab[LOHTab$color == "blue","color"] <- "skyblue"
#LOHTab <- LOHTab[LOHTab$len>=10,]
LOHTab[,15] <- as.character(LOHTab[,15])

for (indS in samples) {
  
  temptab <- LOHTab[LOHTab[,15] == indS & LOHTab[,14] == ref1,]
  temptab$transparency <- 1
  
  temptab[temptab$color != "skyblue","transparency"] <- temptab[temptab$color != "skyblue","len"] / max(temptab[temptab$color != "skyblue","len"])
  
  colnames(temptab)[7] <- "num_mrk"
  colnames(temptab)[11] <- "mrk_den"
  colnames(temptab)[1] <- "chr"
  
  plotallchr <-  ggplot() + 
    geom_rect(temptab,mapping = aes(xmin=first, xmax=last, ymin=ymin,ymax=ymax,"markers"=num_mrk,"mrkdensity"=mrk_den,"chr"=chr,"start"=first,"end"=last),fill=temptab$color) +
    #geom_rect(chrlen, mapping = aes(xmin=0, xmax=V2, ymin=yminimi,ymax=ymassimi),fill="grey99",color="black",size=0.5,alpha=0) +
    ggtitle(paste0(indS)) +
    geom_point(centromeric,mapping = aes(x = ((as.numeric(V1) + as.numeric(V2))/2)), y=((as.numeric(yminimi) + as.numeric(ymassimi))/2),color="black",fill="black") +
    annotate(geom="text", x=as.numeric(-50), y=((as.numeric(yminimi) + as.numeric(ymassimi))/2), label=z ,color="black",hjust=1.5) +
    xlim(-50,chrlen[2,2]) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position="none")
  
  temptab2 <- temptab[temptab$chr == "chrII",]
  centromeric <- centromeric[2,]
  annotation <- read.table("/home/ntellini/GITHUB/rep/Ann/Scc.all_feature.gff")
  annotation <- annotation[annotation$V1 == "chrII" & annotation$V3 == "gene",]
  
  tab <- as.data.frame(matrix(nrow = nrow(annotation), ncol = ncol(temptab2)))
  colnames(tab) <- colnames(temptab2)
  tab$chr <- annotation$V1
  tab$first <- annotation$V4
  tab$last <- annotation$V5
  tab$hybrid <- sub("Name=","",sapply(strsplit(annotation$V9,";"),"[[",2))

  j=1
  ymin <- c(1.75,2.1)
  ymax <- c(2,2.35)
  for (i in 1:nrow(tab)) {
    tab[i,"ymin"] <- ymin[j]
    tab[i,"ymax"] <- ymax[j]
    j=j+1
    if (j>2) {
      j=1
    }
  }
  
  tab$color <- "darkolivegreen4"
  tab$contorno <- "black"
  tab$ymin <- tab$ymin -0.8
  tab$ymax <- tab$ymax -0.8 - 0.12
  temptab2 <- rbind(temptab2,tab)
  
  tab$start <- NULL
  tab$end <- NULL
  colnames(tab)[2] <- "start"
  colnames(tab)[3] <- "end"
  colnames(tab)[13] <- "gene"
  
  p <- ggplot() + 
    geom_rect(temptab2,mapping = aes(xmin=first, xmax=last, ymin=ymin,ymax=ymax,"markers"=num_mrk,"mrkdensity"=mrk_den,"chr"=chr,"start"=first,"end"=last),fill=temptab2$color) +
    geom_rect(tab,mapping = aes(xmin=start, xmax=end, ymin=ymin,ymax=ymax,"chr"=chr,"start"=start,"end"=end,"gene"=gene),fill=tab$color,color=tab$contorno,size=0.1,alpha=0.5) +
    # geom_rect(chrlen, mapping = aes(xmin=0, xmax=V2, ymin=yminimi,ymax=ymassimi),fill="grey99",color="black",size=0.5,alpha=0) +
    # ggtitle(paste0(indS)) +
    #geom_point(centromeric,mapping = aes(x = ((as.numeric(V1) + as.numeric(V2))/2)), y=((as.numeric(yminimi) + as.numeric(ymassimi))/2),color="black",fill="black") +
    annotate(geom="text", x=as.numeric(-50), y=c(0.5+0.75)/2, label="chrII" ,color="black",hjust=1.5) +
    xlim(-50,chrlen[2,2]) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position="none")

  
  
  plot <- ggplotly(p,dynamicTicks = T,height = 500)
  saveWidget(plot, "/home/ntellini/GITHUB/p1.html", selfcontained = T, libdir = "lib")
  
  pPath <- file.path(paste0(BaseDir,"/int/",indS,".allchr.events.pdf"))
  # dir.create(paste0(BaseDir,"/LOH/",indS),showWarnings = F,recursive = T)
  pdf(file = pPath, width = 16, height = 10)
  print(plotallchr)
  dev.off()
  
  pPath <- file.path(paste0(BaseDir,"/int/",indS,".allchr.events.jpeg"))
  jpeg(file = pPath, width = 1200, height = 800,quality = 100,res=100)
  print(plotallchr)
  dev.off()
}


