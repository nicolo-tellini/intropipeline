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

#BaseDir <- "/home/ntellini/TWO-assembly_pipe"

# Libraries ----

library(ggplot2)

# Table ----

ListFiles <- list.files(path =  paste0(BaseDir,"/Introgressions/AllSegments"), pattern = ".RData" )

LOHTab <- data.frame()
for (ind in ListFiles) {
  load(file = paste0(BaseDir,"/Introgressions/AllSegments/",ind))
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

LOHTab[,15] <- as.character(LOHTab[,15])

for (indS in samples) {
  
  temptab <- LOHTab[LOHTab[,15] == indS & LOHTab[,14] == ref1,]
 
  plotallchr <-  ggplot() + 
    geom_rect(chrlen, mapping = aes(xmin=0, xmax=V2, ymin=yminimi,ymax=ymassimi),fill="grey99",color="black", size=.1) +
    geom_rect(temptab,mapping = aes(xmin=first, xmax=last, ymin=ymin,ymax=ymax),fill=temptab$color,alpha=.5) +
    ggtitle(paste0(indS)) +
    geom_point(centromeric,mapping = aes(x = ((as.numeric(V1) + as.numeric(V2))/2)), y=((as.numeric(yminimi) + as.numeric(ymassimi))/2),color="black",fill="black") +
    annotate(geom="text", x=as.numeric(-50), y=((as.numeric(yminimi) + as.numeric(ymassimi))/2), label=z ,color="black",hjust=1.5) +
    xlim(-50,1532000) +
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
  
  pPath <- file.path(paste0(BaseDir,"/Introgressions/",indS,".allchr.events.pdf"))
  # dir.create(paste0(BaseDir,"/LOH/",indS),showWarnings = F,recursive = T)
  pdf(file = pPath, width = 16, height = 10)
  print(plotallchr)
  dev.off()
  
  pPath <- file.path(paste0(BaseDir,"/Introgressions/",indS,".allchr.events.jpeg"))
  jpeg(file = pPath, width = 1200, height = 800,quality = 100,res=100)
  print(plotallchr)
  dev.off()
 }
                      

