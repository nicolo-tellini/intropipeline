# header ------------------------------------------------------------------

# Author: Lorenzo Tattini
# makes segments plots for all samples by chromosome
# and saves all samples data by chromosome (for heat maps)
# [salva anche gli eventi status 2 ma li filtro dopo negli Heat.R]
# creates folder "AllSegments" for output pdfs
# and SE, FL subfolders for plots
# if plotSeg is TRUE
# makes all samples plots
# questa riga non serve a nulla

# very peculiar stuff:
# length(chromStartrefB) -> 1
# length(chromStartrefA) -> 16

rm(list = ls())
options(stringsAsFactors = F)

library(ggplot2)
library(scales)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]

##ref1Label <-"Scc"
##ref2Label <- "CBS432"
##ref1 <- "Scc"
##ref2 <- "EU" 
##baseDir <-'/home/nico/hulk'

# init
# to plot or not to plot? this is the question
plotSeg <- F
# input folder
dirData <- file.path(baseDir, "int")
# main output folder
dirName <- "AllSegments"

outDirHis <- file.path(dirData, dirName, "HistoLength")
outDirData <- file.path(dirData, dirName)
if (plotSeg) {
  outDirFL <- file.path(dirData, dirName, "FL")
  outDirSE <- file.path(dirData, dirName, "SE")
  dir.create(path = outDirFL, showWarnings = F, recursive = T)
  dir.create(path = outDirSE, showWarnings = F, recursive = T)
}
dir.create(path = outDirHis, showWarnings = F, recursive = T)

summaryCsv <- file.path(outDirData, paste("AllEvents", "csv", sep = "."))
if (file.exists(summaryCsv)) file.remove(summaryCsv)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

# centromeri, (sub)telomeri, chromosomes length
chromLenrefA <- read.table(file = file.path(baseDir, "cnv", "GCdata", ref1Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
chromLenrefB <- read.table(file = file.path(baseDir, "cnv", "GCdata", ref2Label, "LenChr.txt"), header = F, sep = "\t")[, 3]

centrSrefA <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann",
                                                     paste0(ref1Label, ".centromere.txt")), header = F)[, 1])
centrErefA <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann",
                                                     paste0(ref1Label, ".centromere.txt")), header = F)[, 2])
centrSrefB <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann",
                                                     paste0(ref2Label, ".centromere.txt")), header = F)[, 1])
centrErefB <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann",
                                                     paste0(ref2Label, ".centromere.txt")), header = F)[, 2])

subTelLrefA <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 1])
subTelRrefA <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 2])
subTelLrefB <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 1])
subTelRrefB <- as.numeric(read.table(file = file.path(baseDir, "ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 2])

# switch to kb
centrSrefA <- centrSrefA / 1000
centrErefA <- centrErefA / 1000
centrSrefB <- centrSrefB / 1000
centrErefB <- centrErefB / 1000
subTelLrefB <- subTelLrefB / 1000
subTelLrefA <- subTelLrefA / 1000
subTelRrefB <- subTelRrefB / 1000
subTelRrefA <- subTelRrefA / 1000
chromLenrefB <- chromLenrefB / 1000
chromLenrefA <- chromLenrefA / 1000

# add chromosome start (for shifting coordinates)
chromStartrefA <- 0.001
chromStartrefB <- 0.001

# centromere position and shift
cPosrefA <- (centrSrefA + centrErefA) / 2
cPosrefB <- (centrSrefB + centrErefB) / 2
shiftrefBrefA <- cPosrefB - cPosrefA

# shifting refA (sub)telomers, chromosome start and end
subTelLrefA <- subTelLrefA + shiftrefBrefA
subTelRrefA <- subTelRrefA + shiftrefBrefA
chromLenrefA <- chromLenrefA + shiftrefBrefA
chromStartrefA <- chromStartrefA + shiftrefBrefA

# shift centromere
cPosrefA <- cPosrefA + shiftrefBrefA

ptmChr <- proc.time()
ptmGlob <- proc.time()

# all segments files path
allFiles <- list.files(pattern = ".Seg.RData$", path = dirData, recursive = T, full.names = T)

# init variable to accumulate events @ create chr*Events
allEvVar <- c()

# write header string
headerString <- c("chr,start,first,last,end,status,len,denES,evrES,distES,denLF,evrLF,distLF,strain,hybrid,color,yVal")
write.table(x = headerString,
            file = summaryCsv, append = F, sep = "", row.names = F, col.names = F, quote = F)

# plotting ----------------------------------------------------------------

# carica i dati di tutti i campioni, un solo chr
chrCount <- 1
for (indC in allChr) {
  allEv <- c()
  myChrFiles <- grep(pattern = paste(".", indC, ".", sep = ""), x = allFiles, fixed = T, value = T)
  # if control sample has one whole-chromosome deletion e.g. chrV, all *.chrV.*.Seg.RData files do not exist
  if (length(myChrFiles) == 0) {
    chrCount <- chrCount + 1
    next()
  }
  
  # files of segments by reference, chromosome, sample in sample folder under Events
  # e.g. Events/A452R77/A452R77.chrI.refB.Seg.RData
  for (indD in myChrFiles) {
    load(indD)
    bName <- basename(indD)
    nameEl <- unlist(strsplit(x = bName, split = "\\."))
    sampID <- nameEl[1]
    refName <- nameEl[3]
    # reconstruct dataframe name
    evDF <- paste("ev", refName, sep = "")
    # get: call an R object using a character string
    nrDF <- nrow(get(evDF))
    appEv <- data.frame(get(evDF), rep(refName, nrDF), rep(sampID, nrDF))
    colnames(appEv)[14:15] <- c("strain", "hybrid")
    allEv <- rbind(allEv, appEv)
  }
  allEv$hybrid <- factor(x = allEv$hybrid, levels = unique(allEv$hybrid))
  # plot
  # add color and yVal
  nrEv <- nrow(allEv)
  color <- character(length = nrEv)
  color[which(allEv$status == 1)] <- "darkgrey"
  color[which(allEv$status == 0 & allEv$strain == ref1)] <- "blue"
  color[which(allEv$status == 2 & allEv$strain == ref1)] <- "red"
  color[which(allEv$status == 0 & allEv$strain == ref2)] <- "red"
  color[which(allEv$status == 2 & allEv$strain == ref2)] <- "blue"
  color <- factor(color)
  
  yrefA <- .25
  yrefB <- .75
  
  yVal <- numeric(length = nrEv)
  yVal[which(allEv$strain == ref1)] <- yrefA
  yVal[which(allEv$strain == ref2)] <- yrefB
  allEv <- data.frame(allEv, color, yVal)
  
  # saving data
  save(allEv, file = file.path(outDirData, paste("Events", indC, "RData", sep = ".")))
  write.table(x = allEv, file = summaryCsv, append = T, sep = ",", row.names = F, col.names = F, quote = F)
  
  # size markers
  szCentr <- 4
  szEvents <- 10
  szChromTel <- 2
  
  # switch data to kb
  allEv$start <- allEv$start / 1000
  allEv$end <- allEv$end / 1000
  allEv$first <- allEv$first / 1000
  allEv$last <- allEv$last / 1000
  
  # shifting refA data
  indrefA <- which(allEv$strain == ref1)
  allEv[indrefA, c(2:5)] <- allEv[indrefA, c(2:5)] + shiftrefBrefA[chrCount]
  
  # plot first-last marker
  if (plotSeg) {
    pAllSeg <- ggplot(allEv) +
      # title
      ggtitle(indC) +
      coord_cartesian(ylim = c(.0, 1.)) +
      scale_y_continuous(breaks = allEv$yVal, labels = allEv$strain) +
      # events
      geom_segment(aes(x = allEv$first, y = allEv$yVal, xend = allEv$last, yend = allEv$yVal),
                   colour = allEv$color, size = szEvents) +
      # whole chromosome
      geom_segment(aes(x = chromStartrefB, y = yrefB, xend = chromLenrefB[chrCount], yend = yrefB),
                   colour = "black", size =  szChromTel) +
      geom_segment(aes(x = chromStartrefA[chrCount], y = yrefA, xend = chromLenrefA[chrCount], yend = yrefA),
                   colour = "black", size =  szChromTel) +
      # (sub)telomers
      geom_segment(aes(x = chromStartrefB, y = yrefB, xend = subTelLrefB[chrCount], yend = yrefB),
                   colour = "orange", size =  szChromTel) +
      geom_segment(aes(x = subTelRrefB[chrCount], y = yrefB, xend = chromLenrefB[chrCount], yend = yrefB),
                   colour = "orange", size =  szChromTel) +
      geom_segment(aes(x = chromStartrefA[chrCount], y = yrefA, xend = subTelLrefA[chrCount], yend = yrefA),
                   colour = "orange", size =  szChromTel) +
      geom_segment(aes(x = subTelRrefA[chrCount], y = yrefA, xend = chromLenrefA[chrCount], yend = yrefA),
                   colour = "orange", size =  szChromTel) +
      # centromere
      geom_point(shape = 1, aes(x = cPosrefB[chrCount], y = yrefB), size = szCentr) +
      geom_point(shape = 1, aes(x = cPosrefA[chrCount], y = yrefA), size = szCentr) +
      # plain background
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            # title style
            plot.title = element_text(size = 16, face = "bold", hjust = 1),
            plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
            # axis style
            axis.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 16),
            # subtitle format
            strip.background = element_rect(fill = "white", colour = "white"),
            strip.text.x = element_text(size = 12, face = "bold"),
            plot.margin = unit(c(1, 1, 1, 1), "cm")) +
      # position shifted to ref2
      xlab(paste0("Position [kb]\n(", ref2, ")")) +
      ylab(NULL) +
      facet_wrap(facets = ~ hybrid, ncol = 4, scales = "free_x") +
      # integer x-axis values
      scale_x_continuous(labels = comma)
    
    pPath <- file.path(outDirFL, paste("AllSeg", indC, "pdf", sep = "."))
    pdf(file = pPath, width = 26, height = 10)
    print(pAllSeg)
    dev.off()
  }
  
  # plot start-end event
  if (plotSeg) {
    pAllSeg <- ggplot(allEv) +
      # title
      ggtitle(indC) +
      coord_cartesian(ylim = c(.0, 1.)) +
      scale_y_continuous(breaks = allEv$yVal, labels = allEv$strain) +
      # events
      geom_segment(aes(x = allEv$start, y = allEv$yVal, xend = allEv$end, yend = allEv$yVal),
                   colour = allEv$color, size = szEvents) +
      # whole chromosome
      geom_segment(aes(x = chromStartrefB, y = yrefB, xend = chromLenrefB[chrCount], yend = yrefB),
                   colour = "black", size =  szChromTel) +
      geom_segment(aes(x = chromStartrefA[chrCount], y = yrefA, xend = chromLenrefA[chrCount], yend = yrefA),
                   colour = "black", size =  szChromTel) +
      # (sub)telomers
      geom_segment(aes(x = chromStartrefB, y = yrefB, xend = subTelLrefB[chrCount], yend = yrefB),
                   colour = "orange", size =  szChromTel) +
      geom_segment(aes(x = subTelRrefB[chrCount], y = yrefB, xend = chromLenrefB[chrCount], yend = yrefB),
                   colour = "orange", size =  szChromTel) +
      geom_segment(aes(x = chromStartrefA[chrCount], y = yrefA, xend = subTelLrefA[chrCount], yend = yrefA),
                   colour = "orange", size =  szChromTel) +
      geom_segment(aes(x = subTelRrefA[chrCount], y = yrefA, xend = chromLenrefA[chrCount], yend = yrefA),
                   colour = "orange", size =  szChromTel) +
      # centromere
      geom_point(shape = 1, aes(x = cPosrefB[chrCount], y = yrefB), size = szCentr) +
      geom_point(shape = 1, aes(x = cPosrefA[chrCount], y = yrefA), size = szCentr) +
      # plain background
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            # title style
            plot.title = element_text(size = 16, face="bold", hjust = 1),
            plot.subtitle = element_text(size = 12, face="bold", hjust = 1),
            # axis style
            axis.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 16),
            # subtitle format
            strip.background = element_rect(fill = "white", colour = "white"),
            strip.text.x = element_text(size = 12, face = "bold"),
            plot.margin = unit(c(1, 1, 1, 1), "cm")) +
      # position shifted to ref2
      xlab(paste0("Position [kb]\n(", ref2, ")")) +
      ylab(NULL) +
      facet_wrap(facets = ~ hybrid, ncol = 4, scales = "free_x") +
      # integer x-axis values
      scale_x_continuous(labels = comma)
    
    pPath <- file.path(outDirSE, paste("AllSeg", indC, "pdf", sep = "."))
    pdf(file = pPath, width = 26, height = 10)
    print(pAllSeg)
    dev.off()
  }
  
  chrCount <- chrCount + 1
  
  # create chr*Events
  evVar <- paste(indC, "Events", sep = "")
  assign(evVar, value = allEv)
  allEvVar <- c(allEvVar, evVar)
} # chromosome loop
proc.time() - ptmChr

# prepare data for events density plots -----------------------------------
concEvents <- paste("globEvents <- rbind.data.frame(", paste(allEvVar, collapse = ", "), ")", sep = "")
# eval -> execute and create globEvents variable
eval(parse(text = concEvents))

# convertion to kb
globEvents$distES <- globEvents$distES / 1000
globEvents$distLF <- globEvents$distLF / 1000

# ref1 plots ----------------------------------------------------------------

plotEvents <- globEvents[which(globEvents$strain == ref1 & globEvents$status == 0), ]

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) +
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain"),
          subtitle = paste("Length < ", smallThres, " kb; ",
                           "N = ", nEv, "; ",
                           "Bin = ", binWsmall, " kb",
                           sep = "")) +
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall,
                 colour = "black", fill = "white") +
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") +
  ylab("Density") +
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1),
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
    # axis style
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.small.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()
pPath <- file.path(outDirHis, paste0(ref1, ".LF.small.jpg"))
jpeg(file = pPath)
print(pLFs)
dev.off()

# rifare filtri sulla distES del dataframe, ora sono solo su indSmall <- which(plotEvents$distLF < smallThres)
# plot histogram with kernel density: small ES
# pESs <- ggplot(plotSmallEvents, aes(x = distES)) +
#   ggtitle("Start-End Position Distance, refA strain",
#           subtitle = paste("Length < ", smallThres, " kb; ",
#                            "N = ", nEv, "; ",
#                            "Bin = ", binWsmall, " kb",
#                            sep = "")) +
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWsmall,
#                  colour = "black", fill = "white") +
#   geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") +
#   ylab("Density") +
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1),
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"),
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "refA.ES.small.pdf")
# pdf(file = pPath)
# print(pESs)
# dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) +
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain"),
          # unicode character â‰¥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ",
                           "N = ", nEv, "; ",
                           "Bin = ", binWlarge, " kb",
                           sep = "")) +
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge,
                 colour = "black", fill = "white") +
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") +
  ylab("Density") +
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1),
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
    # axis style
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.large.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()
pPath <- file.path(outDirHis, paste0(ref1, ".LF.large.jpg"))
jpeg(file = pPath)
print(pFLl)
dev.off()

# rifare filtri sulla distES del dataframe, ora sono solo su indSmall <- which(plotEvents$distLF < smallThres)
# plot histogram with kernel density: large ES
# pESl <- ggplot(plotLargeEvents, aes(x = distES)) +
#   ggtitle("Start-End Position Distance, refA strain",
#           # unicode character â‰¥ \u2265
#           subtitle = paste("Length > ", smallThres, " kb; ",
#                            "N = ", nEv, "; ",
#                            "Bin = ", binWlarge, " kb",
#                            sep = "")) +
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWlarge,
#                  colour = "black", fill = "white") +
#   geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") +
#   ylab("Density") +
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1),
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"),
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "refA.ES.large.pdf")
# pdf(file = pPath)
# print(pESl)
# dev.off()

# ref2 plots ----------------------------------------------------------------

plotEvents <- globEvents[which(globEvents$strain == ref2 & globEvents$status == 0), ]

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) +
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain"),
          subtitle = paste("Length < ", smallThres, " kb; ",
                           "N = ", nEv, "; ",
                           "Bin = ", binWsmall, " kb",
                           sep = "")) +
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall,
                 colour = "black", fill = "white") +
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") +
  ylab("Density") +
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1),
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
    # axis style
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.small.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()
pPath <- file.path(outDirHis, paste0(ref2, ".LF.small.jpg"))
jpeg(file = pPath)
print(pLFs)
dev.off()

# rifare filtri sulla distES del dataframe, ora sono solo su indSmall <- which(plotEvents$distLF < smallThres)
# plot histogram with kernel density: small ES
# pESs <- ggplot(plotSmallEvents, aes(x = distES)) +
#   ggtitle("Start-End Position Distance, refB strain",
#           subtitle = paste("Length < ", smallThres, " kb; ",
#                            "N = ", nEv, "; ",
#                            "Bin = ", binWsmall, " kb",
#                            sep = "")) +
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWsmall,
#                  colour = "black", fill = "white") +
#   geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") +
#   ylab("Density") +
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1),
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"),
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "refB.ES.small.pdf")
# pdf(file = pPath)
# print(pESs)
# dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) +
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain"),
          # unicode character â‰¥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ",
                           "N = ", nEv, "; ",
                           "Bin = ", binWlarge, " kb",
                           sep = "")) +
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge,
                 colour = "black", fill = "white") +
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") +
  ylab("Density") +
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1),
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
    # axis style
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.large.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()
pPath <- file.path(outDirHis, paste0(ref2, ".LF.large.jpg"))
jpeg(file = pPath)
print(pFLl)
dev.off()

# plot histogram with kernel density: large ES
# pESl <- ggplot(plotLargeEvents, aes(x = distES)) +
#   ggtitle("Start-End Position Distance, refB strain",
#           # unicode character â‰¥ \u2265
#           subtitle = paste("Length > ", smallThres, " kb; ",
#                            "N = ", nEv, "; ",
#                            "Bin = ", binWlarge, " kb",
#                            sep = "")) +
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWlarge,
#                  colour = "black", fill = "white") +
#   geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") +
#   ylab("Density") +
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1),
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1),
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"),
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "refB.ES.large.pdf")
# pdf(file = pPath)
# print(pESl)
# dev.off()

proc.time() - ptmGlob
