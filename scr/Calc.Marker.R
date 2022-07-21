# header ------------------------------------------------------------------

# suddivide i marker in 4 categorie:
# traslocazioni/inversioni intra/inter-cromosoma
# e salva 4 tabelle
# calcola le frazioni di sequence & structural divergence
# va lanciato dopo Clrs
# perché legge le posizioni marker filtrate

# OCCHIO
# la tabella IntraChrom.Translocation.Markers contiene anche i marker
# upstream e downstream, tranne quando la traslocazione è alla fine del cromosoma
# (quindi il marker downstream non esiste)

rm(list = ls())
options(stringsAsFactors = F)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]
# commentami
 # baseDir <- "/home/nico/testseqloratorypipeline"
 # ref1 <- "Scc"
 # ref2 <- "EU"
 # ref1Label <- "Scc"
 # ref2Label <- "CBS432"

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

markFile <- file.path(baseDir, "LOH", "Markers", paste0("Markers.", ref1, "-", ref2, ".RData"))
# carica dfGG
load(file = markFile)
# conta tutti i marker
nMarker <- nrow(dfGG)

# output statistics file
fileStat <- file.path(baseDir, "LOH", "Markers", paste0("Stat.Markers.", ref1, "-", ref2, ".txt"))

# output file for intra-chromosome translocation
fileOutIntraTransloc <- file.path(baseDir, "LOH", "Markers", paste0("IntraChrom.Translocation.Markers.", ref1, "-", ref2, ".txt"))
cat(colnames(dfGG), "\n", file = fileOutIntraTransloc, sep = "\t")

chrLenref1 <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref1Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
chrLenref2 <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref2Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
dfChrLen <- data.frame(chrLenref1, chrLenref2)

# inter-chromosome translocation e salva
indInterChromTransloc <- which(dfGG[, 1] != dfGG[, 3] & dfGG[, 7] == 1 & dfGG[, 8] == 1)
nInterChromTransloc <- length(indInterChromTransloc)
fileOutInterTransloc <- file.path(baseDir, "LOH", "Markers", paste0("InterChrom.Translocation.Markers.", ref1, "-", ref2, ".txt"))
write.table(x = dfGG[indInterChromTransloc, ], file = fileOutInterTransloc, 
            col.names = T, row.names = F, quote = F, append = F, sep = "\t")
# inter-chromosome inversion 
indInterChromInversion <- which(dfGG[, 1] != dfGG[, 3] & (dfGG[, 7] == -1 | dfGG[, 8] == -1))
nInterChromInversion <- length(indInterChromInversion)
fileOutInterInversion <- file.path(baseDir, "LOH", "Markers", paste0("InterChrom.Inversion.Markers.", ref1, "-", ref2, ".txt"))
write.table(x = dfGG[indInterChromInversion, ], file = fileOutInterInversion, 
            col.names = T, row.names = F, quote = F, append = F, sep = "\t")

# intra-chromosome inversion 
indIntraChromInversion <- which(dfGG[, 1] == dfGG[, 3] & (dfGG[, 7] == -1 | dfGG[, 8] == -1))
nIntraChromInversion <- length(indIntraChromInversion)
fileOutIntraInversion <- file.path(baseDir, "LOH", "Markers", paste0("IntraChrom.Inversion.Markers.", ref1, "-", ref2, ".txt"))
write.table(x = dfGG[indIntraChromInversion, ], file = fileOutIntraInversion, 
            col.names = T, row.names = F, quote = F, append = F, sep = "\t")

# filter out all
indOut <- c(indInterChromTransloc, indInterChromInversion, indIntraChromInversion)

nInvInterChromTrans <- length(indOut)
if (nInvInterChromTrans != 0) {
  # in dfGGFilt tutti i marker mappano sullo stesso chrom
  dfGGFilt <- dfGG[-indOut, ]
}

# by chromosome, by reference (there's NO inter-chromosome marker!)
counterChr <- 1
nMarkerIntraChromTransloc <- 0
for (indC in allChr) {
  # puoi pescare il chrom giusto sulla base di un solo reference, mitico!
  dfGGFiltChr <- dfGGFilt[which(dfGGFilt[, 1] == indC), ]
  # indici delle colonne delle posizioni nei due reference
  for (indrefCol in c(2, 4)) {
    vctPos <- dfGGFiltChr[, indrefCol]
    # add chromosome start & end
    vcfPos <- c(1, vctPos, dfChrLen[counterChr, indrefCol / 2])
    # trans1: 9,10,11
    # trans2: 19,20,21,22
    # vctPos <- c(1,3,5,8,12,34,36,38,23,26,28,44,48,54,58,100,110,115,82,88,90,92,200)
    # chiaramente becca anche le traslocazioni con un solo marker
    indDiff <- which(diff(vctPos) < 0)
    nTransloc <- length(indDiff)
    if (nTransloc != 0) {
      indStartTransloc <- indDiff + 1
      indEndTransloc <- numeric(length = nTransloc)
      for (indS in c(1:nTransloc)) {
        lastMarker <- which(vctPos > vctPos[indDiff[indS]])[1] - 1
        # se la traslocazione è alla fine del cromosoma lastMarler = NA
        if (is.na(lastMarker)) {
          indEndTransloc[indS] <- vctPos[length(vctPos)]
        } else {
          indEndTransloc[indS] <- lastMarker
        }
        dfIntraTransloc <- dfGGFiltChr[c(c(indStartTransloc[indS] - 1):c(indEndTransloc[indS] + 1)), ]
        # la tabella IntraChrom.Translocation.Markers contiene anche i marker
        # upstream e downstream, tranne quando la traslocazione è alla fine del cromosoma
        # (quindi il marker downstream non esiste)
        write.table(x = dfIntraTransloc, file = fileOutIntraTransloc, 
                    col.names = F, row.names = F, quote = F, append = T, sep = "\t")
      }
      nMarkerIntraChromTransloc <- nMarkerIntraChromTransloc + sum(indEndTransloc - indStartTransloc + 1)
    }
  }
  counterChr <- counterChr + 1
}

nStructural <- nInvInterChromTrans + nMarkerIntraChromTransloc
nSequence <- (nMarker - nStructural)

fStructural <- signif(nStructural / nMarker, digits = 5)
fSequence <- signif(nSequence / nMarker, digits = 5)

cat("Seq: ", fSequence, "\n", file = fileStat, append = F)
cat("Str: ", fStructural, "\n", file = fileStat, append = T)
cat("N Marker: ", nMarker, "\n", file = fileStat, append = T)
cat("Mean Marker Density [1/kb]:", 
    signif((1 / sum(dfChrLen$chrLenref1) + 1 / sum(dfChrLen$chrLenref2)) * nMarker * 1E3 / 2, digits = 5), 
    file = fileStat, append = T)


