# Thu Nov 17 17:56:59 2022 

# Title:
# Author: Lorenzo Tattini
# Modified : Nicolò Tellini
# Status: Draft

# Comments: 
# use intersection of standard nucmer
# all chromosomes at once
# plot density of markers calculated by nucmer (after (sub)telomeric positions filtering)
# filters multiallelic positions here: remove non-matching GTs (i 2 non vengono contemplati)
# filters for markers: (sub)telomeric, quality, genotype, allele, deletions
# save vcf data (before sorting for LOH segment calls!!!)
# calculate segments of consecutive genotypes

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

# Variables ----

# ref1Label <- "Scc"
# ref2Label <- "CBS432"
# ref1 <-  "Scc"
# ref2 <- "EU"
# BaseDir <- "/home/ntellini/GITHUB"

argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
BaseDir <- argsVal[5]

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", 
            "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

polDir <- file.path(BaseDir, "mrk")
cnvDir <- file.path(BaseDir, "cnv")
resDir <- file.path(BaseDir, "int")
bWidth <- 1.0
outDirMUM <- file.path(resDir, "markers")


# Libraries ----

library(ggplot2)
library(scales)
library(data.table)
library(seqinr)


# function(s) ----

SpecDec <- function(x, k) as.numeric(format(round(x, k), nsmall = k))

ptmGlob <- proc.time()
ptmInit <- proc.time()

# body ----

dir.create(path = outDirMUM, showWarnings = F, recursive = T)

# label tables marker density
outLabelA <- "Density.Table.FltSubTel"
# label plots
outLabelB <- "Distance.FltSubTel"

# centromeri, (sub)telomeri, chromosomes length
chrLenWE <- read.table(file = file.path(BaseDir, "cnv", "GCdata", ref1Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
chrLenNA <- read.table(file = file.path(BaseDir, "cnv", "GCdata", ref2Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
centrSWE <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref1Label, ".centromere.txt")), header = F)[, 1])
centrEWE <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref1Label, ".centromere.txt")), header = F)[, 2])
centrSNA <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref2Label, ".centromere.txt")), header = F)[, 1])
centrENA <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref2Label, ".centromere.txt")), header = F)[, 2])
subTelLWE <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 1])
subTelRWE <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 2])
subTelLNA <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 1])
subTelRNA <- as.numeric(read.table(file = file.path(BaseDir, "ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 2])

myMarkers <- paste0(BaseDir,"/rep/mrktab")

# shift
cPosWE <- (centrSWE + centrEWE) / 2
cPosNA <- (centrSNA + centrENA) / 2
shiftNAWE <- cPosNA - cPosWE

# assembly di riferimento per plot label
labAssemblyRiferimento <- ref2

# retrieve marker positions ----

stdMUM <- fread(file = myMarkers,header = T,data.table = F)
stdMUM <- stdMUM[,grep(pattern = paste0(ref1Label,"|",ref2Label), x = colnames(stdMUM))]

scc_genome <- read.fasta(file = paste0(BaseDir,"/ref/Scc.genome.fa"),set.attributes = F,forceDNAtolower = F)
sp <- grep("Scc",list.files(path = paste0(BaseDir,"/ref/"),pattern = ".fa$"),invert = T,value = T)
sp_genome <- read.fasta(file = paste0(BaseDir,"/ref/",sp),set.attributes = F,forceDNAtolower = F)

stdMUM$allele_sc <- ""
stdMUM$allele_sp <- ""
fstdMUM <- data.frame()
for (i in 1:length(allChr)) {
  print(i)
  tstdMUM <- stdMUM[stdMUM[,1] == paste0(allChr[i],"_",ref1Label),]
  tstdMUM[,"allele_sc"] <- scc_genome[[paste0(allChr[i],"_",ref1Label)]][tstdMUM[,2]]
  tstdMUM[,"allele_sp"] <- sp_genome[[paste0(allChr[i],"_",ref2Label)]][tstdMUM[,4]]
  fstdMUM <- rbind(fstdMUM,tstdMUM)
}
rm(tstdMUM,stdMUM,scc_genome,sp_genome)

dfGG <- fstdMUM
rm(fstdMUM)
colnames(dfGG) <- c(paste0("chr", ref1), paste0("pos", ref1), paste0("chr", ref2), paste0("pos", ref2), paste0("ref", ref1), paste0("ref", ref2))

dfGG[, 2] <- as.numeric(dfGG[, 2])
dfGG[, 4] <- as.numeric(dfGG[, 4])

dfGG[,1] <- sapply(strsplit(dfGG[,1] ,"_"),"[[",1)
dfGG[,3] <- sapply(strsplit(dfGG[,3] ,"_"),"[[",1)

# remove chrXIV:1:38,500 known introgression on CBS432
dfGG <- dfGG[!(dfGG$chrScc == "chrXIV" & dfGG$posScc > 1 & dfGG$posScc < 38500),]

# save dfGG for de novo variants filtering
save(dfGG, file = file.path(outDirMUM, paste0("Markers.", ref1, "-", ref2, ".RData")))

# remove MT
outMT <- unique(c(which(dfGG[, 1] == "chrMT"), which(dfGG[, 3] == "chrMT")))
if (length(outMT) > 0 ) {
  dfGG <- dfGG[-outMT, ]
}

# remove (sub)telomeric markers ----

# numero del chr del marker ref1
indChrref1 <- match(dfGG[, 1], allChr)
booWE1 <- dfGG[, 2] > subTelLWE[indChrref1]
booWE2 <- dfGG[, 2] < subTelRWE[indChrref1]
# boolean delle posiizoni dei marker nel core dei chr ref1
booWE <- booWE1 & booWE2
# numero del chr del marker ref2
indChrref2 <- match(dfGG[, 3], allChr)
booNA1 <- dfGG[, 4] > subTelLNA[indChrref2]
booNA2 <- dfGG[, 4] < subTelRNA[indChrref2]
# boolean delle posiizoni dei marker nel core dei chr ref2
booNA <- booNA1 & booNA2
# boolean of non-(sub)telomeric markers
booref1ref2 <- booWE & booNA
# markers not within (sub)telomeric positions
dfGG <- dfGG[booref1ref2, ]
# plot filtered data (pasted from Markers.Dist che fa lo stesso sui dati non filtrati)
dfSNP <- dfGG
# chr and pos c(1, 2) for WE (refAsm1): sort by chr
dfSNPrefAsm1 <- dfSNP[order(match(dfSNP[, 1], allChr)), ][, c(1,2)]
# chr and pos c(3, 4) for NA (refAsm2): sort by chr
dfSNPrefAsm2 <- dfSNP[order(match(dfSNP[, 3], allChr)), ][, c(3,4)]
# header tables marker density
fileDenWE <- file.path(outDirMUM, paste(outLabelA, ref1, "txt", sep = "."))
fileDenNA <- file.path(outDirMUM, paste(outLabelA, ref2, "txt", sep = "."))
headerDensity <- data.frame("chr", "length", "# markers", "density")
fwrite(x = headerDensity, file = fileDenWE, col.names = F, row.names = F, quote = F, append = F, sep = "\t")
fwrite(x = headerDensity, file = fileDenNA, col.names = F, row.names = F, quote = F, append = F, sep = "\t")

accChrWE <- c()
accChrNA <- c()
accDiffWE <- c()
accDiffNA <- c()
counterChr <- 1
for (indC in allChr) {
  
  # refAsm1 pos sorted
  posGG <- sort(dfSNPrefAsm1[which(dfSNPrefAsm1[, 1] == indC), 2])
  lenPos <- length(posGG)
  diffPosWE <- diff(posGG)
  accChrWE <- c(accChrWE, rep(indC, lenPos - 1))
  accDiffWE <- c(accDiffWE, diffPosWE)
  denWE <- lenPos / chrLenWE[counterChr]
  fwrite(x = data.frame(indC, chrLenWE[counterChr], lenPos, denWE),  
         file = fileDenWE, 
         col.names = F, row.names = F, quote = F, append = T, sep = "\t")
  
  # refAsm2 pos sorted
  posGG <- sort(dfSNPrefAsm2[which(dfSNPrefAsm2[, 1] == indC), 2])
  lenPos <- length(posGG)
  diffPosNA <- diff(posGG)
  accChrNA <- c(accChrNA, rep(indC, lenPos - 1))
  accDiffNA <- c(accDiffNA, diffPosNA)
  denNA <- lenPos / chrLenNA[counterChr]
  write.table(x = data.frame(indC, chrLenNA[counterChr], lenPos, denNA), 
              file = fileDenNA, 
              col.names = F, row.names = F, quote = F, append = T, sep = "\t")
  
  counterChr <- counterChr + 1
}

# refAsm1 dataframe for ggplot 
dfGGDist4WE <- data.frame(accChrWE, accDiffWE)
colnames(dfGGDist4WE) <- c("chr", "dist")
# refAsm2 dataframe for ggplot 
dfGGDist4NA <- data.frame(accChrNA, accDiffNA)
colnames(dfGGDist4NA) <- c("chr", "dist")

# make factor to set facet panels order
dfGGDist4WE$chr <- factor(dfGGDist4WE$chr, levels = allChr)
dfGGDist4NA$chr <- factor(dfGGDist4NA$chr, levels = allChr)
dfGGDist4WE <- dfGGDist4WE[order(match(dfGGDist4WE$chr, allChr)), ]
dfGGDist4NA <- dfGGDist4NA[order(match(dfGGDist4NA$chr, allChr)), ]

# refAsm1 markers distant < 0.5 kb
dfGGDist4WE500 <- dfGGDist4WE[which(dfGGDist4WE$dist < 500), ]
bWidth <- 10
nEv <- nrow(dfGGDist4WE500)

# refAsm2 markers distant < 0.5 kb
dfGGDist4NA500 <- dfGGDist4NA[which(dfGGDist4NA$dist < 500), ]
bWidth <- 10
nEv <- nrow(dfGGDist4NA500)

# sample loop ----

# per tutti gli ibridi
lFiles <- basename(list.files(path = polDir, recursive = T, pattern = "vcf\\.RData$"))
sampIDs <- unique(sapply(strsplit(lFiles, split = ".", fixed = T), "[[", 1))

# load already runned samples

# load table ploidy
# sample for which either the ploidy is not available or is greater than 2 must skip the CNVs markers esclusion 
# tabinfostrains <- read.table(paste0(BaseDir,"/seq/ploidy"),header = F)
# tabinfostrains[is.na(tabinfostrains[,2]),2] <- 9999

for (ind1 in sampIDs) {
  
  # make path to output folder
  sampOutDir <- file.path(resDir, ind1)
  dir.create(sampOutDir, showWarnings = F, recursive = T)
  
  # Quality filters ---------------------------------------------------------
  
  # load WE vcf
  pathSampWE <- file.path(polDir, ind1, grep(ref1Label, 
                                             grep(pattern = paste0("^",ind1,".",ref1Label), 
                                                  lFiles, value = T), value = T))
  load(pathSampWE) # load vcfData list
  vcfFiltWE <- list(vcfData$meta, as.data.frame(vcfData$fix), as.data.frame(vcfData$format))
  names(vcfFiltWE) <- c("meta", "fix", "format")
  
  # QUAL calculation
  vcfFiltWE$fix[, 6] <- as.numeric(vcfFiltWE$fix[, 6])
  mWE <- 20
  sWE <- 0
  cat(paste(ind1, " ", ref1, " mean(QUAL): ", mWE, "\n", 
            ind1, " ", ref1, " sd(QUAL): ", sWE, "\n", sep = ""), 
      file = file.path(sampOutDir, paste(ind1, ref1, "QUALdata.txt", sep = ".")))
  indQualWE <- which(vcfFiltWE$fix[, 6] > c(mWE - sWE))
  
  # save low quality vcf
  nCoreVarWE <- length(vcfFiltWE$fix$POS)
  indLowQualWE <- setdiff(c(1:nCoreVarWE), indQualWE)
  lowQualWEvcf <- list(vcfFiltWE$meta, vcfFiltWE$fix[indLowQualWE, ], vcfFiltWE$format[indLowQualWE, ])
  names(lowQualWEvcf) <- c("meta", "fix", "format")
  pathOutLowQualWE <- file.path(sampOutDir, paste(ind1, "LQ", ref1, "RData", sep = "."))
  cat(c("Fraction of", ref1, "low quality variants (/# core positions)", paste(ind1, ":", sep = ""), 
        paste(length(indLowQualWE), nCoreVarWE, sep = "/"), "\n"))
  save(list = c("lowQualWEvcf"), file = pathOutLowQualWE)
  
  # QUAL filtering
  vcfFiltWE$fix <- vcfFiltWE$fix[indQualWE, ]
  vcfFiltWE$format <- vcfFiltWE$format[indQualWE, ]
  
  # load NA vcf
  pathSampNA <- file.path(polDir, ind1, grep(ref2Label, 
                                             grep(pattern = paste0("^",ind1,".",ref2Label), 
                                                  lFiles, value = T), value = T))
  load(pathSampNA) # load vcfData list
  vcfFiltNA <- list(vcfData$meta, as.data.frame(vcfData$fix), as.data.frame(vcfData$format))
  names(vcfFiltNA) <- c("meta", "fix", "format")
  
  # QUAL calculation
  vcfFiltNA$fix[, 6] <- as.numeric(vcfFiltNA$fix[, 6])
  mNA <- 20
  sNA <- 0
  cat(paste(ind1, " ", ref2, " mean(QUAL): ", mNA, "\n", 
            ind1, " ", ref2, " sd(QUAL): ", sNA, "\n", sep = ""), 
      file = file.path(sampOutDir, paste(ind1, ref2, "QUALdata.txt", sep = ".")))
  indQualNA <- which(vcfFiltNA$fix[, 6] > c(mNA - sNA))
  
  # save low quality vcf
  nCoreVarNA <- length(vcfFiltNA$fix$POS)
  indLowQualNA <- setdiff(c(1:nCoreVarNA), indQualNA)
  lowQualNAvcf <- list(vcfFiltNA$meta, vcfFiltNA$fix[indLowQualNA, ], vcfFiltNA$format[indLowQualNA, ])
  names(lowQualNAvcf) <- c("meta", "fix", "format")
  pathOutLowQualNA <- file.path(sampOutDir, paste(ind1, "LQ", ref2, "RData", sep = "."))
  cat(c("Fraction of", ref2, "low quality variants (/# core positions)", paste(ind1, ":", sep = ""), 
        paste(length(indLowQualNA), nCoreVarNA, sep = "/"), "\n"))
  save(list = c("lowQualNAvcf"), file = pathOutLowQualNA)
  
  # QUAL filtering
  vcfFiltNA$fix <- vcfFiltNA$fix[indQualNA, ]
  vcfFiltNA$format <- vcfFiltNA$format[indQualNA, ]
  
  # hic fuit chr loop
  # chrGG <- dfGG[, c(2, 4)]
  # bad guy
  
  chrGG <- data.frame(paste(dfGG[, 1], dfGG[, 2], sep = ""), paste(dfGG[, 3], dfGG[, 4], sep = ""), 
                      dfGG[, 5], dfGG[, 6])
  names(chrGG) <- c("posref1", "posref2", "alleleref1", "alleleref2")
  
  dfVcfWEfix <- vcfFiltWE$fix
  dfVcfNAfix <- vcfFiltNA$fix
  chrVcfWEformat <- vcfFiltWE$format[, 2]
  chrVcfNAformat <- vcfFiltNA$format[, 2]
  
  # determine common chr variants
  # booleans 
  # T se la variante è nel vcf
  # F se non c'è
  pstChrPosWE <- paste(dfVcfWEfix$CHR, dfVcfWEfix$POS, sep = "")
  pstChrPosNA <- paste(dfVcfNAfix$CHR, dfVcfNAfix$POS, sep = "")
  ieWE <- is.element(chrGG$posref1, pstChrPosWE)
  ieNA <- is.element(chrGG$posref2, pstChrPosNA)
  goldBool <- ieWE & ieNA
  # va come una scheggia
  # a <- c("29068", "29125", "29129", "29153")
  # b <- c("29086", "29125", "29129")
  # pos <- c("29064", "29068", "29086", "29125", "29129", "29153")
  # ie1 <- is.element(pos, a)
  # ie2 <- is.element(pos, b)
  # ie1 & ie2
  # mappings present in both vcfs (gold)
  # Good Guys (filtered, i.e. core mapping) Gold
  chrGGG <- chrGG[goldBool, ]
  rm(chrGG)
  # pesca le posizioni GGG nei vcf di WE e NA
  # indici del secondo vettore da prendere: dfGGG vs vcf
  # ricordati che nei non-colineari, o in caso di inversioni
  # questi vettori sono diversi (ad esempio se il marker chrI1
  # del ref1 corrisponde al marker chrX1234, indPosWE[1] sarà 1 ma
  # indPosNA sarà N >> 1)
  indPosWE <- match(chrGGG$posref1, pstChrPosWE)
  indPosNA <- match(chrGGG$posref2, pstChrPosNA)
  # in questi df le varianti sono ordinate secondo la mappa dei marker chrGGG
  dfWEfixGGG <- dfVcfWEfix[indPosWE, ]
  dfNAfixGGG <- dfVcfNAfix[indPosNA, ]
  rm(dfVcfWEfix, dfVcfNAfix)
  
  chrWEformatGGG <- chrVcfWEformat[indPosWE]
  chrNAformatGGG <- chrVcfNAformat[indPosNA]
  rm(chrVcfWEformat, chrVcfNAformat)
  
  # spacchetta il campo sample
  # convert genotype to state ID integer
  # TODO: pre-allocate data structure
  
  dfWEformatGGG <- data.frame(sapply(strsplit(chrWEformatGGG, split = ":"), "[[", 1), 
                              sapply(strsplit(chrWEformatGGG, split = ":"), "[[", 2), 
                              sapply(strsplit(chrWEformatGGG, split = ":"), "[[", 3), 
                              sapply(strsplit(chrWEformatGGG, split = ":"), "[[", 4), 
                              sapply(strsplit(chrWEformatGGG, split = ":"), "[[", 5), 
                              sapply(strsplit(chrWEformatGGG, split = ":"), "[[", 6), 
                              sapply(strsplit(chrWEformatGGG, split = ":"), "[[", 7))
  rm(chrWEformatGGG)
  dfNAformatGGG <- data.frame(sapply(strsplit(chrNAformatGGG, split = ":"), "[[", 1), 
                              sapply(strsplit(chrNAformatGGG, split = ":"), "[[", 2), 
                              sapply(strsplit(chrNAformatGGG, split = ":"), "[[", 3), 
                              sapply(strsplit(chrNAformatGGG, split = ":"), "[[", 4), 
                              sapply(strsplit(chrNAformatGGG, split = ":"), "[[", 5), 
                              sapply(strsplit(chrNAformatGGG, split = ":"), "[[", 6), 
                              sapply(strsplit(chrNAformatGGG, split = ":"), "[[", 7))
  rm(chrNAformatGGG)
  colnames(dfWEformatGGG) <- c("GT", "PL", "DP", "SP", "ADF", "ADR", "AD")
  colnames(dfNAformatGGG) <- c("GT", "PL", "DP", "SP", "ADF", "ADR", "AD")
  
  # GT filters ---------------------------------------------------------
  
  # remove non-matching GTs
  nM1 <- which(dfWEformatGGG$GT == "1/1" & dfNAformatGGG$GT != "0/0")
  nM2 <- which(dfWEformatGGG$GT == "0/0" & dfNAformatGGG$GT != "1/1")
  nM3 <- which(dfWEformatGGG$GT == "0/1" & dfNAformatGGG$GT != "0/1")
  nM4 <- which(dfWEformatGGG$GT != "0/1" & dfNAformatGGG$GT == "0/1")
  nM5 <- which(dfNAformatGGG$GT == "0/0" & dfWEformatGGG$GT != "1/1")
  nM6 <- which(dfNAformatGGG$GT == "1/1" & dfWEformatGGG$GT != "0/0")
  
  nMatchGT <- unique(c(nM1, nM2, nM3, nM4, nM5, nM6))
  
  if (length(nMatchGT) > 0) {
    cat(c("Number of non-matching GT", paste(ind1, ":", sep = ""), length(nMatchGT), "\n"))
    
    unMtWEfix <- dfWEfixGGG[nMatchGT, ]
    dfWEfixGGG <- dfWEfixGGG[-nMatchGT, ]
    unMtNAfix <- dfNAfixGGG[nMatchGT, ]
    dfNAfixGGG <- dfNAfixGGG[-nMatchGT, ]
    
    unMtWEformat <- dfWEformatGGG[nMatchGT, ]
    dfWEformatGGG <- dfWEformatGGG[-nMatchGT, ]
    unMtNAformat <- dfNAformatGGG[nMatchGT, ]
    dfNAformatGGG <- dfNAformatGGG[-nMatchGT, ]
    
    unMtWEvcf <- list(vcfFiltWE$meta, unMtWEfix, unMtWEformat)
    pathOutUMWE <- file.path(sampOutDir, paste(ind1, "UnMatchedGT", ref1, "RData", sep = "."))
    save(list = c("unMtWEvcf"), file = pathOutUMWE)
    
    unMtNAvcf <- list(vcfFiltNA$meta, unMtNAfix, unMtNAformat)
    pathOutUMNA <- file.path(sampOutDir, paste(ind1, "UnMatchedGT", ref2, "RData", sep = "."))
    save(list = c("unMtNAvcf"), file = pathOutUMNA)
    
    unMtGGG <- chrGGG[nMatchGT, ]
    chrGGG <- chrGGG[-nMatchGT, ]
    pathOutUM <- file.path(sampOutDir, paste(ind1, "UnMatchedGT.RData", sep = "."))
    save(unMtGGG, file = pathOutUM)
  }
  
  # filtering non-matching alleles: la versione V4 non funziona per strain non colineari
  # perché nei ref riarrangiati le varianti nei vcf non sono
  # ordinate come nella mappa di nucmer
  # ovvero alla riga N della mappa compaiono due varianti che non necessariamente
  # si trovano alla riga N dei due vcf
  
  # inoltre il check degli alleli non funziona
  # nelle regioni invertite
  # perché BWA mappa tutte le reverse complement
  # quindi l'allele ALT del vcf1 è il complementare
  # dell'allele REF del vcf2
  
  # trick
  # V1 <- c("A","T","T","A","G")
  # V2 <- c("A","C","T","A","G")
  # lP <- length(V1)
  # pV1 <- paste(V1, rep(1:lP), sep = "")
  # pV2 <- paste(V2, rep(1:lP), sep = "")
  # !is.element(pV1, pV2)
  
  # non 0/0 positions in WE
  indWEvAltAll <- which(!is.na(dfWEfixGGG$ALT))
  lPWE <- length(indWEvAltAll)
  # CHR, POS and ALT in WE
  cPosAltWE <- dfWEfixGGG[indWEvAltAll, c(1, 2, 5)]
  # CHR, POS and REF in NA
  cPosrefNA <- dfNAfixGGG[indWEvAltAll, c(1, 2, 4)]
  # it's a character trap!
  alleleBoolWE <- !is.element(paste0(cPosAltWE$ALT, rep(1:lPWE)), 
                              paste0(cPosrefNA$REF, rep(1:lPWE)))
  posNoMatchWE <- paste(cPosAltWE$CHROM[alleleBoolWE], cPosAltWE$POS[alleleBoolWE], sep = "")
  # mSuS1: indici delle varianti non 0/0 in WE il cui allele ALT non
  # matcha con quello REF nella relativa posizione di NA
  mSuS1 <- match(posNoMatchWE, chrGGG$posref1)
  # check for inversions in ref1
  complementAltref1 <- chartr("ATCG", "TAGC", dfWEfixGGG$ALT[mSuS1])
  refAlleleref2 <- dfNAfixGGG$REF[mSuS1]
  # T elements: la base complementare a quella ALT (del parent 1) matcha la base REF (del parent 2)
  booComp <- is.element(paste0(complementAltref1, rep(1:length(complementAltref1))), 
                        paste0(refAlleleref2, rep(1:length(refAlleleref2))))
  # T elements: i marker sono in una inversione (su uno qualsiasi dei reference)
  booInv <- (chrGGG[mSuS1, 5] == -1 | chrGGG[mSuS1, 6] == -1)
  # tieni i marker nelle inversioni con l'allele giusto 
  mSuS1 <- mSuS1[!(booComp & booInv)]
  
  # non 0/0 positions in NA
  indNAvAltAll <- which(!is.na(dfNAfixGGG$ALT))
  lPNA <- length(indNAvAltAll)
  # CHR, POS and ALT in NA
  cPosAltNA <- dfNAfixGGG[indNAvAltAll, c(1, 2, 5)]
  # CHR, POS and REF in WE
  cPosrefWE <- dfWEfixGGG[indNAvAltAll, c(1, 2, 4)]
  # it's a character trap!
  alleleBoolNA <- !is.element(paste0(cPosAltNA$ALT, rep(1:lPNA)), 
                              paste0(cPosrefWE$REF, rep(1:lPNA)))
  posNoMatchNA <- paste(cPosAltNA$CHROM[alleleBoolNA], cPosAltNA$POS[alleleBoolNA], sep = "")
  # mSuS2: indici delle varianti non 0/0 in NA il cui allele ALT non
  # matcha con quello REF della relativa posizione di WE
  mSuS2 <- match(posNoMatchNA, chrGGG$posref2)
  # check for inversions in ref2
  complementAltref2 <- chartr("ATCG", "TAGC", dfNAfixGGG$ALT[mSuS2])
  refAlleleref1 <- dfWEfixGGG$REF[mSuS2]
  # T elements: la base complementare a quella ALT (del parent 2) matcha la base REF (del parent 1)
  booComp <- is.element(paste0(complementAltref2, rep(1:length(complementAltref2))), 
                        paste0(refAlleleref1, rep(1:length(refAlleleref1))))
  # T elements: i marker sono in una inversione (su uno qualsiasi dei reference)
  booInv <- (chrGGG[mSuS2, 5] == -1 | chrGGG[mSuS2, 6] == -1)
  # tieni i marker nelle inversioni con l'allele giusto 
  mSuS2 <- mSuS2[!(booComp & booInv)]
  
  # merge non-matching allele position 
  # and remove lines from chrGGG map, vcf fix and vcf format data
  nMatchAL <- unique(c(mSuS1, mSuS2))
  
  # Alleles filters ---------------------------------------------------------
  
  if (length(nMatchAL) > 0) {
    
    cat(c("Number of non-matching alleles", paste(ind1, ":", sep = ""), length(nMatchAL), "\n"))
    
    unMtWEfix <- dfWEfixGGG[nMatchAL, ]
    dfWEfixGGG <- dfWEfixGGG[-nMatchAL, ]
    unMtNAfix <- dfNAfixGGG[nMatchAL, ]
    dfNAfixGGG <- dfNAfixGGG[-nMatchAL, ]
    
    unMtWEformat <- dfWEformatGGG[nMatchAL, ]
    dfWEformatGGG <- dfWEformatGGG[-nMatchAL, ]
    unMtNAformat <- dfNAformatGGG[nMatchAL, ]
    dfNAformatGGG <- dfNAformatGGG[-nMatchAL, ]
    
    unMtWEvcf <- list(vcfFiltWE$meta, unMtWEfix, unMtWEformat)
    pathOutUMWE <- file.path(sampOutDir, paste(ind1, "UnMatchedAL", ref1, "RData", sep = "."))
    save(list = c("unMtWEvcf"), file = pathOutUMWE)
    
    unMtNAvcf <- list(vcfFiltNA$meta, unMtNAfix, unMtNAformat)
    pathOutUMNA <- file.path(sampOutDir, paste(ind1, "UnMatchedAL", ref2, "RData", sep = "."))
    save(list = c("unMtNAvcf"), file = pathOutUMNA)
    
    unMtGGG <- chrGGG[nMatchAL, ]
    chrGGG <- chrGGG[-nMatchAL, ]
    pathOutUM <- file.path(sampOutDir, paste(ind1, "UnMatchedAL.RData", sep = "."))
    save(unMtGGG, file = pathOutUM)
  }
  
  # CNV filters ---------------------------------------------------------
  
  # load data
  inFileWE <- list.files(path = file.path(cnvDir, "results", ref1Label, ind1), 
                         pattern = "_CNVs.p.value.txt$", full.names = T)
  dfCopyNumberWE <- read.table(file = inFileWE, header = T, sep = "\t")
  inFileNA <- list.files(path = file.path(cnvDir, "results", ref2Label, ind1), 
                         pattern = "_CNVs.p.value.txt$", full.names = T)
  dfCopyNumberNA<- read.table(file = inFileNA, header = T, sep = "\t")
  
  # format chr
  dfCopyNumberWE$chr <- gsub(pattern = "_.*", replacement = "", x = dfCopyNumberWE$chr)
  dfCopyNumberNA$chr <- gsub(pattern = "_.*", replacement = "", x = dfCopyNumberNA$chr)
  
  # keep only robust gain events
  dfCopyNumberWE <- dfCopyNumberWE[which( dfCopyNumberWE$WilcoxonRankSumTestPvalue < 0.01 & 
                                            dfCopyNumberWE$KolmogorovSmirnovPvalue < 0.01), ]
  dfCopyNumberNA <- dfCopyNumberNA[which( dfCopyNumberNA$WilcoxonRankSumTestPvalue < 0.01 & 
                                            dfCopyNumberNA$KolmogorovSmirnovPvalue < 0.01), ]
  
  # if events were detected against both references
  if ( nrow(dfCopyNumberWE) != 0 & nrow(dfCopyNumberNA) != 0) {
    
    # change chromosome encoding in CNV calls
    dfCopyNumberWE$chr <- paste("chr", dfCopyNumberWE$chr, sep = "")
    dfCopyNumberNA$chr <- paste("chr", dfCopyNumberNA$chr, sep = "")
    
    # prepare markers table
    dfGGG4CNV <- data.frame(gsub(pattern = "[[:digit:]].*$", "", chrGGG$posref1), 
                            as.numeric(gsub(pattern = "^.*[[:alpha:]]", "", chrGGG$posref1)), 
                            gsub(pattern = "[[:digit:]].*$", "", chrGGG$posref2), 
                            as.numeric(gsub(pattern = "^.*[[:alpha:]]", "", chrGGG$posref2)))
    
    colnames(dfGGG4CNV) <- c("chrref1", "posref1", "chrref2", "posref2")
    
    # check which markers fall within the events
    stateCopyWE <- character(length = nrow(dfGGG4CNV))
    
    for (indCNV in 1:nrow(dfCopyNumberWE)) {
      
      indLossWE <- which(  dfGGG4CNV$chrref1 == dfCopyNumberWE$chr[indCNV] & 
                             dfGGG4CNV$posref1 >= dfCopyNumberWE$start[indCNV] & 
                             dfGGG4CNV$posref1 <= dfCopyNumberWE$end[indCNV])
      stateCopyWE[indLossWE] <- "remove"
    }
    
    stateCopyNA <- character(length = nrow(dfGGG4CNV))
    for (indCNV in 1:nrow(dfCopyNumberNA)) {
      indLossNA <- which(dfGGG4CNV$chrref2 == dfCopyNumberNA$chr[indCNV] & 
                           dfGGG4CNV$posref2 >= dfCopyNumberNA$start[indCNV] & 
                           dfGGG4CNV$posref2 <= dfCopyNumberNA$end[indCNV])
      stateCopyNA[indLossNA] <- "remove"
    }
    
    # index of markers in the events in one or both references
    
    indLoss <- which(stateCopyWE == "remove" & stateCopyNA == "remove")
    
    if (length(indLoss) > 0) {
      
      pChrPosVcfWE <- paste(dfWEfixGGG$CHROM, dfWEfixGGG$POS, sep = "")
      booWE <- is.element(pChrPosVcfWE, chrGGG$posref1[indLoss])
      inLossWEfixGGG <- dfWEfixGGG[booWE, ]
      dfWEfixGGG <- dfWEfixGGG[!booWE, ]
      inLossWEformatGGG <- dfWEformatGGG[booWE, ]
      dfWEformatGGG <- dfWEformatGGG[!booWE, ]
      
      inLossWEvcf <- list(vcfFiltWE$meta, inLossWEfixGGG, inLossWEformatGGG)
      pathOutLossWE <- file.path(sampOutDir, paste(ind1, "Loss", ref1, "RData", sep = "."))
      names(inLossWEvcf) <- c("meta", "fix", "format")
      save(list = c("inLossWEvcf"), file = pathOutLossWE)
      
      pChrPosVcfNA <- paste(dfNAfixGGG$CHROM, dfNAfixGGG$POS, sep = "")
      booNA <- is.element(pChrPosVcfNA, chrGGG$posref2[indLoss])
      inLossNAfixGGG <- dfNAfixGGG[booNA, ]
      dfNAfixGGG <- dfNAfixGGG[!booNA, ]
      inLossNAformatGGG <- dfNAformatGGG[booNA, ]
      dfNAformatGGG <- dfNAformatGGG[!booNA, ]
      
      inLossNAvcf <- list(vcfFiltNA$meta, inLossNAfixGGG, inLossNAformatGGG)
      pathOutLossNA <- file.path(sampOutDir, paste(ind1, "Loss", ref2, "RData", sep = "."))
      names(inLossNAvcf) <- c("meta", "fix", "format")
      save(list = c("inLossNAvcf"), file = pathOutLossNA)
      
      chrGGG <- chrGGG[-indLoss, ]
    }
  }
  
  # saving data ---------------------------------------------------------
  
  # save chrGGG
  
  # save filtered vcf used for segment detection
  lstWEvcfGGG <- list(vcfFiltWE$meta, dfWEfixGGG, dfWEformatGGG)
  names(lstWEvcfGGG) <- c("meta", "fix", "format")
  pathOutWEvcfGGG <- file.path(sampOutDir, paste(ind1, ref1 ,"VcfFilt.RData", sep = "."))
  save(list = c("lstWEvcfGGG"), file = pathOutWEvcfGGG)
  
  lstNAvcfGGG <- list(vcfFiltNA$meta, dfNAfixGGG, dfNAformatGGG)
  names(lstNAvcfGGG) <- c("meta", "fix", "format")
  pathOutNAvcfGGG <- file.path(sampOutDir, paste(ind1, ref2, "VcfFilt.RData", sep = "."))
  save(list = c("lstNAvcfGGG"), file = pathOutNAvcfGGG)
  
  # segment detection ---------------------------------------------------------
  
  # sort fix & format dataframes by chromosome & position
  sortWEfix <- dfWEfixGGG[order(as.numeric(dfWEfixGGG$POS)), ]
  sortWEformat <- dfWEformatGGG[order(as.numeric(dfWEfixGGG$POS)), ]
  dfWEfixGGG <- sortWEfix[order(match(sortWEfix$CHROM, allChr)), ]
  dfWEformatGGG <- sortWEformat[order(match(sortWEfix$CHROM, allChr)), ]
  
  sortNAfix <- dfNAfixGGG[order(as.numeric(dfNAfixGGG$POS)), ]
  sortNAformat <- dfNAformatGGG[order(as.numeric(dfNAfixGGG$POS)), ]
  dfNAfixGGG <- sortNAfix[order(match(sortNAfix$CHROM, allChr)), ]
  dfNAformatGGG <- sortNAformat[order(match(sortNAfix$CHROM, allChr)), ]
  
  lGen <- length(dfWEformatGGG$GT)
  # apply
  # initialized with 0
  intGenWE <- integer(length = lGen)
  intGenNA <- integer(length = lGen)
  
  for (indG in 1:lGen) {
    intGenWE[indG] <- switch(dfWEformatGGG$GT[indG], 
                             "0/0" = 0, "0/1" = 1, "1/1" = 2)
    intGenNA[indG] <- switch(dfNAformatGGG$GT[indG],
                             "0/0" = 0, "0/1" = 1, "1/1" = 2)
  }
  
  for (indC in allChr) {
    # extract chromosome length
    chrStrWE <- unlist(strsplit(grep(paste(indC, "_", sep = ""), lstWEvcfGGG$meta, value = T), split = ","))[2]
    chromLenWE <- as.numeric(gsub(pattern = "[[:alpha:][:punct:]]", x = chrStrWE, replacement = ""))
    chrStrNA <- unlist(strsplit(grep(paste(indC, "_", sep = ""), lstNAvcfGGG$meta, value = T), split = ","))[2]
    chromLenNA <- as.numeric(gsub(pattern = "[[:alpha:][:punct:]]", x = chrStrNA, replacement = ""))
    
    indVarChrWE <- which(dfWEfixGGG$CHROM == indC)
    intGenChrWE <- intGenWE[indVarChrWE]
    chrWEfixGGG <- dfWEfixGGG[indVarChrWE, ]
    
    indVarChrNA <- which(dfNAfixGGG$CHROM == indC)
    intGenChrNA <- intGenNA[indVarChrNA]
    chrNAfixGGG <- dfNAfixGGG[indVarChrNA, ]
    
    resWE <- rle(intGenChrWE)
    resNA <- rle(intGenChrNA)
    
    nMarkerWE <- sum(resWE$lengths)
    nMarkerNA <- sum(resNA$lengths)
    if (nMarkerWE == 0) {
      cat(paste("No marker", ref1, "found in chromosome", indC, "\n"))
      next()
    }
    if (nMarkerNA == 0) {
      cat(paste("No marker", ref2, "found in chromosome", indC, "\n"))
      next()
    }
    # Het fraction WE
    nHEt <- sum(resWE$lengths[which(resWE$values == 1)])
    ratHH <- nHEt / nMarkerWE
    cat(paste("Het fraction", ind1, indC, paste0(ref1, ":")), SpecDec(ratHH, 5), "\n")
    # mrkRatWE <- nMarkerWE / chrLenWE[counterChr] / vctDenWE[counterChr]
    # cat("Marker ratio", ind1, indC, paste0(ref1, ":"), mrkRatWE, nMarkerWE, length(intGenChrWE), "\n")
    # Het fraction NA
    nHEt <- sum(resNA$lengths[which(resNA$values == 1)])
    ratHH <- nHEt / nMarkerNA
    cat(paste("Het fraction", ind1, indC, paste0(ref2, ":")), SpecDec(ratHH, 5), "\n")
    # mrkRatNA <- nMarkerNA / chrLenNA[counterChr] / vctDenNA[counterChr]
    # cat("Marker ratio", ind1, indC, paste0(ref2, ":"), mrkRatNA, nMarkerNA, length(intGenChrNA), "\n")
    # events map
    lEventsWE <- length(resWE$length)
    lEventsNA <- length(resNA$length)
    
    # make events WE
    evWE <- data.frame(
      chr = rep(indC, lEventsWE), 
      start = numeric(length = lEventsWE), 
      first = numeric(length = lEventsWE), 
      last = numeric(length = lEventsWE), 
      end = numeric(length = lEventsWE), 
      status = resWE$values, 
      len = as.numeric(resWE$lengths), # number of SNPs
      denES = numeric(length = lEventsWE), 
      evrES = numeric(length = lEventsWE), 
      distES = numeric(length = lEventsWE), 
      denLF = numeric(length = lEventsWE), 
      evrLF = numeric(length = lEventsWE), 
      distLF = numeric(length = lEventsWE))
    evWE$start[1] <- 1
    evWE$first[1] <- as.numeric(chrWEfixGGG$POS[1])
    evWE$end[lEventsWE] <- chromLenWE
    evWE$last[lEventsWE] <- as.numeric(chrWEfixGGG$POS[nrow(chrWEfixGGG)])
    
    if (lEventsWE > 1) {
      # indE: index of events
      for (indE in 2:lEventsWE) {
        # ultimo snp dello stato precedente
        evWE$last[indE - 1] <- as.numeric(chrWEfixGGG$POS[sum(resWE$lengths[1:c(indE - 1)])])
        # primo snp dello stato corrente
        evWE$first[indE] <- as.numeric(chrWEfixGGG$POS[sum(resWE$lengths[1:c(indE - 1)]) + 1])
        # end of former state
        evWE$end[indE - 1] <- floor((evWE$first[indE] + evWE$last[indE - 1]) / 2) - 1
        # start of current state
        evWE$start[indE] <- floor((evWE$first[indE] + evWE$last[indE - 1]) / 2)
      }
    }
    # N snps / bp
    evWE$distES <- c(evWE$end - evWE$start)
    evWE$denES <- evWE$len / evWE$distES
    # +1 to take into account events with 1 SNP
    evWE$distLF <- c(evWE$last - evWE$first + 1)
    evWE$denLF <- evWE$len / evWE$distLF
    # bp / N snps
    evWE$evrES <- evWE$distES / evWE$len
    evWE$evrLF <- evWE$distLF / evWE$len
    # format events table
    evWE$denES <- as.numeric(format(evWE$denES, digits = 3))
    evWE$evrES <- round(evWE$evrES)
    evWE$denLF <- as.numeric(format(evWE$denLF, digits = 3))
    evWE$evrLF <- round(evWE$evrLF)
    
    # make events NA
    evNA <- data.frame(
      chr = rep(indC, lEventsNA), 
      start = numeric(length = lEventsNA), 
      first = numeric(length = lEventsNA), 
      last = numeric(length = lEventsNA), 
      end = numeric(length = lEventsNA), 
      status = resNA$values, 
      len = as.numeric(resNA$lengths), # number of SNPs
      denES = numeric(length = lEventsNA), 
      evrES = numeric(length = lEventsNA), 
      distES = numeric(length = lEventsNA), 
      denLF = numeric(length = lEventsNA), 
      evrLF = numeric(length = lEventsNA), 
      distLF = numeric(length = lEventsNA))
    evNA$start[1] <- 1
    evNA$first[1] <- as.numeric(chrNAfixGGG$POS[1])
    evNA$end[lEventsNA] <- chromLenNA
    evNA$last[lEventsNA] <- as.numeric(chrNAfixGGG$POS[nrow(chrNAfixGGG)])
    
    if (lEventsNA > 1) {
      # indE: index of events
      for (indE in 2:lEventsNA) {
        # ultimo snp dello stato precedente
        evNA$last[indE - 1] <- as.numeric(chrNAfixGGG$POS[sum(resNA$lengths[1:c(indE - 1)])])
        # primo snp dello stato corrente
        evNA$first[indE] <- as.numeric(chrNAfixGGG$POS[sum(resNA$lengths[1:c(indE -1)]) + 1])
        # end of former state
        evNA$end[indE - 1] <- floor((evNA$first[indE] + evNA$last[indE - 1]) / 2) - 1
        # start of current state
        evNA$start[indE] <- floor((evNA$first[indE] + evNA$last[indE - 1]) / 2)
      }
    }
    # N snps / bp
    evNA$distES <- c(evNA$end - evNA$start)
    evNA$denES <- evNA$len / evNA$distES
    # +1 to take into account events with 1 SNP
    evNA$distLF <- c(evNA$last - evNA$first + 1)
    evNA$denLF <- evNA$len / evNA$distLF
    # bp / N snps
    evNA$evrES <- evNA$distES / evNA$len
    evNA$evrLF <- evNA$distLF / evNA$len
    # format events table
    evNA$denES <- as.numeric(format(evNA$denES, digits = 3))
    evNA$evrES <- round(evNA$evrES)
    evNA$denLF <- as.numeric(format(evNA$denLF, digits = 3))
    evNA$evrLF <- round(evNA$evrLF)
    
    # save data for AllSeg.R
    pathOutWE <- file.path(sampOutDir, paste(ind1, indC, ref1, "Seg.RData", sep = "."))
    varOutName1 <- paste0("ev", ref1)
    assign(varOutName1, evWE)
    save(list = varOutName1, file = pathOutWE)
    pathOutNA <- file.path(sampOutDir, paste(ind1, indC, ref2, "Seg.RData", sep = "."))
    varOutName2 <- paste0("ev", ref2)
    assign(varOutName2, evNA)
    save(list = varOutName2, file = pathOutNA)
    
    # plotting segments
    indChr <- which(allChr == indC)
    yWE <- .25
    yNA <- .75
    
    # invert!!!
    nrWE <- nrow(evWE)
    color <- character(length = nrWE)
    for (indSw in 1:nrWE) {
      color[indSw] <- switch(as.character(evWE$status[indSw]), 
                             "0"="blue", 
                             # 0/1
                             "1"="darkgrey", 
                             # 1/1
                             "2"="red")
    }
    color <- factor(color)
    yVal <- rep(yWE, nrWE)
    evWE <- data.frame(evWE, color, yVal)
    
    # invert!!!
    nrNA <- nrow(evNA)
    color <- character(length = nrNA)
    for (indSw in 1:nrNA) {
      color[indSw] <- switch(as.character(evNA$status[indSw]), 
                             "0"="red", 
                             # 0/1
                             "1"="darkgrey", 
                             # 1/1
                             "2"="blue")
    }
    color <- factor(color)
    yVal <- rep(yNA, nrNA)
    evNA <- data.frame(evNA, color, yVal)
    # all events
    strain <- c(rep(ref1, nrow(evWE)), rep(ref2, nrow(evNA)))
    allEv <- rbind(evWE, evNA)
    
    # single sample single chromosome plot
    
    # size markers
    szCentr <- 8
    szEvents <- 10
    szChrom <- 2
    
    # shifting WE data
    indWE <- which(allEv$yVal == .25)
    allEv[indWE, c(2:5)] <- allEv[indWE, c(2:5)] + shiftNAWE[indChr]
    
    # plot events from first SNP to last SNP
    pGG <- ggplot(allEv) + 
      # title
      ggtitle(ind1, subtitle = indC) + 
      coord_cartesian(ylim = c(.0, 1.)) + 
      scale_y_continuous(breaks = NULL, labels = NULL) + 
      annotate("text", x = 8, y = 0.95, label = ref2, color = "red", size = theme_get()$text[["size"]] / 1.5) + 
      annotate("text", x = 8, y = 0.45, label = ref1, color = "blue", size = theme_get()$text[["size"]] / 1.5) + 
      # whole chromosome
      geom_segment(aes(x = 1, y = yNA, xend = chromLenNA, yend = yNA), 
                   colour = "black", size =  szChrom) + 
      geom_segment(aes(x = shiftNAWE[indChr], y = yWE, xend = (chromLenWE + shiftNAWE[indChr]), yend = yWE), 
                   colour = "black", size =  szChrom) + 
      # (sub)telomers
      geom_segment(aes(x = 1, y = yNA, xend = subTelLNA[indChr], yend = yNA), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = subTelRNA[indChr], y = yNA, xend = chromLenNA, yend = yNA), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = shiftNAWE[indChr], y = yWE, xend = (subTelLWE[indChr] + shiftNAWE[indChr]), yend = yWE), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = (subTelRWE[indChr] + shiftNAWE[indChr]), y = yWE, 
                       xend = (chromLenWE + shiftNAWE[indChr]), yend = yWE), 
                   colour = "orange", size =  szChrom) + 
      # events
      geom_segment(aes(x = first, y = allEv$yVal, xend = last, yend = allEv$yVal), 
                   colour = allEv$color, size = szEvents) + 
      # centromere
      geom_point(shape = 1, aes(x = (centrSNA[indChr] + centrENA[indChr]) / 2, y = yNA), size = szCentr) + 
      geom_point(shape = 1, aes(x = (centrSWE[indChr] + centrEWE[indChr]) / 2 + shiftNAWE[indChr], y = yWE), 
                 size = szCentr) + 
      # plain background
      theme(plot.margin = unit(c(1, 2, 1, 1), "cm"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(), 
            panel.background = element_blank(), 
            # title style
            plot.title = element_text(size = 16, face = "bold", hjust = 1), 
            plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
            # axis style
            axis.title = element_text(size = 16, face = "bold"), 
            axis.text = element_text(size = 16)) + 
      xlab(paste0("(", labAssemblyRiferimento, ") Position [bp]")) + 
      ylab(NULL) + 
      # integer x-axis values
      scale_x_continuous(labels = comma)
    pPath <- file.path(sampOutDir, paste(ind1, indC, "SegFL.pdf", sep = "."))
    pdf(file = pPath, width = 14, height = 4)
    print(pGG)
    dev.off()
    
  } # chromosome loop
} # sample loop

cat("Global time", "\n")
proc.time() - ptmGlob
