# header ------------------------------------------------------------------

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

# intersect MUMUmer results
# and make final table

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
baseDir <- argsVal[3]

dirRes12 <- file.path(baseDir, "MUMmer", paste(ref1Label, ref2Label, sep = "_"))
file12 <- list.files(dirRes12, pattern = "prt\\.snps$", full.names = T)

dirRes21 <- file.path(baseDir, "MUMmer", paste(ref2Label, ref1Label, sep = "_"))
file21 <- list.files(dirRes21, pattern = "prt\\.snps$", full.names = T)

fileOut <- gsub(pattern = "prt", "intersect", file12)

data12 <- read.table(file12, header = F, skip = 4)
headerData <- c("pos1", "allele1", "allele2", "pos2", "buff", "dist", "lenR", "lenQ", "frm", "strand", "chrom1", "chrom2")
colnames(data12) <- headerData

data21 <- read.table(file21, header = F, skip = 4)
colnames(data21) <- headerData

# correct data 21
# index of inversion 
indInv <- which(data21$strand == -1)

# esempio stringhe in una inversione
# data ref1 vs ref2
# 500676	T	C	521101	3	202314	707288	723414	1	-1	chrX	chrX
# data ref2 vs ref1
# 521101	G	A	500676	3	202314	723414	707288	1	-1	chrX	chrX

# testina
# posizioni di una bella inversione in YPS138 vs YPS128 più altre 1000 per avere dei marker standard
# which(data21$pos1 == 511293 & data21$chrom1 == "chrX")
# which(data12$pos1 == 500665 & data12$chrom1 == "chrX")
# test21 <- data21[c(573609:574609), ]
# test12 <- data12[c(574388:575388), ]
# indInv <- which(test21$tag == -1)
# test21$allele2[indInv] <- chartr("ATCG", "TAGC", test21$allele2[indInv])
# test21$allele1[indInv] <- chartr("ATCG", "TAGC", test21$allele1[indInv])
# string12 <- paste(test12$chrom1, test12$pos1, test12$allele1, test12$chrom2, test12$pos2, test12$allele2, sep = "_")
# string21 <- paste(test21$chrom2, test21$pos2, test21$allele2, test21$chrom1, test21$pos1, test21$allele1, sep = "_")
# booMarkers <- is.element(string12, string21)
# intersectMarker <- test12[booMarkers, ]

data21$allele2[indInv] <- chartr("ATCG", "TAGC", data21$allele2[indInv])
data21$allele1[indInv] <- chartr("ATCG", "TAGC", data21$allele1[indInv])

# interseco anche gli alleli, così filtro un marker se cade in una inversione in 12 e in una posizione standard in 21  
string12 <- paste(data12$chrom1, data12$pos1, data12$allele1, data12$chrom2, data12$pos2, data12$allele2, sep = "_")
string21 <- paste(data21$chrom2, data21$pos2, data21$allele2, data21$chrom1, data21$pos1, data21$allele1, sep = "_")

# prendo i marker di 12
booMarkers <- is.element(string12, string21)
intersectMarker <- data12[booMarkers, ]
# number of marker reteined / number of markers from data12
cat("number of marker retained / number of markers from ", basename(file12), ": ", nrow(intersectMarker), "/", nrow(data12), "\n", sep = "")

# append to file initialized by MUMmer.sh
write.table(file = fileOut, x = intersectMarker, append = T, col.names = F, row.names = F, quote = F, sep = "\t")



