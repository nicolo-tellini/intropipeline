# Thu Jun  1 10:36:00 2023 

# Title:
# Author: Nicolò T.
# Status: Draft
# Input files : The data are stored in NAS
# NAS path: 

## !!! A differenza di quanto potrebbe essere riportato nelle righe di sotto il ranking viene usato SOLO per semplificare lo step di filtering e NON per la costruzione dei blocchi.
# Comments:
# Differenze con ClrS8
# - for improving performances only lists and data.table are allowed (unlsess a task can be more easily performed using a DF),
# - i DT sono modificati on place cosicchè non si generano nuovi oggetti e quando un oggetto (particolarmente grande) non serve più si rimuove,
# - carica una tabella di markers già precompilata
# - non rimuove i subtelomeri (perchè sono già esclusi nella tabella di partenza)
# - rm i markers che stanno nella regone chrXIV pos < 38500 (intro spar-->scer)
# - non fa QUAL filtering (perchè è gia fatto in fase di mrk genotyp. con bcftools) (funfact: anche senza filtrare per la QUAL i risultati non cambiano drasticamente)
# - tanti plots sono spariti, la pipeline per le intro non ne ha necessariamente bisogno (vedi AB etc...)
# - non ci sono steps che salvano i mrk scartati forse in futuro potrà essere utile reintegrarlo

# COME LAVORA:
# - upload di vcfs ridotti (senza QUAL, FILT, ID cols) (usa fread e non vcfR, vcfR è lento)
# - rimuove i mrk all'inizio di chrXIV
# - butta i mrk i cui raking non sono shared across i due vcfs (quello contro scc e quello contro CBS432)
# - per ogni ALT controlla che l'allele sia quello REF sull'altro genoma, oppure ".", se ciò non accede la posizione è scartata in ambedue i vcfs (usando come strategia quella di bool in ClrS che funziona molto bene) 
# - il GT è controllato come in ClrS8 ma, invece di usare gli index del data.frame, uso i valori rank (non cambia molto la pappa è la stessa)
# A questo punto il filtering è fatto, con fwrite salvo i mrk filtrati per entrambi scc e cbs432 (aggungendo .gz al nome del file fwrite li salva zippati)
# - posso permettermi di tenere uno solo dei 2 vcfs per il plot (in futuro questo può essere cambiato, per ora di default tiene ref1label che è sempre scc assembly)
# - qui si generano i blocchi e si annotano con le informazioni geniche e qui è la svolta, con Granges sono 3 righe di codice
# - si calcolano le solite info: num mrk, density and so on ...
# - BLUE-RED plot che alla gente piace tanto 
# - salvataggio tabella blocchi
# END

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

# Variables ----

 # ref1Label <- "Scc"
 # ref2Label <- "CBS432"
 # baseDir <- "/home/ntellini/intropipeline"

argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
baseDir <- argsVal[3]

allChr <- c("chrI", "chrII", "chrIII", "chrIV", 
            "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX","chrX", "chrXI", "chrXII", 
            "chrXIII", "chrXIV", "chrXV", "chrXVI")

mrkDir <- file.path(baseDir, "mrk")
intDir <- file.path(baseDir, "int")

dir.create(intDir)
dir.create(paste0(intDir,"/mrk/"))

# Libraries ----

library(data.table)
library(seqinr)
library(R.filesets)
library(GenomicRanges)
library(ggplot2)
library(purrr)
library(dplyr)

# functions ----

str_split_custom <- function(x,y,z) sapply(strsplit(x,split = z),"[[",y) # frequent strsplit

block_ann <- function(x){
  
  ann_l <- list()
  
  spp <- c(ref1Label)
  
  spp_load <- map(spp,function(x) paste0(baseDir,"/ref/Ann/",x,".all_feature.gff"))
  
  ann_l <- lapply(spp_load,function(x) fread(x, data.table = F,header = F) )  
  
  ann_l <- lapply(ann_l, function(x) x[x[,3] == "gene",] )
  
  gene_names <- lapply(ann_l, function(x){sub("Name=","",grep("Name", unlist(base::strsplit(x[,9],";")),value = T))} )
  
  ann_l <- mapply(cbind, ann_l, "genes"=gene_names, SIMPLIFY=F)
  
  rm(gene_names)
  
  setDF(x)
  
  names(ann_l) <- ref1Label
  
  fun_ann <- function(x){
    sp <- ref1Label
    chr <- as.character(x[4])
    st <- as.integer(x[2])
    en <- as.integer(x[3])
    
    a <- ann_l[[sp]][ann_l[[sp]][,1] == chr & ann_l[[sp]][, 5] >= st & ann_l[[sp]][, 4] <= en, c(10, 4, 5)]
    
    x[12] <- paste(a[,1], sep = "_", collapse = "_")
    x[13] <- paste(a[,2], sep = "_", collapse = "_")
    x[14] <- paste(a[,3], sep = "_", collapse = "_")
    return(t(x))
  }
  ann_b <- as.data.table(t(apply(x,1,fun_ann)))
  
  return(ann_b)
}

# PLOT common ----

centromeric <- read.delim(paste0(baseDir,"/rep/Ann/",ref1Label,".centromere.txt"), header=FALSE) # centromers

chrlen <- read.delim(paste0(baseDir,"/rep/Ann/",ref1Label,".chrs.txt"), header=FALSE) # chrLEN

if (nrow(chrlen) == 17) {chrlen <- chrlen[-nrow(chrlen),]} # remove chrMT    

centromeric$chr <- allChr

ymin <- 0
ymax <- 0.25
for (i in allChr) {
  
  centromeric[centromeric[,"chr"] == i,"ymin"] <- ymin
  centromeric[centromeric[,"chr"] == i,"ymax"] <- ymax
  
  ymin <- ymin + 0.5
  ymax <- ymax + 0.5
}

colnames(chrlen) <- c("chr","len")

ymin <- 0
ymax <- 0.25
for (i in allChr) {
  
  chrlen[chrlen[,"chr"] == i,"ymin"] <- ymin
  chrlen[chrlen[,"chr"] == i,"ymax"] <- ymax
  
  ymin <- ymin + 0.5
  ymax <- ymax + 0.5
}


# body ----

# load mrk and assign allele ----

refs <- c(ref1Label,ref2Label)

mrkdir <- file.path(baseDir, "mrk")

sites <- paste0(baseDir,"/rep/mrktab")

ver_mrk <- fread(file = sites,header = T,data.table = F)

ver_mrk <- ver_mrk[,grep(pattern = paste0(ref1Label,"|",ref2Label), x = colnames(ver_mrk))]

sc_genome <- read.fasta(file = paste0(baseDir,"/ref/",ref1Label,".genome.fa"),set.attributes = F,forceDNAtolower = F)

sp_genome <- read.fasta(file = paste0(baseDir,"/ref/",ref2Label,".genome.fa"),set.attributes = F,forceDNAtolower = F)

ver_mrk$allele_sc <- ""

ver_mrk$allele_sp <- ""

ver_mrk[,5] <- as.character(mapply(function(x,y) sc_genome[[x]][[y]], ver_mrk[,1], ver_mrk[,2]))

ver_mrk[,6] <- as.character(mapply(function(x,y) sp_genome[[x]][[y]], ver_mrk[,3], ver_mrk[,4]))

rm(sc_genome,sp_genome)

mrk_samp_files <- unique(str_split_custom(list.files(mrkdir,pattern = "tab.gz"),1,"\\."))

for (f in mrk_samp_files) {
  
  # load per-sample data ----
  
  print(f)
  
  files <- paste(f,refs,"tab.gz",sep = ".")
  
  mrk_samp <- lapply(files,function(x)fread(paste0(mrkdir,"/",x),sep = "\t",data.table = T))
  
  names(mrk_samp) <- refs
  
  # remove the mrk with rank corresponding to the intro Spar --> Scer ----
  
  rm_rank <- mrk_samp[[ref1Label]][chr == "chrXIV" & POS <= 38500,rank]
  
  mrk_samp[[ref1Label]] <-  mrk_samp[[ref1Label]][!(chr == "chrXIV" & POS <= 38500),]
  
  mrk_samp[[ref2Label]] <- mrk_samp[[ref2Label]][!(is.element(mrk_samp[[ref2Label]][["rank"]],rm_rank)),]
  
  rm(rm_rank)
  
  # butta i mrk i cui raking non sono shared ----
  
   common_ranked_pos <- intersect( mrk_samp[[ref1Label]][,rank],mrk_samp[[ref2Label]][,rank])
   
   mrk_samp[[ref1Label]] <-  mrk_samp[[ref1Label]][match(common_ranked_pos,mrk_samp[[ref1Label]][,rank]),]
   
   mrk_samp[[ref2Label]] <-  mrk_samp[[ref2Label]][match(common_ranked_pos,mrk_samp[[ref2Label]][,rank]),]
   
   rm(common_ranked_pos)
  
  # check, for each mrk position, the allele matches the one on virgin mrk (ALLELES IN THE REFERENCES) ----
  
  setDT(ver_mrk)
  
  ver_mrk[,rank := c(1:nrow(ver_mrk))]
  
  mrk_samp[[ref1Label]][,ASSEMBLYALLELE := ver_mrk[match(mrk_samp[[ref1Label]][,rank],ver_mrk[,rank]),allele_sc] ] 
  
  mrk_samp[[ref2Label]][,ASSEMBLYALLELE := ver_mrk[match(mrk_samp[[ref2Label]][,rank],ver_mrk[,rank]),allele_sp] ] 
  
  boolALTref1 <-  mrk_samp[[ref1Label]][,ALT] == mrk_samp[[ref2Label]][,ASSEMBLYALLELE] | mrk_samp[[ref1Label]][,ALT] == "."
  
  boolALTref2 <-  mrk_samp[[ref2Label]][,ALT] == mrk_samp[[ref1Label]][,ASSEMBLYALLELE] | mrk_samp[[ref2Label]][,ALT] == "."
  
  bool <- boolALTref1 & boolALTref2
  
  mrk_samp[[ref1Label]] <-  mrk_samp[[ref1Label]][bool,]
  
  mrk_samp[[ref2Label]] <-  mrk_samp[[ref2Label]][bool,]
  
  rm(bool,boolALTref1,boolALTref2)
  
  # GT CHECK consistency across couple of assemblies 0/0 must match 1/1 and 0/1 must match 0/1 ---- 
  
  # homo ref1Label
  nM1 <- mrk_samp[[ref1Label]][mrk_samp[[ref1Label]][[5]] == "1/1" & mrk_samp[[ref2Label]][[5]] != "0/0",rank] # to be removed
  nM2 <- mrk_samp[[ref1Label]][mrk_samp[[ref1Label]][[5]] == "0/0" & mrk_samp[[ref2Label]][[5]] != "1/1",rank] # to be removed
  
  # het
  nM3 <- mrk_samp[[ref1Label]][mrk_samp[[ref1Label]][[5]] == "0/1" & mrk_samp[[ref2Label]][[5]] != "0/1",rank] # to be removed
  nM4 <- mrk_samp[[ref1Label]][mrk_samp[[ref1Label]][[5]] != "0/1" & mrk_samp[[ref2Label]][[5]] == "0/1",rank] # to be removed
  
  # homo ref2Label
  nM5 <-mrk_samp[[ref1Label]][mrk_samp[[ref1Label]][[5]] != "1/1" & mrk_samp[[ref2Label]][[5]] == "0/0",rank]  # to be removed
  nM6 <-mrk_samp[[ref1Label]][mrk_samp[[ref1Label]][[5]] != "0/0" & mrk_samp[[ref2Label]][[5]] == "1/1",rank]  # to be removed
  
  nM <- unique(c(nM1,nM2,nM3,nM4,nM5,nM6))
  
  mrk_samp[[ref1Label]] <- mrk_samp[[ref1Label]][-c(match(nM,mrk_samp[[ref1Label]][,rank])),]
  mrk_samp[[ref2Label]] <- mrk_samp[[ref2Label]][-c(match(nM,mrk_samp[[ref2Label]][,rank])),]
  
  # THE FILTERING IS DONE; SAVE ALL THE MARKES and SAMPLES ---- 
  
  out <- paste0(intDir,"/mrk/",f,".",ref1Label,".mrk.filt.gz")
  
  fwrite(x = mrk_samp[[ref1Label]],file = out,append = F,quote = F,sep = "\t",row.names = F,col.names = T,nThread = 2)
  
  out <- paste0(intDir,"/mrk/",f,".",ref2Label,".mrk.filt.gz")
  
  fwrite(x = mrk_samp[[ref2Label]],file = out,append = F,quote = F,sep = "\t",row.names = F,col.names = T,nThread = 2)
  
  # block generation and annotation; final: PLOT ----
  
  mrk_samp <- mrk_samp[[ref1Label]] # keep only one reference! 
  
  # THREE LINE BELOW TAKE ADVANTAGE OF THE RANKING FOR BLOCKING ----
  # spp in questo caso è il genotipo 0/0 e non la species 
  # qui l'idea è: rompi il blocco in caso di diverso chromosoma o GT diverso;
  
  colnames(mrk_samp) <- paste0("col",1:ncol(mrk_samp))

  setDT(mrk_samp)[, grp := rleid(col1, col5)]
  
  df <- mrk_samp[, .(
    start = min(col2),
    end   = max(col2),
    chr   = first(col1),
    spp   = first(col5),
    count_mrk = length(col1)
  ), by = grp]
  
  
  df[df$spp == "0/0","color"] <- "lightblue2"
  df[df$spp == "0/1","color"] <- "gray50"
  df[df$spp == "1/1","color"] <- "red"
  
  df$ymin <- "" 
  df$ymax <- "" 
  
  df$chr <- factor(df$chr, levels=allChr)
  
  ymin <- 0
  ymax <- 0.25
  setDF(df)
  for (i in allChr) {
    
    df[df[,"chr"] == i,"ymin"] <- ymin
    df[df[,"chr"] == i,"ymax"] <- ymax
    
    ymin <- ymin + 0.5
    ymax <- ymax + 0.5
  }
  
  rm(ymin,ymax)
  
  setDT(df)
  
  # colnames(mrk_samp)[6] <- "start"
  # 
  # df[mrk_samp, on = 'start', pSTART := mrk_samp[,2]]
  # 
  # colnames(mrk_samp)[6] <- "end"
  # 
  # df[mrk_samp, on = 'end', pEND := mrk_samp[,2]]
  # 
  # colnames(mrk_samp)[6] <- "rank"
   
  df[,diff:= .(end-start+1)] 
  
  df[,dens:= .(round(abs(count_mrk/diff),digits = 2))] # DENSITY MARKER INSIDE A BLOCK (nMARKERS/PHYSICAL LEN)  
  
  setDF(df)

  blocks <- block_ann(df) # ANNOTATE BLOCK TABLE
  
  blocks <- as.data.frame(apply(blocks, 2, function(x) gsub('\\s+', '', x)))
  
  colnames(blocks) <- c(colnames(df),"gene","gene_start","gene_end")
  
  rm(df)
  
  blocks <- blocks %>% mutate_at(c('diff', 'dens','ymin','ymax'), as.numeric)
  
  blocks <- blocks %>% mutate_at(c('start','end','diff','gene_start','gene_end'), as.integer)
  
  blocks <- blocks %>% arrange(factor(blocks[,4], levels = allChr),start)
  
  # BLUE-RED plot ----
  p <-  ggplot() +
    geom_rect(blocks,mapping=aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax),fill=blocks$color) +
    geom_rect(chrlen, mapping = aes(xmin=0, xmax=len, ymin=ymin,ymax=ymax),fill="grey99",color="black", linewidth=.3,alpha=0.000001) +
    geom_point(centromeric,mapping = aes(x = (centromeric[,1] + centromeric[,2])/2, y=(centromeric[,4] + centromeric[,5])/2),shape=21) +
    annotate(geom="text", x=as.numeric(-50), y=(centromeric[,4] + centromeric[,5])/2, label=allChr ,color="black",hjust=1.5) +
    ggtitle(f) +
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
  
  pPath <- file.path(paste0(baseDir,"/int/",f,".allchr.events.pdf"))
  pdf(file = pPath, width = 16, height = 10)
  print(p)
  dev.off()
  
  blocks[blocks$color == "lightblue2","color"] <- "sc"
  blocks[blocks$color == "red","color"] <- "sp"
  blocks[blocks$color == "gray50","color"] <- "het"
  
  blocks$ymin <- NULL
  
  blocks$ymax <- NULL
  
  colnames(blocks) <- c("grp","first","last","chr","state","count_mrk","species","len","mrk_dens","gene","gene_start","gene_end")
  
  fwrite(x = blocks,file = paste0(baseDir,"/int/",f,".blocks.txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  
  rm(blocks)
  
}
