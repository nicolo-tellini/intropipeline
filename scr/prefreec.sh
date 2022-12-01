#!/bin/bash

# preparing data for GC content normalization

ref1=$1
ref2=$2
MainDir=$3

ref="$ref1 $ref2"
MyChr="chrI chrII chrIII chrIV chrIX chrV chrVI chrVII chrVIII chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI"
WorkDir=$MainDir"/ref/"

cd $WorkDir

for Ind1 in $ref
do
  PathrefDir=$MainDir"/cnv/GCdata/"$Ind1
  if [[ ! -d $$PathrefDir"/Chrref" ]]; then
    mkdir -p $PathrefDir"/Chrref"
  fi
  # prepara il file LenChr.txt
  grep -v "chrMT" $Ind1.genome.fa.fai | awk 'BEGIN {FS="\t"} {OFS="\t"} {print NR, $1, $2}' > $PathrefDir"/LenChr.txt"

  for Ind2 in $MyChr
  do
    # prepara chrFiles
    samtools faidx $Ind1".genome.fa" $Ind2"_"$Ind1 > $PathrefDir"/Chrref/"$Ind2"_"$Ind1".fa"
  done
done
