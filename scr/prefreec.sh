#!/bin/bash

# preparing data for GC content normalization

ref1=$1
ref2=$2
baseDir=$3

ref="$ref1 $ref2"
chrs="chrI chrII chrIII chrIV chrIX chrV chrVI chrVII chrVIII chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI"
refDir=$baseDir"/ref/"

cd $refDir

for i in $ref
do
  prefDir=$baseDir"/cnv/GCdata/"$i
  if [[ ! -d $$prefDir"/Chrref" ]]; then
    mkdir -p $prefDir"/Chrref"
  fi
  # prepara il file LenChr.txt
  grep -v "chrMT" $i.genome.fa.fai | awk 'BEGIN {FS="\t"} {OFS="\t"} {print NR, $1, $2}' > $prefDir"/LenChr.txt"

  for j in $chrs
  do
    # prepara chrFiles
    samtools faidx $i".genome.fa" $j"_"$i > $prefDir"/Chrref/"$j"_"$i".fa"
  done
done
