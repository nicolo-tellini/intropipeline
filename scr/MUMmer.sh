#!/bin/bash

ref1Label=$1
ref2Label=$2
GenosStruct=$3
BaseDir=$4

refDir=${BaseDir}"/ref"

# reference 1 paths
refPair12=$ref1Label"_"$ref2Label
OutDir12=${BaseDir}"/MUMmer/$refPair12"
# reference 2 paths
refPair21=$ref2Label"_"$ref1Label
OutDir21=${BaseDir}"/MUMmer/$refPair21"

if [[ $GenosStruct == "collinear" ]]; then
  # estrae i cromosomi SOLO dal reference 1
  MyChromosomes=$(grep "^>" $refDir/$ref1Label.genome.fa | sed 's|.||' )
  
  # reference 1 calculations
  
  if [[ ! -d $OutDir12 ]]; then
    mkdir -p $OutDir12
  fi
  
  cd $OutDir12
  
  if [[ -e $refPair12.prt.all ]]; then
      rm -f $refPair12.prt.all
  fi
  
  Flag=0
  for IndC in $MyChromosomes
  do
    samtools faidx $refDir/$ref1Label.genome.fa $IndC > $refDir/$ref1Label.genome.$IndC.fa
    samtools faidx $refDir/$ref2Label.genome.fa $IndC > $refDir/$ref2Label.genome.$IndC.fa
    nucmer --mum --prefix=$refPair12.$IndC $refDir/$ref1Label.genome.$IndC.fa $refDir/$ref2Label.genome.$IndC.fa
    if [[ $Flag == 0 ]]; then
      show-snps -ClrT $refPair12.$IndC.delta >> $refPair12.prt.all
      Flag=1
    else
      show-snps -ClrT $refPair12.$IndC.delta | sed '1,4d' >> $refPair12.prt.all
    fi
  done
  
  # separate SNM from indels
  # questi marker (refPair12.prt.all) li uso per filtrare le varianti simulate
  head -4 $refPair12.prt.all > $refPair12.prt.snps
  sed '1,4d' $refPair12.prt.all | awk 'BEGIN{FS="\t"}{if ($2 != "." && $3 != ".") print $0}' >> $refPair12.prt.snps
  head -4 $refPair12.prt.all > $refPair12.prt.indel
  sed '1,4d' $refPair12.prt.all | awk 'BEGIN{FS="\t"}{if ($2 == "." || $3 == ".") print $0}' >> $refPair12.prt.indel
  
  # make header file for Marker.Intersect.R (used for appending results)
  head -4 $refPair12.prt.all > $refPair12.intersect.snps
  sed 's|NUCMER|NUCMER intersection|' $refPair12.intersect.snps > temp.snps
  mv temp.snps $refPair12.intersect.snps
  
  # reference 2 calculations
  
  if [[ ! -d $OutDir21 ]]; then
    mkdir -p $OutDir21
  fi
  
  cd $OutDir21
  
  if [[ -e $refPair21.prt.all ]]; then
    rm -f $refPair21.prt.all
  fi
  
  Flag=0
  for IndC in $MyChromosomes
  do
    nucmer --prefix=$refPair21.$IndC $refDir/$ref2Label.genome.$IndC.fa $refDir/$ref1Label.genome.$IndC.fa
    if [[ $Flag == 0 ]]; then
      show-snps -ClrT $refPair21.$IndC.delta >> $refPair21.prt.all
      Flag=1
    else
      show-snps -ClrT $refPair21.$IndC.delta | sed '1,4d' >> $refPair21.prt.all
    fi
  done
  
  # separate SNM from indels
  # queste indels non servono a una sega
  head -4 $refPair21.prt.all > $refPair21.prt.snps
  sed '1,4d' $refPair21.prt.all | awk 'BEGIN{FS="\t"}{if ($2 != "." && $3 != ".") print $0}' >> $refPair21.prt.snps
  head -4 $refPair21.prt.all > $refPair21.prt.indel
  sed '1,4d' $refPair21.prt.all | awk 'BEGIN{FS="\t"}{if ($2 == "." || $3 == ".") print $0}' >> $refPair21.prt.indel

fi

cd $BaseDir"/scr"
Rscript Marker.Intersect.R $ref1Label $ref2Label $BaseDir


