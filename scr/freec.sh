#!/bin/bash

baseDir=$1
nThreads=$2

tfile=$baseDir"/scr/Control-FREEC.Config.Template.txt"
p2lib=$baseDir"/scr/"
mapDir=$baseDir"/map"
DataExt=".bam"

# retrieve read len from seq/reand.len for plotting

for ind in $(cat $baseDir/seq/read.len | cut -f1 | sed 's+_1++g' | sed 's+_2++g' | sort -u)
 do

 ploidy=2
 readlen=$( grep $ind $baseDir/seq/read.len | cut -f2 | sort -u)

if [[ ! -s $tfile ]]; then
  echo "Missing template file"
  exit 1
fi

cd $mapDir

for indB in $(ls $ind"."*${DataExt})
do
 if (( i % nThreads == 0 )); then
  wait
 fi
 ((i++))
 (
 reference=$(echo $indB | cut -d"." -f2)
 refPath=$(echo $baseDir"/ref/"$reference".genome.fa")
 SampleName=$(echo $indB | cut -d "." -f 1)

  OutDir=$baseDir"/cnv/results/"$reference/$SampleName
  if [[ ! -d $OutDir ]]; then
    mkdir -p $OutDir
  fi

  Config=$baseDir"/cnv/config/"$reference"/"$SampleName
  if [[ ! -d $Config ]]; then
    mkdir -p $Config
  fi

  # crea il file config_GC.txt specifico per il campione/reference
  sed "s|STRINGchrLenFile|$baseDir/cnv/GCdata/$reference/LenChr.txt|" < $tfile | sed "s|STRINGchrFiles|$baseDir/cnv/GCdata/$reference/Chrref|" | sed "s|STRINGmateFile|$baseDir/map/$SampleName.$reference.bam|" | sed "s|STRINGoutputDir|$OutDir|" | sed "s|NumPloidy|$ploidy|" | sed "s|STRINGgemMappabilityFile|$baseDir/cnv/mappability/$reference/$readlen/$reference.mappability|" | sed -e "s|Sambamba|$baseDir/scr/sambamba|" > $Config"/CF."$SampleName".Input.GC.txt"

  # running Control-FREEC
  freec -conf $Config"/CF."$SampleName".Input.GC.txt" &> $OutDir"/"$SampleName".log"

  # calculate significance
  Rscript $baseDir/scr/assess_significance.R $OutDir"/"$SampleName"."$reference".bam_ratio.txt" $OutDir"/"$SampleName"."$reference".bam_CNVs"

  AllChr=$(cut -f2 $baseDir"/cnv/GCdata/"$reference"/LenChr.txt" | sed 's|chr||')
  Nchrom=$(echo ${AllChr} | wc -w)
  cat $p2lib"/makeplotcnv.R" | R --slave --args $ploidy $OutDir"/"$SampleName"."$reference".bam_ratio.txt" $Nchrom $AllChr
  ) &
 done
done

wait
