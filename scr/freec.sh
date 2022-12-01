#!/bin/bash

BaseDir=$1

TemplateFile=$BaseDir"/scr/Control-FREEC.Config.Template.txt"
Path2CFlib=$BaseDir"/scr/"

Nthreads=4
mapDir=$BaseDir"/map"
DataExt=".rmd.bam"

# retrieve read len from seq/reand.len for plotting

for indS in $(cat $BaseDir/seq/read.len | cut -f1)
 do

 ploidy=2
 readlen=$(grep -w $indS $BaseDir/seq/read.len | sort | uniq | cut -f2)

if [[ ! -s $TemplateFile ]]; then
  echo "Missing template file"
  exit 1
fi

cd $mapDir

for indB in $(ls $indS"."*${DataExt})
do
 if (( i % Nthreads == 0 )); then
  wait
 fi
 ((i++))
 (
 reference=$(echo $indB | cut -d"." -f2)
 refPath=$(echo $BaseDir"/ref/"$reference".genome.fa")
 SampleName=$(echo $indB | cut -d "." -f 1)

  OutDir=$BaseDir"/cnv/results/"$reference/$SampleName
  if [[ ! -d $OutDir ]]; then
    mkdir -p $OutDir
  fi

  Config=$BaseDir"/cnv/config/"$reference"/"$SampleName
  if [[ ! -d $Config ]]; then
    mkdir -p $Config
  fi

  # crea il file config_GC.txt specifico per il campione/reference
  sed "s|STRINGchrLenFile|$BaseDir/cnv/GCdata/$reference/LenChr.txt|" < $TemplateFile | sed "s|STRINGchrFiles|$BaseDir/cnv/GCdata/$reference/Chrref|" | sed "s|STRINGmateFile|$BaseDir/map/$SampleName.$reference.srt.rmd.bam|" | sed "s|STRINGoutputDir|$OutDir|" | sed "s|NumPloidy|$ploidy|" | sed "s|STRINGgemMappabilityFile|$BaseDir/cnv/mappability/$reference/$readlen/$reference.mappability|" | sed -e "s|Sambamba|$BaseDir/scr/sambamba|" > $Config"/CF."$SampleName".Input.GC.txt"

  # running Control-FREEC
  freec -conf $Config"/CF."$SampleName".Input.GC.txt" &> $OutDir"/"$SampleName".log"

  # calculate significance
  cat $BaseDir/scr/assess_significance.R | R --slave --args $OutDir"/"$SampleName"."$reference".srt.rmd.bam_CNVs" $OutDir"/"$SampleName"."$reference".srt.rmd.bam_ratio.txt"

  AllChr=$(cut -f2 $BaseDir"/cnv/GCdata/"$reference"/LenChr.txt" | sed 's|chr||')
  Nchrom=$(echo ${AllChr} | wc -w)
  cat $Path2CFlib"/makeplotcnv.R" | R --slave --args $ploidy $OutDir"/"$SampleName"."$reference".srt.rmd.bam_ratio.txt" $Nchrom $AllChr
  ) &
 done
done

wait
