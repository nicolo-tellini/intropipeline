#!/bin/bash

MainDir=$1

TemplateFile=$MainDir"/scr/Control-FREEC.Config.Template.txt"
Path2CFlib=$MainDir"/scr/"

Nthreads=4
BamDir=$MainDir"/map"
DataExt=".rmd.bam"

# retrieve ploidy from seq/ploidy for plotting
# retrieve read len from seq/reand.len for plotting

for indS in $(cat $MainDir/seq/ploidy | cut -f1)
 do 
 
 Ploidy=$(grep -w $indS $MainDir/seq/ploidy | sort | uniq | cut -f2)
 readlen=$(grep -w $indS $MainDir/seq/read.len | sort | uniq | cut -f2)

 if [[ "$Ploidy" == "NA"  ]]
 then
   continue      
 fi

if [[ ! -s $TemplateFile ]]; then
  echo "Missing template file"
  exit 1
fi

cd $BamDir

for indBAM in $(ls $indS"."*${DataExt})
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
 (
   reference=$(echo $indBAM | cut -d"." -f2)
   refPath=$(echo $MainDir"/ref/"$reference".genome.fa")
   SampleName=$(echo $indBAM | cut -d "." -f 1)

  OutDir=$MainDir"/CNV/results/"$reference/$SampleName
  if [[ ! -d $OutDir ]]; then
    mkdir -p $OutDir
  fi
  
  Config=$MainDir"/CNV/config/"$reference"/"$SampleName
  if [[ ! -d $Config ]]; then
    mkdir -p $Config
  fi
  
  # crea il file config_GC.txt specifico per il campione/reference
  sed "s|STRINGchrLenFile|$MainDir/CNV/GCdata/$reference/LenChr.txt|" < $TemplateFile | sed "s|STRINGchrFiles|$MainDir/CNV/GCdata/$reference/Chrref|" | sed "s|STRINGmateFile|$MainDir/map/$SampleName.$reference.srt.rmd.bam|" | sed "s|STRINGoutputDir|$OutDir|" | sed "s|NumPloidy|$Ploidy|" | sed "s|STRINGgemMappabilityFile|$MainDir/CNV/mappability/$reference/$readlen/$reference.mappability|" | sed -e "s|Sambamba|$MainDir/scr/sambamba|" > $Config"/CF."$SampleName".Input.GC.txt" 

  # running Control-FREEC
  freec -conf $Config"/CF."$SampleName".Input.GC.txt" &> $OutDir"/"$SampleName".log"
  
  # calculate significance
  cat $MainDir/scr/assess_significance.R | R --slave --args $OutDir"/"$SampleName"."$reference".srt.rmd.bam_CNVs" $OutDir"/"$SampleName"."$reference".srt.rmd.bam_ratio.txt"
  
  # plot: in primis, trovo i nomi dei cromosomi, poi li passo allo script R
  AllChr=$(cut -f2 $MainDir"/CNV/GCdata/"$reference"/LenChr.txt" | sed 's|chr||')
  Nchrom=$(echo ${AllChr} | wc -w)
  cat $Path2CFlib"/makeGraph.R" | R --slave --args $Ploidy $OutDir"/"$SampleName"."$reference".srt.rmd.bam_ratio.txt" $Nchrom $AllChr
  ) &
 done
done

wait
