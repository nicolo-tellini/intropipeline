#!/bin/bash

ref1Label=$1
ref2Label=$2
BaseDir=$3

mrkDir=$BaseDir/mrk

if [[ ! -d $mrkDir ]]; then
  mkdir $mrkDir
fi


refs="$ref1Label $ref2Label"

for indR in $refs
do
awk -F'\t' -vcols=chr_$indR,pos_$indR '(NR==1){n=split(cols,cs,",");for(c=1;c<=n;c++){for(i=1;i<=NF;i++)if($(i)==cs[c])ci[c]=i}}{for(i=1;i<=n;i++)printf "%s" FS,$(ci[i]);printf "\n"}' $BaseDir/rep/mrktab | sed 1d > $mrkDir/Sites.$indR.bed
done

BamDir=$BaseDir/map
DataExt=.srt.rmd.bam

Nthreads=1
Nruns=4

cd $BamDir

for Ind1 in $(ls *${DataExt})
do
   if (( i % Nruns == 0 )); then
    wait
    fi
  ((i++)) 
  (
   refName=$(echo $Ind1 | cut -d"." -f2)
   refPath=$BaseDir/ref/$refName.genome.fa
   SampleName=$(echo $Ind1 | cut -d "." -f 1)
   samtools mpileup --positions $mrkDir/Sites.$refName.bed -u -min-MQ5 --output-tags AD,ADF,ADR,DP,SP --skip-indels --redo-BAQ -f $refPath $Ind1 | bcftools call -c -Oz > $mrkDir"/"$SampleName"."$refName".vcf.gz"
  ) &
done

wait
