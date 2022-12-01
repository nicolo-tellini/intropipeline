#!/bin/bash

BaseDir=$1
Nthreads=2
Nruns=2

cd $BaseDir

OutDir="map"

if [[ ! -d $OutDir ]]; then
  mkdir $OutDir
fi

for IndS in $(ls seq/*gz | cut -d"." -f1 | cut -d"/" -f2 | sort | uniq)
do
  for IndR in $(ls ${BaseDir}"/ref/"*.genome.fa)
  do
   if (( i % Nruns == 0 )); then
   wait
   fi
  ((i++)) 
  (
    refID=$(basename $IndR | cut -d "." -f 1)
    bwa mem -M -t $Nthreads $IndR seq/$IndS".R1"*".gz" seq/$IndS".R2"*".gz" 2>/dev/null | samtools view -bS -F 1036 > $OutDir/$IndS.$refID.bam  
    samtools sort -o $OutDir/$IndS.$refID.srt.bam -@ $Nthreads $OutDir/$IndS.$refID.bam
    rm $OutDir/$IndS.$refID.bam
    samtools rmdup $OutDir/$IndS.$refID.srt.bam $OutDir/$IndS.$refID.srt.rmd.bam
    rm $OutDir/$IndS.$refID.srt.bam
    samtools index $OutDir/$IndS.$refID.srt.rmd.bam
  )&
done
done
wait
