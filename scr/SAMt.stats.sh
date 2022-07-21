#!/bin/bash

# Calculate depth of coverage statistics and check REF name consistency
# between file name and bam file header

BaseDir=/home/nico/expl_pipeline_V2

BamDir=${BaseDir}"/MapDdp"
Depth=${BaseDir}"/Depth"
Nthreads=1

if [[ ! -d $Depth ]]; then
  mkdir -p $Depth
fi

cd $BamDir

for Ind1 in *bam
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  AlnName=$(echo $Ind1 | cut -d "." -f 1,2)
  samtools stats $Ind1 > $Depth/$AlnName.stats.txt
  gzip $Depth/$AlnName.depth.txt
  ) &
  
done
