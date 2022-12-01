#!bin/bash

BaseDir=$1

cd $BaseDir/seq/

# run fastqc to take the read length 
for ind in $(ls *R1*)
do
  fastqc $ind
done

wait
