#!/bin/bash

BaseDir=$1
nSamples=$2 # number of samples
nThreads=$3 # per-sample number of threads 

cd $BaseDir

cpsDir=$BaseDir"/cps"

mapDir=$BaseDir"/map"

tmpDir=$BaseDir"/tmp"

if [[ ! -d $mapDir ]]; then
  mkdir $mapDir
fi

if [[ ! -d $cpsDir ]]; then
  mkdir $cpsDir
fi

if [[ ! -d $tmpDir ]]; then
  mkdir $tmpDir
fi

if [[ ! -f $cpsDir/cps.txt ]]
then
    echo "cps.txt does not exist."
    touch $cpsDir/cps.txt
    echo "cps.txt created."
    else
    echo "cps.txt exists, nothing to do."
    echo "continue ..."
fi

echo "short read mapping"
for IndS in $(ls -Sr ./seq/*.fastq.gz | cut -d"/" -f3 | cut -d"." -f1 | sed 's+_1++g' | sed 's+_2++g' | sort -u)
do
  # take the sample name.
  # check if it is already stored in cps/cps.txt
  cps=$(grep -w $IndS $BaseDir/cps/cps.txt)
  
  # if $cps is empty the sample will be processed, skipped otherwise. 
  if [ -z "$cps" ]
   then
    for IndR in $(ls ${BaseDir}"/ref/"*.genome.fa)
     do
     ((i++))
     (
     ## check for the presence of paired files 
     revers=$(ls ./seq/$IndS"_2."*"gz")
     refID=$(basename $IndR | cut -d "." -f 1)
     ## run BWA according to the fact the reads are either single or paired end 
     if [ ! -z "$revers" ]
     then 
      #bwa mem -t $nThreads $BaseDir/ref/$refID.genome.fa $BaseDir/seq/$IndS"_1"*".gz" $BaseDir/seq/$IndS"_2"*".gz" -o $mapDir/$IndS"."$refID".sam" 2> /dev/null
      minimap2 -t $nThreads -a -x sr $BaseDir/ref/$refID.genome.fa $BaseDir/seq/$IndS"_1"*".gz" $BaseDir/seq/$IndS"_2"*".gz" -o $mapDir/$IndS"."$refID".sam" 2> /dev/null
     else
      #bwa mem -t $nThreads $BaseDir/ref/$refID.genome.fa $BaseDir/seq/$IndS*".gz" -o $mapDir/$IndS"."$refID".sam" 2> /dev/null
      minimap2 -t $nThreads -a -x sr $BaseDir/ref/$refID.genome.fa $BaseDir/seq/$IndS"_1"*".gz" -o $mapDir/$IndS"."$refID".sam" 2> /dev/null
     fi
     
     ## samtools workflow (http://www.htslib.org/workflow/fastq.html)
     samtools fixmate -@$nThreads -O bam,level=1 -m $mapDir/$IndS"."$refID".sam" $mapDir/$IndS"."$refID".fix.bam"
      
     samtools sort -l 1 -@$nThreads $mapDir/$IndS"."$refID".fix.bam" -T $tmpDir/$IndS"."$refID".tmp.bam" -o $mapDir/$IndS"."$refID".fix.srt.bam"	 
    
     samtools markdup -@$nThreads -O bam,level=1 $mapDir/$IndS"."$refID".fix.srt.bam" $mapDir/$IndS"."$refID".fix.srt.mrk.bam"	 
      
     samtools view -@$nThreads $mapDir/$IndS"."$refID".fix.srt.mrk.bam" -o $mapDir/$IndS"."$refID".bam"
      
     samtools index $mapDir/$IndS"."$refID".bam"
     
     ## comment the line below if you want to keep intermediate files 
     rm -r $mapDir/$IndS"."$refID".sam"  $mapDir/$IndS"."$refID".fix.bam" $mapDir/$IndS"."$refID".fix.srt.bam" $mapDir/$IndS"."$refID".fix.srt.mrk.bam"
     )&
  
  # Job controller allows the efficient use of the threads required. 
  # This is obtained by an if statement that permits the for loop to proceed to the next sample as soon as the number of running samples moves down, below the upper limit imposed by the variable $nSamples.
  if (( i % nSamples == 0 )); then
   wait -n
   i=$(($nSamples-1))
  fi 
  done
  echo $IndS >> $BaseDir/cps/cps.txt
  fi
done
