#!/bin/bash

ref1=$1
ref2=$2
BaseDir=$3

Nthreads=2

CopyNumDir=${BaseDir}"/cnv"
MapDir=${BaseDir}"/cnv/mappability"
ModDir=${BaseDir}"/ref"

Allref="${ref1} ${ref2}"

if [[ ! -d ${MapDir} ]]; then
  mkdir -p ${MapDir}
fi

for Indref in ${Allref}
do
  WorkDir=$MapDir"/"$Indref
  if [[ ! -d ${WorkDir} ]]; then
    mkdir -p ${WorkDir}
  fi

  cp $ModDir"/"$Indref".genome.fa" $WorkDir
  cd $WorkDir

 for readLen in $(cat ${BaseDir}/seq/read.len | cut -f2 | uniq)
  do
   
    OutDir=$WorkDir/$readLen
   
    if [[ ! -d $OutDir ]]; then
       mkdir -p $OutDir
    fi

    cd $WorkDir/$readLen
    gem-indexer -T $Nthreads -c dna -i $WorkDir/$Indref".genome.fa" -o $WorkDir/$readLen/$Indref
    gem-mappability -T $Nthreads -I  $WorkDir/$readLen/$Indref".gem" -l $readLen -o  $WorkDir/$readLen/$Indref

 done
  wait
 
 rm -f $Indref".genome.fa"
done
