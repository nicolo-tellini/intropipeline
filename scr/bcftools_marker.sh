#!/bin/bash

ref1Label=$1
ref2Label=$2
BaseDir=$3
nSamples=$4 
nThreads=$5


mrkDir=$BaseDir/mrk
mapDir=$BaseDir/map

if [[ ! -d $mrkDir ]]; then
  mkdir $mrkDir
fi


refs="$ref1Label $ref2Label"

for indR in $refs
do
awk -F'\t' -vcols=chr_$indR,pos_$indR '(NR==1){n=split(cols,cs,",");for(c=1;c<=n;c++){for(i=1;i<=NF;i++)if($(i)==cs[c])ci[c]=i}}{for(i=1;i<=n;i++)printf "%s" FS,$(ci[i]);printf "\n"}' $BaseDir/rep/mrktab | sed 1d > $mrkDir/Sites.$indR.bed
done

BamDir=$BaseDir/map
DataExt=.bam

cd $BamDir

for Ind1 in $(ls *${DataExt})
do
  ((i++)) 
  (

   refName=$(echo $Ind1 | cut -d"." -f2)
   refPath=$BaseDir/ref/$refName.genome.fa
   SampleName=$(echo $Ind1 | cut -d "." -f 1)
   
   bcftools mpileup -E -Ou -q 5 -a DP -f ${BaseDir}/ref/${refName}.genome.fa ${mapDir}/${SampleName}.${refName}.bam | bcftools call -mO z | bcftools view --max-alleles 2 --exclude-types indels | bcftools annotate -x FORMAT/PL,ID,INFO,FILTER | bgzip > $mrkDir"/"${SampleName}"."${refName}".vcf.gz"
   vcftools --gzvcf $mrkDir"/"${SampleName}"."${refName}".vcf.gz" --positions ${mrkDir}/Sites.${refName}.bed --recode --recode-INFO-all --stdout | bgzip -c  > $mrkDir"/"${SampleName}"."${refName}".mrk.vcf.gz"

   rm $mrkDir"/"${SampleName}"."${refName}".vcf.gz"

   mv $mrkDir"/"${SampleName}"."${refName}".mrk.vcf.gz" $mrkDir"/"${SampleName}"."${refName}".vcf.gz"

   bcftools index $mrkDir"/"${SampleName}"."${refName}".vcf.gz"

  ) &

  if (( i % nSamples == 0 )); then
   wait -n
   i=$(($nSamples-1))
  fi

done

wait
