BaseDir=$1

cd $BaseDir/seq

## Get seq len from R1 fastq
touch $BaseDir/seq/read.len
rm $BaseDir/seq/read.len
touch $BaseDir/seq/read.len

for i in $(ls *.R1.fastq.gz)
do 
name=$(echo $i | cut -d"." -f1)
len=$(zcat $i | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | awk 'BEGIN{max=0}{if(($2)>max)  max=($2)}END {print $1}')
echo -e $name'\t'$len >> read.len
done 

wait
