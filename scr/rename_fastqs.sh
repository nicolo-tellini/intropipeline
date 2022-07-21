BaseDir=$1

cd $BaseDir/seq 

for i in $(ls *_1.fastq.gz)
do 

name=$(echo $i | cut -d"_" -f1)
mv $i $name".R1.fastq.gz"

done 

for i in $(ls *_2.fastq.gz)
do 

name=$(echo $i | cut -d"_" -f1)
mv $i $name".R2.fastq.gz"

done 
