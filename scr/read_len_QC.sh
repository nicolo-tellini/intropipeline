#!bin/bash

# 1) mv fastqc output in FASTQC
# 2) unzip the fastqc output
# 3) for each R1 collects the read length
# 4) conserves the results in $BaseDir/seq/read.len
# 5) parser read.len output (keep the longest in case of ranges)

BaseDir=$1

cd $BaseDir/seq/

mkdir FASTQC
wait

mv *zip ./FASTQC/
mv *html ./FASTQC/
wait

cd ./FASTQC/
wait

# unzip
for ind in $(ls *zip)
do
 ( unzip $ind )&
done
wait

touch $BaseDir/seq/read.len
rm $BaseDir/seq/read.len
touch $BaseDir/seq/read.len

wait

# entra nelle cartelle unzippate e raccogli la length
for ind in $(ls -d *R1_fastqc)
do 
 sample=$(grep Filename $BaseDir/seq/FASTQC/$ind/fastqc_data.txt | cut -f2 | cut -d"." -f 1)
 readlen=$(grep length $BaseDir/seq/FASTQC/$ind/fastqc_data.txt | cut -f2)
 echo -e $sample"\t"$readlen >> $BaseDir/seq/read.len
done
wait 

# Rscript parser 
Rscript $BaseDir/scr/parser.read.len.r $BaseDir > $BaseDir/scr/Logs/intro.allchr.plots.out 2> $BaseDir/scr/Logs/Time.intro.allchr.plots.err
