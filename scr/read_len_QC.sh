#!bin/bash

# 1) mv fastqc output in FASTQC
# 2) unzip the fastqc output
# 3) for each R1 collects the read length
# 4) conserves the results in $BaseDir/seq/read.len
# 5) parser read.len output (keep the longest in case of ranges)

BaseDir=$1

cd $BaseDir/seq/

mkdir -p $BaseDir/seq/FASTQC
wait

mv -f $BaseDir/seq/*zip $BaseDir/seq/FASTQC
mv -f $BaseDir/seq/*html $BaseDir/seq/FASTQC
wait

cd $BaseDir/seq/FASTQC

wait

# unzip
for ind in $(ls *zip)
do
  unzip -o $ind
done

wait

touch $BaseDir/seq/read.len
rm $BaseDir/seq/read.len
touch $BaseDir/seq/read.len

wait

# entra nelle cartelle unzippate e raccogli la len
for ind in $(ls -d *R1_fastqc)
do 
 sample=$(grep Filename $BaseDir/seq/FASTQC/$ind/fastqc_data.txt | cut -f2 | cut -d"." -f 1)
 readlen=$(grep length $BaseDir/seq/FASTQC/$ind/fastqc_data.txt | cut -f2)
 echo -e $sample"\t"$readlen >> $BaseDir/seq/read.len
done
wait 

# Rscript parser 
Rscript $BaseDir/scr/parser_read_len.r $BaseDir > $BaseDir/scr/Logs/parser.read.len.out 2> $BaseDir/scr/Logs/Time.parser.read.len.err
