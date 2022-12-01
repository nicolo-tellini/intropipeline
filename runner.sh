#!/bin/bash

#####################
### user settings ###
#####################

## S. paradoxus reference assembly

ref2Label="CBS432" ## the Spar assembly you think better fit the origin of the introgressions

ref2="EU" ## choose a short name for Spar

## STEP 1
fastqQC="yes" ## fastqc control (required) ("yes","no" or "-" the last is skip)

## STEP 2
shortReadMapping="yes" ## ("yes","no")

## STEP 3
mrkgeno="yes" ## ("yes","no")

## STEP 4
cnv="yes" ## ("yes","no")

## STEP 5
intro="yes" ## ("yes","no")

#####################
### settings' end ###
#####################

BaseDir=$(pwd)

ref1="Scc"

ref1Label="Scc"

# check Logs folder
if [[ ! -d ${BaseDir}/scr/Logs ]]; then mkdir ${BaseDir}/scr/Logs; fi

mkdir -p ${BaseDir}/ref/Ann

cp ${BaseDir}/rep/Asm/${ref1Label}.genome.fa ${BaseDir}/ref &

cp ${BaseDir}/rep/Asm/${ref2Label}.genome.fa ${BaseDir}/ref &

cp ${BaseDir}/rep/Ann/${ref1Label}.* ${BaseDir}/ref/Ann &

cp ${BaseDir}/rep/Ann/${ref2Label}.* ${BaseDir}/ref/Ann &

wait

echo "Preparing the assemblies..."

for i in $(ls $BaseDir/ref/*fa)
do
  genome=$(echo $i | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
  sed -i "/^>/s/$/_$genome/" $i
  bwa index $i 2>/dev/null
  samtools faidx $i 2>/dev/null
done

wait

if [ $fastqQC == "no" ] 
then
	echo "Preparing fastqs ..."
	/usr/bin/time -v bash ${BaseDir}/scr/get_reads_len.sh "${BaseDir}" > ${BaseDir}/scr/Logs/get_reads_len.out 2> ${BaseDir}/scr/Logs/Time.get_reads_len.err
        wait
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.get_reads_len.err | cut -d":" -f2 | tr -d "[:blank:]")
	if [ $check == 1 ]
	then
	echo "Exit at get_reads_len.sh"
	exit
	else
        echo "WARNING: required only read length extraction whitout fastQC report"
	echo "get_reads_len ... OK"
	fi
elif [ $fastqQC == "yes" ]
then
	echo "Preparing fastqs ..."
       	/usr/bin/time -v bash ${BaseDir}/scr/fastQC.sh "${BaseDir}" > ${BaseDir}/scr/Logs/fastQC.out 2> ${BaseDir}/scr/Logs/Time.fastQC.err
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.fastQC.err | cut -d":" -f2 | tr -d "[:blank:]")
	if [ $check == 1 ]
	then
	echo "Exit at fastQC.sh"
	exit
	else
        echo "fastQC reports ... OK"
	fi
       	wait
	
	/usr/bin/time -v bash ${BaseDir}/scr/read_len_QC.sh "${BaseDir}" > ${BaseDir}/scr/Logs/read_len_QC.out 2> ${BaseDir}/scr/Logs/Time.read_len_QC.err
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.read_len_QC.err | cut -d":" -f2 | tr -d "[:blank:]")
	if [ $check == 1 ]
	then
	echo "Exit at read_len_QC.sh"
	exit
	else
        echo "read_len_QC ... OK"
	fi
elif [ $fastqQC == "-" ]
then
	echo "WARNING: read length extraction AND fastQC report ... skipped"
fi

wait

if [ $shortReadMapping == "yes" ] 
then
	echo "short-read mapping ..."

	/usr/bin/time -v bash ${BaseDir}/scr/bwa.sh "${BaseDir}" > ${BaseDir}/scr/Logs/bwa.out 2> ${BaseDir}/scr/Logs/Time.bwa.err
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.bwa.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.bwa.err)
	if [ $check == 1 ]
	then
	echo "Exit at bwa.sh"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "bwa failed steps, take a look to scr/Logs/Time.bwa.err"
        exit
	else
	echo "BWA mapping ... OK"	
	fi

else
echo "BWA mapping ... skipped"
fi

wait

if [ $mrkgeno == "yes" ] 
then
	echo "Genotyping markers ..."
	/usr/bin/time -v bash ${BaseDir}/scr/samtools_marker.sh "${ref1Label}" "${ref2Label}" "${BaseDir}" > ${BaseDir}/scr/Logs/samtools_marker.out 2> ${BaseDir}/scr/Logs/Time.samtools_marker.err
        /usr/bin/time -v Rscript ${BaseDir}/scr/parser_marker.R "${BaseDir}" > ${BaseDir}/scr/Logs/parser_marker.out 2> ${BaseDir}/scr/Logs/Time.parser_marker.err
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.samtools_marker.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.samtools_marker.err)
	if [ $check == 1 ]
	then
	echo "Exit at samtools_marker.sh"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "samtools mpileup/call failed steps, take a look to scr/Logs/Time.samtools_marker.err"
        exit
	else
	echo "markers geno ... OK"	
	fi
else
echo "markers geno ... skipped"
fi

wait

if [ $cnv == "yes" ]
then
	echo "Preparing CNVs ..."
        /usr/bin/time -v bash ${BaseDir}/scr/gem.sh $ref1Label $ref2Label $BaseDir > ${BaseDir}/scr/Logs/gem.out 2> ${BaseDir}/scr/Logs/Time.gem.err

	check=$(grep Exit ${BaseDir}/scr/Logs/Time.gem.err | cut -d":" -f2 | tr -d "[:blank:]")
        failed=$(grep failed ${BaseDir}/scr/Logs/Time.gem.err)
        if [ $check == 1 ]
        then
        echo "Exit at gem.sh"
        exit
        elif [[ ! -z "$failed" ]]
        then
        echo "gem failed, take a look to scr/Logs/Time.gem.err"
        exit
        else
        echo "Mappability ... OK"       
        fi
wait
/usr/bin/time -v bash ${BaseDir}/scr/prefreec.sh $ref1Label $ref2Label $BaseDir > ${BaseDir}/scr/Logs/prefreec.out 2> ${BaseDir}/scr/Logs/Time.prefreec.err

	check=$(grep Exit ${BaseDir}/scr/Logs/Time.prefreec.err | cut -d":" -f2 | tr -d "[:blank:]")
        failed=$(grep failed ${BaseDir}/scr/Logs/Time.prefreec.err)
        if [ $check == 1 ]
        then
        echo "Exit at prefreec.sh.sh"
        exit
        elif [[ ! -z "$failed" ]]
        then
        echo "pre control freec failed, take a look to scr/Logs/Time.prefreec.sh.err"
        exit
        else
        echo "preFREEC ... OK"       
        fi
wait
echo "Detecting CNVs ..."
/usr/bin/time -v bash ${BaseDir}/scr/freec.sh $BaseDir > /dev/null 2> ${BaseDir}/scr/Logs/Time.freec.err
        
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.freec.err | cut -d":" -f2 | tr -d "[:blank:]")
        failed=$(grep failed ${BaseDir}/scr/Logs/Time.freec.err)
        if [ $check == 1 ]
        then
        echo "Exit at freec.sh"
        exit
        elif [[ ! -z "$failed" ]]
        then
        echo "freec failed, take a look to scr/Logs/Time.freec.err"
        exit
        else
        echo "FREEC ... OK"       
        fi
else
echo "CNVs ... skipped"
fi

if [ $intro == "yes" ] 
then
        echo "Detecting segments ..."
        /usr/bin/time -v Rscript ${BaseDir}/scr/ClrS_V8.R $ref1Label $ref2Label $ref1 $ref2 $BaseDir > ${BaseDir}/scr/Logs/ClrS_V8.out 2> ${BaseDir}/scr/Logs/Time.ClrS_V8.err
	
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.ClrS_V8.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.ClrS_V8.err)
	if [ $check == 1 ]
	then
	echo "Exit at ClrS_V8"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "ClrS.V8 failed, take a look to scr/Logs/Time.ClrS_V8.err"
        exit
	else
	echo "ClrS ... OK"	
	fi
        wait
	
	/usr/bin/time -v Rscript ${BaseDir}/scr/AllSeg.R $ref1Label $ref2Label $ref1 $ref2 $BaseDir > ${BaseDir}/scr/Logs/AllSeg.out 2> ${BaseDir}/scr/Logs/Time.AllSeg.err
	
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.AllSeg.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.AllSeg.err)
	if [ $check == 1 ]
	then
	echo "Exit at AllSeg.R"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "AllSeg.R failed, take a look to scr/Logs/Time.AllSeg.err"
        exit
	else
	echo "AllSeg ... OK"	
	fi
	wait
	
	/usr/bin/time -v Rscript ${BaseDir}/scr/intro_allchr_plots.R ${BaseDir}  > ${BaseDir}/scr/Logs/intro_allchr_plots.out 2> ${BaseDir}/scr/Logs/Time.intro_allchr_plots.err
	
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.intro_allchr_plots.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.intro_allchr_plots.err)
	if [ $check == 1 ]
	then
	echo "Exit at intro.allchr.plots.r"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "intro failed, take a look to scr/Logs/Time.intro.allchr.plots.err"
        exit
	else
	echo "red-blue plots ... OK"	
	fi
else
echo "Inrogression ... skipped"
fi

if [ $heatmap == "yes" ] 
then
	
	/usr/bin/time -v Rscript ${BaseDir}/scr/heatmap_events.R ${BaseDir} ${ref1} > ${BaseDir}/scr/Logs/heatmap_events.out 2> ${BaseDir}/scr/Logs/Time.heatmap_events.err
	
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.heatmap_events.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.heatmap_events.err)
	if [ $check == 1 ]
	then
	echo "Exit at heatmap_events.R"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "heatmap failed, take a look to scr/Logs/Time.heatmap_events.err"
        exit
	else
	echo "Heatmap ... OK"	
	fi
else
echo "Heatmap ... skipped"
fi

echo "Thank you for using the pipeline."
