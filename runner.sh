#!/bin/bash

#####################
### user settings ###
#####################

## S. paradoxus reference assembly

ref2Label="CBS432"

## short labels (used to name file)

ref2="EU"

# CHUNCKS of code.
# Suggested default with a few samples (<50)
# Suggested default several samples (>100)

fastqRename="no" # from samp_1.fastq.gz/samp_2.fastq.gz to samp.R1.fastq.gz/samp.R2.fastq.gz
fastqQC="no" # Detailed fastqc control
shortReadMapping="no"
mrkgeno="no"
cnv="no"
intro="yes"

#####################
### settings' end ###
#####################

### part 0

# base folder
BaseDir=$(pwd)

ref1="Scc"
ref1Label="Scc"

# check Logs folder
if [[ ! -d ${BaseDir}/scr/Logs ]]; then mkdir ${BaseDir}/scr/Logs; fi

### part I: references & annotations initialization

#mkdir -p ${BaseDir}/ref/Ann

#cp ${BaseDiir}/rep/Asm/${ref1Label}.genome.fa ${BaseDir}/ref &

#cp ${BaseDir}/rep/Asm/${ref2Label}.genome.fa ${BaseDir}/ref &

#cp ${BaseDir}/rep/Ann/${ref1Label}.* ${BaseDir}/ref/Ann &

#cp ${BaseDir}/rep/Ann/${ref2Label}.* ${BaseDir}/ref/Ann &

#wait

#echo "Preparing the assemblies..."

#for i in $(ls $BaseDir/ref/*fa)
#do
#genome=$(echo $i | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
#sed -i "/^>/s/$/_$genome/" $i
#bwa index $i 2>/dev/null
#samtools faidx $i 2>/dev/null
#done

#wait

## short reads dataset preparation ---- 

if [ $fastqRename == "yes" ] 
then
	echo "Processing the fastqs ..."
	/usr/bin/time -v bash ${BaseDir}/scr/rename_fastqs.sh "${BaseDir}" > ${BaseDir}/scr/Logs/rename_fastqs.out 2> ${BaseDir}/scr/Logs/Time.rename_fastqs.sh.err
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.rename_fastqs.sh.err | cut -d":" -f2 | tr -d "[:blank:]")
	if [ $check == 1 ]
	then
	echo "Exit at rename_fastqs.sh"
	exit
	else
	echo "fastqs rename ... OK"
	fi
elif [ $fastqRename == "no" ]
then
echo "fastqs rename ... skipped"
fi

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

	/usr/bin/time -v bash ${BaseDir}/scr/BWA.Mapping.sh "${BaseDir}" > ${BaseDir}/scr/Logs/BWA.Mapping.out 2> ${BaseDir}/scr/Logs/Time.BWA.Mapping.err
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.BWA.Mapping.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.BWA.Mapping.err)
	if [ $check == 1 ]
	then
	echo "Exit at BWA.mapping.sh"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "bwa failed steps, take a look to scr/Logs/Time.BWA.Mapping.err"
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
	/usr/bin/time -v bash ${BaseDir}/scr/SAMt.Marker.sh "${ref1Label}" "${ref2Label}" "${BaseDir}" > ${BaseDir}/scr/Logs/SAMt.Marker.out 2> ${BaseDir}/scr/Logs/Time.SAMt.Marker.err
        /usr/bin/time -v Rscript ${BaseDir}/scr/Marker.Parser.R "${BaseDir}" > ${BaseDir}/scr/Logs/Marker.Parser.out 2> ${BaseDir}/scr/Logs/Time.Marker.Parser.err
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.SAMt.Marker.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.SAMt.Marker.err)
	if [ $check == 1 ]
	then
	echo "Exit at SAMt.Marker.sh"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "samtools mpileup/call failed steps, take a look to scr/Logs/Time.BWA.Mapping.err"
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
        /usr/bin/time -v bash ${BaseDir}/scr/Mappability.GEM.sh $ref1Label $ref2Label $BaseDir > ${BaseDir}/scr/Logs/Mappability.GEM.out 2> ${BaseDir}/scr/Logs/Time.Mappability.GEM.err

	check=$(grep Exit ${BaseDir}/scr/Logs/Time.Mappability.GEM.err | cut -d":" -f2 | tr -d "[:blank:]")
        failed=$(grep failed ${BaseDir}/scr/Logs/Time.Mappability.GEM.err)
        if [ $check == 1 ]
        then
        echo "Exit at Mappability.GEM.sh"
        exit
        elif [[ ! -z "$failed" ]]
        then
        echo "gem failed, take a look to scr/Logs/Time.Mappability.GEM.err"
        exit
        else
        echo "Mappability ... OK"       
        fi
wait
/usr/bin/time -v bash ${BaseDir}/scr/PreControl-FREEC.sh $ref1Label $ref2Label $BaseDir > ${BaseDir}/scr/Logs/PreControl-FREEC.out 2> ${BaseDir}/scr/Logs/Time.PreControl-FREEC.err

	check=$(grep Exit ${BaseDir}/scr/Logs/Time.PreControl-FREEC.err | cut -d":" -f2 | tr -d "[:blank:]")
        failed=$(grep failed ${BaseDir}/scr/Logs/Time.PreControl-FREEC.err)
        if [ $check == 1 ]
        then
        echo "Exit at PreControl-FREEC.sh"
        exit
        elif [[ ! -z "$failed" ]]
        then
        echo "pre control freec failed, take a look to scr/Logs/Time.PreControl-FREEC.err"
        exit
        else
        echo "preFREEC ... OK"       
        fi
wait
echo "Detecting CNVs ..."
/usr/bin/time -v bash ${BaseDir}/scr/Control-FREEC.sh $BaseDir > /dev/null 2> ${BaseDir}/scr/Logs/Time.Control-FREEC.err
        
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.Control-FREEC.err | cut -d":" -f2 | tr -d "[:blank:]")
        failed=$(grep failed ${BaseDir}/scr/Logs/Time.Control-FREEC.err)
        if [ $check == 1 ]
        then
        echo "Exit at Control-FREEC.sh"
        exit
        elif [[ ! -z "$failed" ]]
        then
        echo "freec failed, take a look to scr/Logs/Time.Control-FREEC.err"
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
        /usr/bin/time -v Rscript ${BaseDir}/scr/ClrS.V8.R $ref1Label $ref2Label $ref1 $ref2 $BaseDir > ${BaseDir}/scr/Logs/ClrS.V8.out 2> ${BaseDir}/scr/Logs/Time.ClrS.V8.err
	
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.ClrS.V8.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.ClrS.V8.err)
	if [ $check == 1 ]
	then
	echo "Exit at ClrS.V8"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "ClrS.V8 failed steps, take a look to scr/Logs/Time.ClrS.V8.err"
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
	echo "AllSeg.R failed steps, take a look to scr/Logs/Time.AllSeg.err"
        exit
	else
	echo "AllSeg ... OK"	
	fi
	wait
	
	/usr/bin/time -v Rscript ${BaseDir}/scr/intro.allchr.plots.r $BaseDir  > ${BaseDir}/scr/Logs/intro.allchr.plots.out 2> ${BaseDir}/scr/Logs/Time.intro.allchr.plots.err
	
	check=$(grep Exit ${BaseDir}/scr/Logs/Time.intro.allchr.plots.err | cut -d":" -f2 | tr -d "[:blank:]")
	failed=$(grep failed ${BaseDir}/scr/Logs/Time.intro.allchr.plots.err)
	if [ $check == 1 ]
	then
	echo "Exit at intro.allchr.plots.r"
	exit
        elif [[ ! -z "$failed" ]]
        then	
	echo "intro chrs plot failed steps, take a look to scr/Logs/Time.intro.allchr.plots.err"
        exit
	else
	echo "IntroGPlot ... OK"	
	fi
else
echo "Inrogression ... skipped"
fi

echo "Thank you for using this pipeline."
