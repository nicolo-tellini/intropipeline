# introgression <- cambia nome

[]logospace

## Description

An automated computational framework for detecting *S.par* introgressions in *S.cer* diploid strains from paired-end illumina sequencing.

**Scer <-- Spar Introgressions** sheme

Pipeline used in Tellini et al 20xx for detecting *S.par* introgressions in *S.cer* diploid strains. bla bla bla 

## Download
 
:octocat:
  
```sh
git clone --recurse https://github.com/nicolo-tellini/introgression.git
cd introgression
```

## Github Content

:open_file_folder: :

```{bash}
.
├── rep
│   ├── Ann
│   └── Asm
├── runner.sh
├── scr
└── seq

5 directories 1 file
```

- ```rep``` : repository with assemblies, annotations and pre-computed marker table,</br>
- ```runner.sh``` : the script you have to edit and run,</br>
- ```scr``` : scripts,</br>
- ```seq``` : put the FASTQs files here,</br>

### About the fastqs 

Paired-end FASTQs data **must** be gziped and suffixed with **.R1.fastq.gz** and **.R2.fastq.gz**.

### How to run

Edit ```runner.sh``` :page_with_curl: 

```{bash}
#!/bin/bash

#####################
### user settings ###
#####################

## S. paradoxus reference assembly

ref2Label="CBS432" ## choose the Spar assembly you think better fit the origin of your samples

## short labels (used to name file)

ref2="EU" ## choose a short name for Spar

# STEP 1
fastqQC="yes" ## fastqc control (required) ("yes","no" or "-" the last is skip)

# STEP 2
shortReadMapping="yes" ## ("yes","no")

# STEP 3
mrkgeno="yes" ## ("yes","no")

# STEP 4
cnv="yes" ## ("yes","no")

# STEP 5
intro="yes" ## ("yes","no")

#####################
### settings' end ###
#####################
```

Why do I have Y/N options ?

Because if one step fails you want to restart the run from the step the pipeline failed.

Run ```runner.sh``` :runner: 

```{bash}
nohup bash runner.sh &
```


## Dependencies

### Softwares
* [FastQC](https://github.com/s-andrews/FastQC/releases/tag/v0.11.9) v. 0.11.9
* [bwa](https://github.com/lh3/bwa/releases/tag/v0.7.17) v. 0.7.17-r1198-dirty
* [samtools](https://github.com/samtools/samtools/releases/tag/1.9) v. 1.9
* [GEM](https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/) v. 1.315 (beta)
!! The GEM version used for the analyses is 1.759 (not available anymore). 
* [Control-FREEC](https://github.com/BoevaLab/FREEC/releases/tag/v11.6) v. 11.6; makeGraph.R script was renamed makeplotcnv.R; A copy of all the scripts in [FREEC/scripts/](https://github.com/BoevaLab/FREEC) is in scr. Nevertheless freec has to be installed
* A copy of [sambamba](https://github.com/biod/sambamba/releases/tag/v0.6.5) v. 0.6.5 is provided with the pipeline (no installation required)

### R libraries

* [data.table](https://rdocumentation.org/packages/data.table/versions/1.14.2) v. 1.14.2
* [ggplot2](https://github.com/tidyverse/ggplot2/releases/tag/v3.3.5) v. 3.3.5
* [vcfR](https://github.com/knausb/vcfR/releases/tag/v1.12.0) v. 1.12.0
* [scales](https://cran.r-project.org/src/contrib/Archive/scales/) v. 1.1.1
* [rtracklayer](http://www.bioconductor.org/packages/3.11/bioc/html/rtracklayer.html) v. 1.48.0
* [seqinr](https://cran.r-project.org/src/contrib/Archive/seqinr/) v. 4.2-8

## Find out more about 
bla bla bla

# FAQ

## Citations

## Release history

* v1.0.0 released on 2023

## TO-DO list

### short-term updates
- rename columns names ClrS
- contest
- fig. scheme and res example

### long-term updates
- simplify ClrS
- traffic light to solve main bottlenecks 
- heatmap events
- implement [interactive plots](https://nbviewer.org/github/nicolo-tellini/introspect/blob/loaded/introgression-interactive.html)
- software update(?)
