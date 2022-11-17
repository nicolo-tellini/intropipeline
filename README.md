# introgression

An automated computational framework for detecting *Saccahromyces paradoxus* introgressions in *Saccahromyces cerevisiae* diploid strains from paired-end illumina sequencing.

<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/logo2.png" alt="Sublime's custom image"/>
</p>

## Description

Pipeline described in Tellini et al 20xx for detecting *S.par* introgressions in *S.cer* diploid strains.

## Download
 
:octocat: :
  
```sh
git clone --recursive https://github.com/nicolo-tellini/introgression.git
cd introgression
```

## Content

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

Run ```runner.sh``` :runner: 

```{bash}
nohup bash runner.sh &
```
## The result



## How to interprer the result

The interpretation of the results is a process that require the integration of different data the pipeline produces.
Blue-Red plots provides an overview of potential introgressed DNA across the genome :

<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/res1.png" alt="Sublime's custom image"/>
</p>

:exclamation: Reminder: blocks are defined as consecutive markers besring the same genomic infos (Homo S.cer, Homo S.par, Het).

How are markers distributed inside the *S.par* block?

A couple of possible scenarious: 

**Case 1**: abundant markers suporting the block
<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/res2.png" alt="Sublime's custom image"/>
</p>
:exclamation: Note: Only a few markers in the figure above are represented in the cartoon; 

**Case 2**: *not* so abundant markers suporting the block
<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/res3.png" alt="Sublime's custom image"/>
</p>

:exclamation: Note: Nevertheless, you should *not* exclude the possibility that a large events is supported by a low number of markers as in the example. 

The number of markers supporting the blocks, the marker density and the info concerning the genotype are stored in ```int``` and ```int/AllSegments```. 

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

## FAQ

**Q**: Why do I have Y/N options in runner.sh ? </br>
**A**: Because if one step fails you want to restart the run from the step failed.

**Q**: Why do I need to specify the *S.par* assembly in the runner.sh if the markers are pre-computed ? </br>
**A**: The pipeline performs 2 mappings: one against *S.c* consensus and one against *S.par* assembly; choosing the *S.par* assembly that fits better the introgressions facilitates the mapping.   

**Q**: How do I interpreter the results ? </br>
**A**: [see here](https://github.com/nicolo-tellini/introspect#how-to-interprer-the-result)

**Q**: What is the reason because there are blocks such as the one reported in **Case 2** ? </br>
**A**: There are several reasons because of these blocks:
   - there are *no* marker positions to genotype because pseudogenes/Ty/tRNA elements in the assemblies generated umbigous alignment,
   - the markers have been genotyped but discarderd because of : 
      - QUAL (QUAL < 20),
      - discordance between the call against *S.c* consensus and *S.par* assembly (different GT or allele),
      - CNV (markers in CNVs are discarded).

**Q**: How do I deal with blocks supported by a scarse number of markers ? Can I trust them ? </br>
**A**: You can filter out the blocks supported by a few consecutive markers. In Tellini et. all 20xx we retaned only blocks supported by 5 consecutive markers. Higher the threshold more relialable the results **but** keep in mind that increasing the threshold you may miss biological relevant introgressions. 

**Q**: I see a single large *S.par* introgression overlapping the subtelomeric region, is it an HGT? </br>
**A**: You cannot discriminate between introgressions and HGTs. Nevertheless, HGTs occur mainly in subtelomeric regions while genome-spread *S.par* blocks may indicate they are the result of the introgression process ([Melania et al. 2018, Nature](https://www.nature.com/articles/s41586-020-2889-1)).

**Q**: Any ways for validating doubt but relevant *S.par* blocks? (there is my favorite gene down there)  </br>
**A**: Some : 
  - Checking both the mappings provides a good indication of what is happening (annotations are in ```ref/Ann```),
  - Check if the markers were available (```rep/mrktab.txt```) and, eventually, at what STEP they were discarded, 
  - Competitive mapping *S.c-S par*, 
  - Phylogenetic methods,
  - further evalutaion with long-read data.

**Q**: Are subelomeric and telomeric regions evaluated ? (there is my favorite gene down there)  </br>
**A**: No, subetelomeric and telomeric regions (as defined in [Jia-Xing Yue et al. 2017, Nature Genetics](https://www.nature.com/articles/ng.3847)) are excluded. Because of the *S.cer* introgression at the beginning of chrXIV ([Liti et al. 2006, Genetics](https://academic.oup.com/genetics/article/174/2/839/6061582#326019337)) in the European *S.par* also this region is excluded. 

## Find out more

*S.cer* consensus assembly [link method paper]
Marker definition [link method paper]
 
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
