# introgression <- cambia nome

Pipeline used in Tellini et al 20xx for detecting *S.par* introgressions in *S.cer* diploid strains.

**Scer <-- Spar Introgressions**

An automated computational framework for detecting *S.par* introgressions in *S.cer* diploid strains from large datasets of paired-end illumina reads.

## Description

## Download

```sh
https://github.com/nicolo-tellini/introgression.git
cd introgression
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

## Citations

## Release history

* v1.0.0 released on 2023
