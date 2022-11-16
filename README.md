# SppComp

![Mozilla Add-on](https://img.shields.io/amo/dw/nicolo-tellini?label=Downloads&logo=GitHub)

**SppComp: Saccharomyces species composition**

An automated modular computational framework for a rapid glympse to the species composition of *Saccharomyces* yeasts from large datasets of paired-end illumina reads adtapt to deal different computational resourses.

## Description
  
The species composition of *Saccharomyces* strains plays a major role in biological studies providing valuable insights in the evolutionary history of the *genus* while being exploited for improving indutrial phenotypes. Short-read sequencing are the most popular choice for large-scale genomics projects due to their rapid processesing and affordable prices. As a leading model organisms, the *Saccharomyces* yeasts, heve been massively sequenced with short reads Illumina platfroms. However, integrating several short-read sequencing data from different projects may generate workflow bottlenecks slowing down genomic analysises. Here we present **SppComp**, a modular computational pipeline that takes advantage of the chromosome-level end-to-end genome assemblies produced/reannotated with [LRSDAY](https://github.com/yjx1217/LRSDAY) and competitive short read mapping, as implemented and described in [MuLo-YDH](https://bitbucket.org/lt11/muloydh/src/master/), to assess the species composition of *Saccharomyces* yeasts from large datasets of paired-end illumina reads. **SppComp** is written in bash and R. By means the implementation of state-of-the-art softwares, [functional programming](http://adv-r.had.co.nz/Functional-programming.html), [vectorised code](https://adv-r.hadley.nz/perf-improve.html#vectorise) and effective silly-billy bash tricks, **SppComp** reduces computational slowdowns, allowing a full control of the processes along with the possibility to skip steps, where appropriate.

The mantra is:
> "Make it work, make it right, make it fast."<br />

cit. Kent Beck

<p align="center">
  <img src="https://github.com/nicolo-tellini/sppComb/blob/main/sppComp.flow.png" alt="sppComb flow"/>
</p>

## Download

```sh
git clone https://github.com/nicolo-tellini/sppComb.git
cd sppComb
```
Find out HOW-TO run, [HERE](https://github.com/nicolo-tellini/sppComb/blob/main/howto.md). 

## THE OUTPUT

[CBS2834](https://github.com/nicolo-tellini/sppComb/blob/main/CBS2834_profile.png)

## FIND OUT more about 

:heavy_check_mark: [The Assemblies](https://github.com/nicolo-tellini/sppComb/blob/main/assemblies.md)

:heavy_check_mark: [Competitive mapping](https://academic.oup.com/mbe/article/36/12/2861/5545524?login=false)

:heavy_check_mark: [The Implementation](https://github.com/nicolo-tellini/sppComb/blob/main/deepintothecode.md)

:heavy_check_mark: [Testing](https://github.com/nicolo-tellini/sppComb/blob/main/testing.md)

:heavy_check_mark: [Organisation of the directories](https://github.com/nicolo-tellini/sppComb/blob/main/dirtree.md)

## Release history

* v1.0.0 Released on 2023

## Dependencies

### Softwares
* [samtools](https://github.com/samtools/samtools/releases) v1.14 
* [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)

### R libraries
* [Biostrings](https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
* [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
* [data.table](https://rdocumentation.org/packages/data.table/versions/1.14.2)
* [ggplot2](https://ggplot2.tidyverse.org/)

## Citations
