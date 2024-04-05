<p align="center">
  <img src="https://github.com/nicolo-tellini/intropipeline/blob/loaded/img/intropipeline.logo.png" alt="intropipeline logo" width="450" height="300"/> 
</p>

NEWS:

:rocket:  A v.1.1 with several improvements in stability, speed and memory consumption has been released.

# intropipeline

[![Licence](https://img.shields.io/github/license/nicolo-tellini/intropipeline?style=plastic)](https://github.com/nicolo-tellini/intropipeline/blob/main/LICENSE)
[![Release](https://img.shields.io/github/v/release/nicolo-tellini/intropipeline?style=plastic)](https://github.com/nicolo-tellini/intropipeline/releases/tag/v.1.0.0)
[![release date](https://img.shields.io/github/release-date/nicolo-tellini/intropipeline?color=violet&style=plastic)](https://github.com/nicolo-tellini/intropipeline/releases/tag/v.1.0.0)
[![commit](https://img.shields.io/github/last-commit/nicolo-tellini/intropipeline?color=yellow&style=plastic)](https://github.com/nicolo-tellini/intropipeline/graphs/commit-activity)

An automated computational framework for detecting *Saccharomyces paradoxus* introgressions in *Saccharomyces cerevisiae* strains from paired-end illumina sequencing.

<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/logo2.png" alt="Sublime's custom image"/>
</p>

## Description

v1.0. is described in Tellini, et al. 2024 Nat. EcoEvo, for detecting *S.par* introgressions in *S.cer* strains. <br>

v1.1. contains the following implementations and changes:
- ```minimap2``` replaced ```bwa mem``` almost halving the running time (see [Heng Li 2018, Bioinformatics](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778?login=true)) achieving comparable results; <br>

  sample: ERR3010122 <br>
  
  threads: 2 <br>
  
  Architecture: x86_64 <br>
  
  CPU: Intel(R) Core(TM) i9-9900K CPU @ 3.60GHz <br>
  

  | script  | Elapsed Time | Maximum resident set size (GB) |
  | ------------- | ------------- | ------------- |
  | bwa mem + samtools (v1) | 6:21 (m:ss) | 1.3   |
  | minimap2 + samtools (v1.1) | 3:36 (m:ss) | 1.3  |
  
- improved the reproducibility of the mapping by implementing the standard samtools workflow according to [samtools' guideline](http://www.htslib.org/workflow/fastq.html)
- improved the roboustness of the mapping by appending the name of the strain to a checkpoint (cps) file (```./cps/cps.txt```). The strains which names are stored in ```./cps/cps.txt``` will not be mapped again.
- introduced ```data.table```, ```lapply``` and custom function for large file manipulation for reducing runtime and RAM load.
  example:
  
  sample: ERR3010122 <br>
  
  threads: 2 <br>
  
  Architecture: x86_64 <br>
  
  CPU: Intel(R) Core(TM) i9-9900K CPU @ 3.60GHz <br>

  | script  | Elapsed Time (s) | Maximum resident set size (GB) |
  | ------------- | ------------- | ------------- |
  | parser_marker.r (v1)  | 0:17 s | 0.8 |
  |  parser_marker.r (v1.1)  | 0:06 s | 0.5 |
  | clrs.r (v1) | 0:49 s  | 1.9 |
  | clrs.r (v1.1) | 0:17 s  | 0.7 |
  
- introduced the variables ```nSamples``` and ```nThreads``` inside ```runner.sh```. The first variable controls the number of samples to run in paralell and the second the per-samples number of threads. ```nSamples``` guarantees a contant number of samples running in parallel; as soon as the count drop of one sample an other will start to run. The definition of these variables affect the scripts ```minimap2.sh``` (which replaces ```bwa.sh```), ```bcftools_markers.sh``` (which replaces ```samtools_marker.sh```) and ```freec.sh```;
- corrected an error that prevented the detection of the CNVs;
- Added a new approach for merging markers in blocks:
  
  In v1 the markers are (1) genotyped, (2) filtered and (3) joined as long as they are consecutive and carry the same information. In v1.1 this does not change.

  In v1.1 the markers are (1) ranked, (2) genotyped, (3) filtered, (4) joined as long as they are consecutive in the **ranking** and carry the same information. v1 did not use the ranking.
  Inevitably, this results in a more fragmented signal but provides a more realistic and faithful representation of the introgression reflecting regions where the genotyping was either discordant or failed.
  The ranking also represents the strategy that allowed the speedup of ```clrs.r``` (the script that generates the blocks). 
  
  <p align="center">
  <img src="https://github.com/nicolo-tellini/intropipeline/blob/loaded/img/mrkstrategy.png" alt="Sublime's custom image"/>
</p>

## Download
 
:octocat: :
  
```sh
git clone --recursive https://github.com/nicolo-tellini/intropipeline.git
```

## Content

:open_file_folder: :

```{bash}
.
├── rep
│   ├── Ann
│   └── Asm
├── runner.sh
├── scr
└── seq

5 directories 1 file
```

- ```rep``` : repository with assemblies, annotations and pre-computed marker table,</br>
- ```runner.sh``` : the script you edit and run,</br>
- ```scr``` : scripts,</br>
- ```seq``` : put the FASTQs files here,</br>

### Before starting 

``` gzip -d ./rep/mrktab.gz ```

``` gzip -d ./rep/Asm/*gz```

### About the fastqs 

Move the FASTQs inside ```./seq/```

Paired-end FASTQs data **must** be gziped and suffixed with **.R1.fastq.gz** and **.R2.fastq.gz**.

### Default 

```./scr/bwa.sh``` uses 2 thread for sample (n.samples = 2).

```./scr/samtools_markers.sh``` uses 1 thread for sample (n.samples = 4).

```./scr/gem.sh``` uses 2 threads.

```./scr/freec.sh``` uses 4 threads.

these values can be changed editing the scripts.

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

The results concerning the introgressions are stored in ```./int```

Ex. 

An Alpechin strain:

<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/example-result.png" alt="res"/>
</p>


## How to interprer the result

Blue-Red plots provides an overview of potential introgressed DNA across the genome.
The interpretation of the results is a process that require the integration of different data the pipeline produces.

<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/res1.png" alt="Sublime's custom image"/>
</p>

:exclamation: Reminder: blocks are defined as consecutive markers besring the same genomic info (Homo S.cer, Homo S.par, Het).

<br />

How are markers distributed inside the *S.par* block?

A couple of possible scenarious: 

**Case 1**: abundant markers suporting the block
<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/res2.png" alt="Sublime's custom image"/>
</p>

:exclamation: Note: Only a few markers in the figure above are represented in the cartoon; 


**Case 2**: *not* so abundant markers suporting the block
<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/res3.png" alt="Sublime's custom image"/>
</p>

:exclamation: Note: you should *not* exclude the possibility that a large events is supported by a low number of markers as in the example. 

The number of markers supporting the blocks, the marker density and the info concerning the genotype are stored in ```int``` and ```int/AllSegments```. 

## Dependencies

### Softwares

* [FastQC](https://github.com/s-andrews/FastQC/releases/tag/v0.11.9) 
* minimap2
* samtools
* bcftools
* [GEM](https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/) v. 1.315 (beta)
!! The GEM version used for the analyses is 1.759 (not available anymore). 
* [Control-FREEC](https://github.com/BoevaLab/FREEC/releases/tag/v11.6) v. 11.6; makeGraph.R script was renamed makeplotcnv.R; A copy of all the scripts in [FREEC/scripts/](https://github.com/BoevaLab/FREEC) is in scr. Nevertheless freec has to be installed
* A copy of [sambamba](https://github.com/biod/sambamba/releases/tag/v0.6.5) v. 0.6.5 is provided with the pipeline (no installation required)

### R libraries

* [data.table](https://rdocumentation.org/packages/data.table/versions/1.14.2) 
* [ggplot2](https://github.com/tidyverse/ggplot2/releases/tag/v3.3.5) 
* [rtracklayer](http://www.bioconductor.org/packages/3.11/bioc/html/rtracklayer.html) 
* R.filesets
* GenomicRanges
* purrr
* dplyr
* R.utilis

## Find out more

Marker definition [Methods](https://www.nature.com/articles/s41559-024-02352-5)
 
## Citations

Please cite this paper when using intropipeline for your publications.

> Ancient and recent origins of shared polymorphisms in yeast </br>
> Nicolò Tellini, Matteo De Chiara, Simone Mozzachiodi, Lorenzo Tattini, Chiara Vischioni, Elena S. Naumova, Jonas Warringer, Anders Bergström & Gianni Liti </br>
> Nature Ecologya and Evolution, 2024, https://doi.org/10.1038/s41559-024-02352-5

```
@article{tellini2024ancient,
  title={Ancient and recent origins of shared polymorphisms in yeast},
  author={Tellini, Nicol{\`o} and De Chiara, Matteo and Mozzachiodi, Simone and Tattini, Lorenzo and Vischioni, Chiara and Naumova, Elena S and Warringer, Jonas and Bergstr{\"o}m, Anders and Liti, Gianni},
  journal={Nature Ecology \& Evolution},
  pages={1--16},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```

## Release history

* v1.0 released in 2023
* v1.1 released in 2024
