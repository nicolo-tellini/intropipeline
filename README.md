<p align="center">
  <img src="https://github.com/nicolo-tellini/intropipeline/blob/loaded/img/intropipeline.logo.png" alt="intropipeline logo" width="450" height="300"/> 
</p>

## 🚀 New Release 2025-10-29

v1.2. contains the following implementations and changes:
- small updates to the runner
- introduced a Mamba environment to simplify the installation of most tools; only GEM requires manual installation
- the generation of blocks has been reverted to v1.0 (see figure below panel on the left), while the ranking from v1.1 is retained solely to facilitate the filtering phase. Briefly, block boundaries are determined only by chromosome changes and genotype, not by ranking. The rationale behind this choice is that the original strategy provides a more realistic and faithful representation of introgression, preserving regions where genotyping failed due to excessive divergence between species. 
- added a script to generate a heatmap of the introgressed blocks,
- simplified outputs in ```int```: for each sample, a PDF and a TXT file are generated containing relevant information about the blocks and their overlap with genes

# intropipeline

[![License](https://img.shields.io/github/license/nicolo-tellini/intropipeline?style=plastic)](https://github.com/nicolo-tellini/intropipeline/blob/main/LICENSE)
[![Release](https://img.shields.io/github/v/release/nicolo-tellini/intropipeline?style=plastic)](https://github.com/nicolo-tellini/intropipeline/releases/tag/v.1.0.0)
[![release date](https://img.shields.io/github/release-date/nicolo-tellini/intropipeline?color=violet&style=plastic)](https://github.com/nicolo-tellini/intropipeline/releases/tag/v.1.0.0)
[![commit](https://img.shields.io/github/last-commit/nicolo-tellini/intropipeline?color=yellow&style=plastic)](https://github.com/nicolo-tellini/intropipeline/graphs/commit-activity)

An automated computational framework for detecting *Saccharomyces paradoxus* introgression in *Saccharomyces cerevisiae* strains from paired-end illumina sequencing.

<p align="center">
  <img src="https://github.com/nicolo-tellini/intropipeline/blob/master/img/scheme.png" alt="scheme" width=1000/>
</p>

## Quick Start

1. Clone the repository:
   ```sh
   git clone --recursive https://github.com/nicolo-tellini/intropipeline.git
   ```
2. Install dependencies:
   ```sh
   mamba create -f intropipeline.yml
   mamba activate intropipeline-env
   ```
3. Install GEM
4. Place your gzipped paired-end FASTQ files in the `seq` directory, named as `*_1.fastq.gz` and `*_2.fastq.gz`.
5. Configure options in `runner.sh`.
6. Run the pipeline:
   ```sh
   nohup bash runner.sh &
   ```
7. Find results in the `int` directory.

## Older versions

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
- improved the robustness of the mapping by appending the name of the strain to a checkpoint (cps) file (```./cps/cps.txt```). The strains which names are stored in ```./cps/cps.txt``` will not be mapped again.
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
  
- introduced the variables ```nSamples``` and ```nThreads``` inside ```runner.sh```. The first variable controls the number of samples to run in paralell and the second the per-samples number of threads. ```nSamples``` guarantees a constant number of samples running in parallel; as soon as the count drop of one sample an other will start to run. The definition of these variables affect the scripts ```minimap2.sh``` (which replaces ```bwa.sh```), ```bcftools_markers.sh``` (which replaces ```samtools_marker.sh```) and ```freec.sh```;
- corrected an error that prevented the detection of the CNVs;
- Added a new approach for merging markers in blocks:
  
  In v1 the markers are (1) genotyped, (2) filtered and (3) joined as long as they are consecutive and carry the same information. In v1.1 this does not change.

  In v1.1 the markers are (1) ranked, (2) genotyped, (3) filtered, (4) joined as long as they are consecutive in the **ranking** and carry the same information. v1 did not use the ranking.
  Inevitably, this results in a more fragmented signal. The ranking also represents the strategy that allowed the speedup of ```clrs.r``` (the script that generates the blocks). 
  
  <p align="center">
  <img src="https://github.com/nicolo-tellini/intropipeline/blob/master/img/mrkstrategy.png" alt="mrk stategy"/>
</p>

v1.0. is described in Tellini, et al. 2024 Nat. EcoEvo, for detecting *S.par* introgressions in *S.cer* strains. <br>

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

### About the fastqs 

Move the FASTQs inside ```./seq/```

Paired-end FASTQs data **must** be gziped and suffixed with **_1.fastq.gz** and **_2.fastq.gz**.

### How to run

Edit ```runner.sh``` :page_with_curl: 

```{bash}
############################
# User configuration
############################

# References
ref1Label="Scc"
ref2Label="CBS432"   # S. paradoxus reference assembly label
ref2="EU"            # Short name for reference

# Resources
nSamples=5
nThreads=8

# Switches
indexing="yes"          # yes|no
fastqQC="yes"           # yes|no|-
shortReadMapping="yes"  # yes|no
mrkgeno="yes"           # yes|no
cnv="yes"               # yes|no
intro="yes"             # yes|no
heatmap="yes"           # yes|no
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


## How to interpret the result

Blue-Red plots provides an overview of potential introgressed DNA across the genome.
The interpretation of the results is a process that require the integration of different data the pipeline produces.

<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/res1.png" alt="Sublime's custom image"/>
</p>

:exclamation: Reminder: blocks are defined as consecutive markers besring the same genomic info (Homo S.cer, Homo S.par, Het).

<br />

How are markers distributed inside the *S.par* block?

A couple of possible scenarious: 

**Case 1**: abundant markers supporting the block
<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/res2.png" alt="Sublime's custom image"/>
</p>

:exclamation: Note: Only a few markers in the figure above are represented in the cartoon; 

**Case 2**: *not* so abundant markers suporting the block
<p align="center">
  <img src="https://github.com/nicolo-tellini/introspect/blob/loaded/img/res3.png" alt="Sublime's custom image"/>
</p>

:exclamation: Note: you should *not* exclude the possibility that a large events is supported by a low number of markers as in the example. 

The number of markers supporting the blocks, the marker density and the info concerning the genotype are stored in the TXT in ```int```. 

## Dependencies

The dependencies are now stored inside ```intropipeline.yml```. [GEM](https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/) can be found and installed at the link. ```sambamba``` is provided in scr dir.

To install the env: ```mamba create -f intropipeline.yml``` and ```mamba activate intropipeline-env```.

## Find out more

Marker definition [Methods](https://www.nature.com/articles/s41559-024-02352-5)
 
## Citations

Please cite this paper when using intropipeline for your publications.
Also, do not forget to cite the papers of the tools used. You can find them inside ```intropipeline.yml``` plus GEM and sambamba.

> Ancient and recent origins of shared polymorphisms in yeast </br>
> Nicolò Tellini, Matteo De Chiara, Simone Mozzachiodi, Lorenzo Tattini, Chiara Vischioni, Elena S. Naumova, Jonas Warringer, Anders Bergström & Gianni Liti </br>
> Nature Ecology and Evolution, 2024, https://doi.org/10.1038/s41559-024-02352-5

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
* v1.2 released in 2025
