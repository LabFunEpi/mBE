<!-- ![image](https://user-images.githubusercontent.com/5864081/173095749-711e59c1-b70c-4bde-8b5d-4adc9217e3ad.png) -->

Read my post about our work here: 
https://cancercommunity.nature.com/posts/macroh2a-histone-variants-modulate-enhancer-activity-to-repress-oncogenic-programs-and-cellular-reprogramming

Mohammed Ismail, W., Mazzone, A., Ghiraldini, F.G. et al. MacroH2A histone variants modulate enhancer activity to repress oncogenic programs and cellular reprogramming. Commun Biol 6, 215 (2023). https://doi.org/10.1038/s42003-023-04571-1

# Pipeline for clustering regulatory elements based on open chromatin and epigenetic marks

## Requirements:
- R 4.0.3 or higher
- R packages: tidyverse, magrittr, cowplot, patchwork, ggpubr
- Bedtools 2.27.1 or higher
- deepTools computeMatrix 3.5.0 or higher

## Installation: 
git clone https://github.com/LabFunEpi/mBE.git

## Running the program: 
Fill in all the parameters (paths to input BED and bigWig files) in parameters.txt and then run mBE_pipeline.sh from Unix command prompt as follows: 
```
./mBE_pipeline.sh parameters.txt
```

## Parameters description
The following parameters are required: 
```
ATAC_peaks=<ATAC peaks in BED file format>
H3K4me1_peaks=<H3K4me1 peaks in BED file format>
blacklist=<ENCODE blacklist file in BED file format; can be downloaded from https://github.com/Boyle-Lab/Blacklist/tree/master/lists>
ENCODE_cCRE=<ENCODE candidate cis-regulatory elements in BED file format; included in this repo>

# The following parameters are ChIP-seq signal files in bigWig format for the 8 epigenetic markers as indicated by the parameter name. 
H3K27ac=
H3K4me1=
H3K4me3=
H2Az=
H3K27me3=
mH2A1=
mH2A2=
CTCF=

k=<number of clusters>
outdir=<directory to store the output files>
```
