# Pipeline for clustering ATAC peaks based on ChIP-seq signals of epigenetic markers

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


