#!/usr/bin/env python3
import sys
import argparse
import re
import pandas as pd
import os
from pathlib import Path
import numpy as np
import glob
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import normaltest
import scipy.stats as stats
from sklearn.model_selection import train_test_split
import math
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 


#
# Parse Arguments
#
parser = argparse.ArgumentParser()

args = parser.parse_args()
# End arguments


# Define input files
variants = ["gwas1", "gwas2"]

# cellLines = ["HMEC", "231L", "MCF7", "NHM", "HepG2"]
cellLines = ["HMEC", "231L", "MCF7"]

# Get size of total genome
totalGenomeBPs = 3137161264	# http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics


# Create empty dataframes to store results of Fisher Test
pvalsDF = pd.DataFrame(columns = ["GWAS1", "GWAS2"], index=cellLines)
oddsRatiosDF = pd.DataFrame(columns = ["GWAS1", "GWAS2"], index=cellLines)


# Loop over the variant files
for variant in variants:
    # We'll need the name in uppercase later
    variantName = variant.upper()

    # Read in the bed file as a dataframe in order to clean it up
    variantDF = pybedtools.BedTool(variant+"_hg19.bed").to_dataframe()

    
    # Drop superfluous columns
    for col in variantDF.columns:
        if(not (col in ["chrom", "start", "end", "name"])):
           variantDF.drop(columns=[col], inplace=True)

    # Convert back to a bed object
    variantBED = pybedtools.BedTool.from_dataframe(variantDF)

    # Sum the variants. Technically, each variant is a single
    # base-pair, so we could just use the number of lines. But this
    # doesn't take that long and it generalizes better.
    #variantBPs = variantDF["end"].sum()-variantDF["start"].sum()
    variantBPs = variantDF.shape[0]

    # Loop over the cell lines
    for cellLine in cellLines:

        # Read in the bed file with the enhancers and covert to a dataframe. (We'll need both the BED and the DF)
        enhancerBED = pybedtools.BedTool(cellLine + "_enhancers_hg19.bed")
        enhancerDF = enhancerBED.to_dataframe()

        
        # Count the number of base-pairs that are enhancers
        enhancerBPs = enhancerDF["end"].sum() - enhancerDF["start"].sum()

        # Intersect the variants with the enhancer locations
        intersectBED = variantBED.intersect(enhancerBED, wa=True, wb=True)

        if(intersectBED.count()):
            # Convert to a dataframe and drop and rename various
            # columns. The latter isn't really required, but if I'm going
            # to do any more with this later, I'll want to have tidied up.
            intersectDF = intersectBED.to_dataframe()
            intersectDF = intersectDF.drop(columns=["score", "strand", "thickStart", "itemRgb"]).rename({"chrom":"Chrom", "start":"Start", "end":"End", "name":"Name", "thickEnd":"EnhancerPeak", "blockCount":"Type", "blockSizes":"Class"}, axis=1)


            # Count the base-pairs variants that are in enhancers. Again, this is really just th elength of the dataframe...
            variantsInEnhancersBPs = intersectDF["End"].sum() - intersectDF["Start"].sum()
        else:
            variantsInEnhancersBPs = 0
            

        # Compute the other elements of the Fisher Test grid.
        variantsNotInEnhancersBPs = variantBPs - variantsInEnhancersBPs
        enhancersNotVariantsBPs = enhancerBPs - variantsInEnhancersBPs

        nonEnhancerNonVariantBPs = totalGenomeBPs - variantBPs - enhancersNotVariantsBPs


        # Create the input table and run the Fisher Test
        inputTable = [[variantsInEnhancersBPs, variantsNotInEnhancersBPs], [enhancersNotVariantsBPs, nonEnhancerNonVariantBPs]]
        
        print(variant, cellLine, inputTable)
        (oddsRatio, pval) = stats.fisher_exact(inputTable)


        # Store the results in the dataframe
        oddsRatiosDF.loc[cellLine, variantName] = oddsRatio
        pvalsDF.loc[cellLine, variantName] = pval
        
        





# Save the results to CSV files.
pvalsDF.to_csv("gwasFisherTestPvals.csv")
oddsRatiosDF.to_csv("gwasFisherTestOddsRatios.csv")
