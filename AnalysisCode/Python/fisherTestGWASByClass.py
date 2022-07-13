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
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import normaltest
import scipy.stats as stats
from sklearn.model_selection import train_test_split
import math
import warnings
import useful

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 


#
# Parse Arguments
#
parser = argparse.ArgumentParser()
parser.add_argument("--cellLine", help="Cell line to test. Options are HMEC, MCF7, HepG2, NHM, and 231L. Defaults to HMEC.", type=str, default="HMEC")
parser.add_argument("--pCrit", help="Critical value of (unadjusted) p-value for annotating the heatmap. Default =0.01.", type=float, default=0.01)
parser.add_argument("--plotRange", help="Extreme (absolute value) values to plot in log heatmap. Default=2.5", type=float, default=2.5)
parser.add_argument("--allEnhancers", help="When set, tests (run separately) for enrichment in all enhancers collectively are also shows.", action="store_true")
# parser.add_argument("--separateGWAS", help="Separately plot heatmaps of the two GWASes rather than together.", action="store_true")

args = parser.parse_args()
# End arguments


# Define input files
variants = ["gwas1", "gwas2"]

# Get size of total genome
totalGenomeBPs = 3137161264	# http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics




# Read in the bed file with the enhancers and covert to a dataframe. (We'll need both the BED and the DF)
enhancerDF = pybedtools.BedTool(args.cellLine + "_enhancers_hg19.bed").to_dataframe().drop(columns=["name", "score", "strand"]).rename({"thickStart":"name"}, axis=1)
# Get the classes
enhancerClasses = enhancerDF["name"].unique()

fileLabel = args.cellLine

if(args.allEnhancers):
    enhancerClasses = np.append(enhancerClasses, ["All Enhancers"])

    allPvalsDF = pd.read_csv("gwasFisherTestPvals.csv", index_col=0)
    allOddsRatiosDF = pd.read_csv("gwasFisherTestOddsRatios.csv", index_col=0)

    fileLabel += "AllEnhancers"


# Create empty dataframes to store results of Fisher Test
pvalsDF = pd.DataFrame(columns = ["GWAS1", "GWAS2"], index=enhancerClasses, dtype=np.float)
oddsRatiosDF = pd.DataFrame(columns = ["GWAS1", "GWAS2"], index=enhancerClasses, dtype=np.float)
logOddsRatiosDF =  pd.DataFrame(columns = ["GWAS1", "GWAS2"], index=enhancerClasses, dtype=np.float)


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
    for enhancerClass in enhancerClasses:
        if(enhancerClass == "All Enhancers"):
            pval = allPvalsDF.loc[args.cellLine, variantName]
            oddsRatio = allOddsRatiosDF.loc[args.cellLine, variantName]
        else:
        
            # Get the sub-DF of just the one enhancer class and turn it into a BED
            enhancerClassDF = enhancerDF[enhancerDF["name"]==enhancerClass]
            enhancerClassBED = pybedtools.BedTool.from_dataframe(enhancerClassDF)


            # Count the number of base-pairs that are enhancers
            enhancerClassBPs = enhancerClassDF["end"].sum() - enhancerClassDF["start"].sum()


            # Intersect the variants with the enhancer locations
            intersectBED = variantBED.intersect(enhancerClassBED, wa=True, wb=True)

            if(intersectBED.count()):
                # Convert to a dataframe and drop and rename various
                # columns. The latter isn't really required, but if I'm going
                # to do any more with this later, I'll want to have tidied up.
                intersectDF = intersectBED.to_dataframe()
                intersectDF = intersectDF.drop(columns=["score", "strand", "thickStart"]).rename({"chrom":"Chrom", "start":"Start", "end":"End", "name":"Name", "thickEnd":"Class"}, axis=1)


                # Count the base-pairs variants that are in enhancers. Again, this is really just th elength of the dataframe...
                variantsInEnhancerClassBPs = intersectDF["End"].sum() - intersectDF["Start"].sum()
            else:
                variantsInEnhancerClassBPs = 0

            # Compute the other elements of the Fisher Test grid.
            variantsNotInEnhancerClassBPs = variantBPs - variantsInEnhancerClassBPs
            enhancerClassNotVariantsBPs = enhancerClassBPs - variantsInEnhancerClassBPs

            nonEnhancerClassNonVariantBPs = totalGenomeBPs - variantBPs - enhancerClassNotVariantsBPs


            # Create the input table and run the Fisher Test
            inputTable = [[variantsInEnhancerClassBPs, variantsNotInEnhancerClassBPs], [enhancerClassNotVariantsBPs, nonEnhancerClassNonVariantBPs]]

            (oddsRatio, pval) = stats.fisher_exact(inputTable)

        # Store the results in the dataframe
        #### This is necessary to get annotations to show up. No idea why.
        if(np.abs(np.log10(oddsRatio))<0.01):
            oddsRatio = 1.0000001


        oddsRatiosDF.loc[enhancerClass, variantName] = oddsRatio
        logOddsRatiosDF.loc[enhancerClass, variantName] = np.log2(oddsRatio)
        pvalsDF.loc[enhancerClass, variantName] = pval
        
        




pvalsDF.to_csv("gwasFisherTestPvals"+fileLabel+".csv")
oddsRatiosDF.to_csv("gwasFisherTestOddsRatios"+fileLabel+".csv")



# if(args.separateGWAS):
    # colsList = ["GWAS1", "GWAS2"]
# else:
    # colsList = [["GWAS1", "GWAS2"]]


# for (i, col) in enumerate(colsList):
    # plt.figure(figsize=(6,14))

    # if(args.separateGWAS):
        # colTmp = [col]
        # data = logOddsRatiosDF.loc[:,[col]]
        # fileLabelTmp = fileLabel+col
    # else:
        # colTmp = col
        # data = logOddsRatiosDF
        # fileLabelTmp = fileLabel

    # annotationsArray = pvalsDF[colTmp].applymap(lambda x:  f"{x:.2E}" if(x<1E-3) else f"{x:.3f}").to_numpy()

    # ax = sns.heatmap(data, cmap='seismic', vmin=-args.plotRange,
                     # vmax=args.plotRange, linewidths=0.0, cbar_kws={'label': 'Log2(Odds Ratio)'},
                     # annot = annotationsArray, fmt='s', annot_kws={"color":"black"})
    # plt.tight_layout()
    
    # # Add black rectangle
    # for i in range(data.shape[0]):
        # for j in range(data.shape[1]):
            # if(pvalsDF.iloc[i,j] < args.pCrit):
                
                # ax.add_patch(Rectangle((j,i), 1, 1, fill=False, edgecolor='black', lw=5))

    # useful.writePDFPlot("gwasFisherTestHeatmap"+fileLabelTmp+"pCrit"+str(args.pCrit)+".pdf")
                
                
