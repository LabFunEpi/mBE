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
import matplotlib
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import normaltest
import scipy.stats as stats
from sklearn.model_selection import train_test_split
import math
import warnings
# from useful import writePDFPlot
# from bioinfokit import visuz 

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 


#
# Parse Arguments
#
parser = argparse.ArgumentParser()
parser.add_argument("--pCrit", help="Critical value of (unadjusted) p-value for annotating the heatmap. Default =0.01.", type=float, default=0.01)
parser.add_argument("--plotRange", help="Extreme (absolute value) values to plot in log heatmap. Default=2.5", type=float, default=2.5)
parser.add_argument("--allEnhancers", help="When set, tests (run separately) for enrichment in all enhancers collectively are also shows.", action="store_true")
parser.add_argument("--yRange", help="Maximum of p-value range. (In log units.)", type=float, default=None)
args = parser.parse_args()
# End arguments


# Define input files
gwases = ["GWAS1", "GWAS2"]
# cellLines = ["HMEC", "MCF7", "231L", "HepG2", "NHM" ]
cellLines = ["HMEC", "MCF7", "231L"]
gwasName = {"GWAS1":"GWAS-1396", "GWAS2":"GWAS-313"}


markers={"HMEC":"o", "MCF7":"+", "231L":"s", "HepG2":"t", "NHM":"*"}
colors = {"mBE":"red", "Active":"Green", "Inactive":"Orange", "ATAConly":"blue", "H3K4me3":"purple", "All Enhancers":"black"}




for gwas in gwases:
    oddsRatiosDF =  pd.DataFrame(columns=["CellLine", "EnhancerClass", "OddsRatio", "pValue"])
    for cellLine in cellLines:
        
        if(args.allEnhancers):
            pvalsLoadDF = pd.read_csv("gwasFisherTestPvals"+cellLine+"AllEnhancers.csv", index_col=0)
            oddsRatiosLoadDF = pd.read_csv("gwasFisherTestOddsRatios"+cellLine+"AllEnhancers.csv", index_col=0)
        else:
            pvalsLoadDF = pd.read_csv("gwasFisherTestPvals"+cellLine+".csv", index_col=0)
            oddsRatiosLoadDF = pd.read_csv("gwasFisherTestOddsRatios"+cellLine+".csv", index_col=0)
        # print(oddsRatiosLoadDF)
        oddsRatiosLoadDF = oddsRatiosLoadDF.drop(["All Enhancers"])
            
            
        #pvalsDF[cellLine] = pvalsLoadDF[gwas]
        #oddsRatiosDF[cellLine] = oddsRatiosLoadDF[gwas].apply(lambda x: np.log2(x))

        for enhancerClass in oddsRatiosLoadDF.index:
            oddsRatiosDF.loc[len(oddsRatiosDF.index)] = [cellLine, enhancerClass, np.log2(oddsRatiosLoadDF.loc[enhancerClass, gwas]), -np.log10(pvalsLoadDF.loc[enhancerClass, gwas])]

        oddsRatiosDF =oddsRatiosDF.replace(-np.inf, -args.plotRange)

    oddsRatiosDF.replace("H3K4me3", "APL", inplace=True)

    #print(oddsRatiosDF)
    plt.figure(figsize=(10,10))


    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 12}
    
    #matplotlib.rc('font', **font)
    
    xrange = 1.1*args.plotRange
    ax = sns.scatterplot(data=oddsRatiosDF, x="OddsRatio", y="pValue", hue="EnhancerClass", style="CellLine",
                              palette=["#19417C", "#FB3D20", "#E375B6", "#278C34", "#e0ad69"], s=300)
    #                         palette=["#19417C", "#FB3D20", "#E375B6", "#278C34", "#e0ad69", "gray"], s=300)
    #                         palette=["darkslateblue", "orangered", "palevioletred", "green", "#e0ad69", "gray"], s=300)
    plt.plot(xrange*np.array([-1,1]), -np.log10(args.pCrit)*np.array([1,1]), 'k-', lw=2)
    plt.title(gwasName[gwas])
    plt.xlim([-xrange, xrange])
    plt.yticks([0, 5, 10, 15])
    yranges = ax.get_ylim()
    if(not args.yRange):
         args.yRange = yranges[1]
    plt.ylim([yranges[0], args.yRange])
    plt.legend(loc='upper left')
    #plt.ylim([0, 10])
    plt.xlabel("log2(Odds Ratio)")
    plt.ylabel("-log10(p-value)")

    
    if(args.allEnhancers):
        plt.savefig("fisherTestVolcano"+gwasName[gwas]+"Pcrit"+str(args.pCrit)+"AndAll.pdf")
        # writePDFPlot("fisherTestVolcano"+gwasName[gwas]+"AndAll.pdf", title="Fisher Test Results from all Cell Lines, "+gwasName[gwas])
    else:
        plt.savefig("fisherTestVolcano"+gwasName[gwas]+"Pcrit"+str(args.pCrit)+".pdf")
        # writePDFPlot("fisherTestVolcano"+gwasName[gwas]+".pdf", title="Fisher Test Results from all Cell Lines, "+gwasName[gwas])

                                       
