#!/usr/bin/env python3
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
import glob
import pybedtools
import warnings
import argparse
import os
import sys
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
# from useful import writePDFPlot

import scipy.stats as stats

# Same as fisherTestMus (fisher Test on the mice genome), but now I'm
# going to try to collect all three clusters into one plot, just
# showing the KO.

#
# Argument stuff
#

class ExtendAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)
parser = argparse.ArgumentParser()
parser.register('action', 'extend', ExtendAction)

# Parse input arguments
parser.add_argument('--sortBy', help="Cluster to sort by", type=str, default="0")
parser.add_argument('--cellLines', help="Cell line(s)", type=str, action="extend", nargs="+")
parser.add_argument('--cellLineBase', help="Base name for cell lines, pre-appended to all cell lines.", type=str)
parser.add_argument('--bedDir', help="Root dir where BED files are store. Pre-appended to filenames.", type=str, default="/")
parser.add_argument('--tfBed', help="BED file with the TF binding sites", type=str, default="remap2022_nr_macs2_mm10_v1_0.bed")
parser.add_argument('--filebase', help = "Base filename")
#parser.add_argument('--threshold', help="Threshold for number of peaks of a given TF in order to include it.", type=int, default=50)
parser.add_argument('--merge1', help="First set of clusters to merge into one group", type=str, action="extend", nargs="+")
parser.add_argument('--merge1Name', help='Name to apply to all clusters in the merge1 group. Default=merge1.', type=str, default='merge1')
parser.add_argument('--TFs', help="List of transcription factors to include. Default is all.", type=str, action="extend", nargs="+")
parser.add_argument('--plotOrder', help="Order of clusters on the x-axis. Can be the entire set of clusters or a subset. If the latter, the unspecified clusters will be at the end, randomly ordered. Default is whatever order they are spit out in.", type=str, action="extend", nargs="+")
parser.add_argument("--dropClusters", help="If present, a list of clusters to drop from the plot", type=str, action ="extend", nargs="+")
parser.add_argument('--extremes', help="Only plot the extreme N highest and lowest TFs", type=int)
parser.add_argument('--tfLimit', help="Limit on how many TFs is too many to label; default = 30", type=int, default=30)
parser.add_argument('--noLogLabel', help="If set, supresses the 'Log2' part of the labels", action="store_true")
parser.add_argument("--trimTFNames", help="Text to remove from start of TF names.", type=str, default=None)
parser.add_argument("--readRawIntersects", help="Read the raw intersects files rather than computing them. Defaults to off.", action="store_true")
parser.add_argument("--plotOnly", help="Read in the stats file and only redo the plot. Supercedes all of the non-plot keywords, including --readRawIntersects.", action="store_true")
parser.add_argument("--pvalueThreshold", help="Threshold p-value to require for all clusters to plot. If set, all three clusters must have p-values below this threshold to be included. Defaults to none and p-value is ignored.", type=float, default=None)
parser.add_argument("--nClusters", help="Number of clusters", type=int, default=3)
args = parser.parse_args()


# Define a few things
cellTypes = np.array(["WT", "KO"])
clusters = np.array([str(x) for x in range(args.nClusters)])
nCells = cellTypes.size



#
# Create a filename
#


# Create the basic filename base from the input files and add the sorting column.
if(args.filebase):
    outbase = "fishers_" +str(args.nClusters)+"Clusters_"+ args.filebase
else:

    outbase  = "fishers_"+str(args.nClusters)+"Clusters"

# outbase += "_sortOn"+args.sortBy


    # If celllines is set, selected the cell lines.
# if(args.cellLines):
    # outbase += "_cellLines"
    # outbase += "".join(cellLines)
# else:
    # outbase +="_allLines"

print("Outbase: " + outbase, flush=True)

# Bookkeeping. Add the csv to the basename.
outfile = outbase + ".csv"


if(args.plotOnly):
    resultDF = pd.read_csv(outfile, index_col="TF")
    resultDF.dropna(inplace=True)
else:
    if(not args.readRawIntersects):
        tfFilename = os.path.join(args.bedDir, args.tfBed)

        # Open the TF binding sites files. Extract the cell line from the name column. Drop the superflous columns
        tfDF = pybedtools.BedTool(tfFilename).to_dataframe().rename(columns={"chrom":"Chromosome", "start":"Start", "end":"End"})
        if(tfDF['name'].str.match('(^\w*)\.(\w*)\.([\w_\-]*)$')[0]):
            tfDF = tfDF.merge(tfDF['name'].str.extract('(^\w*)\.(\w*)\.([\w_\-]*)$'), how='left', left_index=True, right_index=True).rename(columns={0:'Source', 1:'TF', 2:"Line"})
            tfDF.drop(columns=['name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'], inplace=True)
        elif(tfDF['name'].str.match('(^\w*)\:([\w_\-]*)$')[0]):
            tfDF = tfDF.merge(tfDF['name'].str.extract('(^\w*)\:([\w_\-]*)$'), how='left', left_index=True, right_index=True)
            tfDF.rename(columns={0:'TF', 1:"Line"}, inplace=True)
            tfDF.drop(columns=['name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'], inplace=True)
        elif((tfDF['name'].str.match('^[\w]*_\w*_\w*$')[0]) or (tfDF['name'].str.match('^\w*_\w*$')[0])):
            tfDF.rename(columns={'name':"TF"}, inplace=True)
            tfDF['Line'] = "DF"
        else:
            print("Name column doesn't meet an exisit format structure. Help?")
            sys.exit()




        if(args.TFs):
            tfDF = tfDF[tfDF['TF'].isin(args.TFs)]

            # If celllines is set, selected the cell lines.
            if(args.cellLines):
                cellLines = args.cellLines


                # This bit adds a basic start to all the cell lines, if set.
                if(args.cellLineBase):
                    cellLines = [args.cellLineBase + lineName for lineName in cellLines]
                    tfDF=tfDF[tfDF['Line'].isin(cellLines)]
                    # Add to the output base to make it clear.
                    if(args.cellLines[0][0] in ["-", "_"]):
                        cellLines=[x[1:] for x in args.cellLines]
                    else:
                        cellLines = args.cellLines



        # Get the list of TF names and clusters available. Allows flexibility
        tfs = np.array(tfDF['TF'].unique(), dtype="<U100")
        #tfs = tfs[np.logical_not(np.isnan(tfs))]
        tfs = [x for x in tfs if ((x != "nan") and (x != ""))]
        





    # We have three possible ways to order these clusters. One is just to
    # use whatever Python generated above. This is the default if neither
    # sortBy nor plotOrder are keywords are set.  The second is to take an
    # order from the command line. The third is, if the sortBy is set, to
    # move that one to the front of the list and leave the rest as
    # default.

    if(args.plotOrder):
        for (i, cellType) in enumerate(args.plotOrder):
            # Tricky swap here, needs to record the original location of the swapped element
            oldClusterLoc = np.where(cellTypes==cellType)[0][0]
            # No need to do anything if the cluster is already in the right spot.
            if(i != oldCellLoc):
                cellTmp = cellTypes[i] # Record what was in the location where moving the cluster to
                cellTypes[i] = cellType           # Move the cluster to it's location
                cellTypes[oldCellLoc] = cellTmp # Put what was there into the cluster old location

    elif(any(np.isin(np.char.find(cellTypes, args.sortBy), 0))):
        # Re-sort the clusters to put the sortBy column first

        # Slightly tricky swap here, gonna need a dummy variable for a few lines, plus a stored location for where the sortBy was
        oldSortByLoc = np.where(cellTypes==args.sortBy)[0][0]
        # No need to do anything if this cluster was already first
        if(oldSortByLoc != 0):
            cellTmp = cells[0]
            cells[0] = args.sortBy
            cells[oldSortByLoc] = cellsTmp

    oddsRatioTitles = ["Cluster " + x + " Odds Ratio" for x in clusters]
    log2OddsRatioTitles = ["Log2("+x+")" for x in oddsRatioTitles]
    pvalueTitles = ["Cluster " + x+ " p-value" for x in clusters]

    if(args.readRawIntersects):
        resultDF = None
    else:
        resultDF = pd.DataFrame(index = tfs, columns = np.array(oddsRatioTitles + log2OddsRatioTitles + pvalueTitles))



    # Open the three cluster files
    for cluster in clusters:

        if(args.readRawIntersects):
            df = pd.read_csv(outbase+"Cluster"+cluster+"RawIntersects.csv", index_col=0, dtype={"index":str, "WT":int, "KO":int})
            df.dropna(inplace=True)
            if("nan" in df.index):
                df.drop(index=["nan"], inplace=True)
            tfs = df.index

            if(resultDF is None):
                resultDF = pd.DataFrame(index = tfs, columns = np.array(oddsRatioTitles + log2OddsRatioTitles + pvalueTitles))

        else:
            clusterDF = None

            for cellType in cellTypes:

                # Open file
                newClusterDF = pybedtools.BedTool(os.path.join(args.bedDir, "mus"+str(args.nClusters)+"Clusters."+cluster+"."+cellType+".bed")).to_dataframe()
                newClusterDF["CellType"]=cellType

                if(clusterDF is not None):
                    clusterDF = clusterDF.append(newClusterDF)
                else:
                    clusterDF = newClusterDF


            # Empty array to handle the items in the loop
            df = pd.DataFrame(index=tfs, columns = cellTypes)




            for tf in tfs:
                tfBED = pybedtools.BedTool.from_dataframe(tfDF[tfDF['TF']==tf])

                for cellType in cellTypes:
                    clusterBED = pybedtools.BedTool.from_dataframe(clusterDF[clusterDF['CellType']==cellType])
                    df.loc[tf, cellType] = tfBED.intersect(clusterBED, wa=True, u=True).count()





            df.to_csv(outbase+"Cluster"+cluster+"RawIntersects.csv")




        allBindings = df.to_numpy().sum()



        for tf in tfs:


            allBindingsOfTF = df.loc[tf, :].sum() # Sum across Row



            for cellType in ["KO"]:
                allBindingsInCellType = df[cellType].sum()

                tfBindingsInCellType = df.loc[tf, cellType]
                notTFBindingsInCellType = allBindingsInCellType - tfBindingsInCellType
                tfBindingsInOtherCellTypes = allBindingsOfTF - tfBindingsInCellType
                notTFBindingsInOtherCellTypes = allBindings - allBindingsInCellType - tfBindingsInOtherCellTypes

                inputTable = [[tfBindingsInCellType, tfBindingsInOtherCellTypes],
                              [notTFBindingsInCellType, notTFBindingsInOtherCellTypes]]

                print(cluster, tf, cellType, inputTable, flush=True)

                (oddsRatio, pValue) = stats.fisher_exact(inputTable)
                # print("\t {cellType}\n\t\t {oddsratio: 0.2f}, {pvalue:0.2E}".format(cellType=cellType, oddsratio=oddsRatio, pvalue=pValue), flush=True)


                resultDF.loc[tf, "Cluster " + cluster+" Odds Ratio"] = oddsRatio
                resultDF.loc[tf, "Log2(Cluster "+cluster+" Odds Ratio)"] = np.log2(oddsRatio)
                resultDF.loc[tf, "Cluster " + cluster+" p-value"] = pValue


    resultDF = resultDF.astype(float)


    # Do the stuff for the log(odds ratio) things.

    if(args.noLogLabel):
        plotTitles = oddsRatioTitles
    else:
        plotTitles = log2OddsRatioTitles

    # Filter out low-count TFs. Do this after the above to avoid the mismatch between the oddsRatio array and the dataframe
    #resultDF = resultDF[resultDF['All TF Peaks']>=args.threshold]


    # If we're using k3, we need a less... specific... sort column
    if(args.sortBy in clusters):
        sortCol = "Cluster "  + args.sortBy + " Odds Ratio"
    else:
        sortCol = plotTitles[0]

    resultDF.sort_values(sortCol, inplace=True)


    # Save to csv file
    resultDF.to_csv(outfile, index_label='TF')


# if(args.pvalueThreshold):
    # for cluster in clusters:
        # resultDF = resultDF[resultDF["Cluster " + cluster + " p-value"] <= args.pvalueThreshold]
    
    # outbase += "_pvalueThresh"+str(args.pvalueThreshold)    
# # If we've set the extremes parameter, select the top and bottom N values
# if(args.extremes):
    # if(2*args.extremes <= resultDF.shape[0]):
        # resultDF = resultDF.iloc[list(range(args.extremes)) + list(range(-args.extremes, 0))]


# # Get the maximum absolute value of the reference ratios; we'll use that to set the color bar range.
# maxRange = np.max(np.abs(resultDF["Cluster " + args.sortBy+ " Odds Ratio"]))

# print(maxRange)
# if(args.trimTFNames):
    # tmpIndex = []
    # for idx in resultDF.index:
        # tmpIndex.append(re.sub(args.trimTFNames, "", idx))
        
    # resultDF.index = tmpIndex



    
# # Figure out how many TFs there are and if too many, don't label.
# numTFs = resultDF.index.nunique()

# if(numTFs > args.tfLimit):
    # yLabelVal = False
# else:
    # yLabelVal = resultDF.index


# # Plot the figure and save it
# plt.figure(figsize=(6,14))

# sns.heatmap(resultDF[plotTitles], vmin = -maxRange, vmax=maxRange,
            # cmap='seismic', linewidths=0.0, yticklabels=yLabelVal, cbar_kws={'label': 'Log2(Odds Ratio)'}, rasterized=True)
# plt.xticks(rotation=90) 
# plt.tight_layout()
# # If we're only taking extreme values, add that to the filename
# if(args.extremes):
    # outbase += "_extreme"+str(args.extremes)


# writePDFPlot(outbase+"HeatMap.pdf", title="New fig")


