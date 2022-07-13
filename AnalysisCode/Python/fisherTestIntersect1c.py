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
# import dask.dataframe as dd
import sys
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
# from useful import writePDFPlot

import scipy.stats as stats


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
parser.add_argument('--sortBy', help="Enhancer class to sort by", type=str, default="mBEi")
parser.add_argument('--cellLines', help="Cell line(s)", type=str, action="extend", nargs="+")
parser.add_argument('--cellLineBase', help="Base name for cell lines, pre-appended to all cell lines.", type=str)
parser.add_argument('--bedDir', help="Root dir where BED files are store. Pre-appended to filenames.", type=str, default="/")
parser.add_argument('--tfBed', help="BED file with the TF binding sites", type=str, default="remap2020_MCF-7_all_macs2_hg38_v1_0.bed")
parser.add_argument('--enhancerBed', help="BED file with the enhancer classes", type=str, default="MCF7_K4_sorted.bed")
parser.add_argument('--filebase', help = "Base filename")
#parser.add_argument('--threshold', help="Threshold for number of peaks of a given TF in order to include it.", type=int, default=50)
parser.add_argument('--merge1', help="First set of enhancers to merge into one group", type=str, action="extend", nargs="+")
parser.add_argument('--merge1Name', help='Name to apply to all enhancers in the merge1 group. Default=merge1.', type=str, default='merge1')
parser.add_argument('--TFs', help="List of transcription factors to include. Default is all.", type=str, action="extend", nargs="+")
parser.add_argument('--plotOrder', help="Order of enhancers on the x-axis. Can be the entire set of enhancer classes or a subset. If the latter, the unspecified classes will be at the end, randomly ordered. Default is whatever order they are spit out in.", type=str, action="extend", nargs="+")
parser.add_argument("--dropEnhancers", help="If present, a list of enhancers to drop from the plot", type=str, action ="extend", nargs="+")
parser.add_argument('--extremes', help="Only plot the extreme N highest and lowest TFs", type=int)
parser.add_argument('--tfLimit', help="Limit on how many TFs is too many to label; default = 30", type=int, default=30)
parser.add_argument('--noLogLabel', help="If set, supresses the 'Log2' part of the labels", action="store_true")
parser.add_argument("--trimTFNames", help="Text to remove from start of TF names.", type=str, default=None)
parser.add_argument("--nonCRE", help="If set, includes all of the TF binding sites outside of our CRE Sites", action="store_true")
args = parser.parse_args()


pybedtools.set_tempdir('/tmp/')

#  Input filenames
tfFilename = args.bedDir + args.tfBed
enhancerFilename = args.bedDir + args.enhancerBed

# Open the file with the enhancer locations
enhancerDF = pybedtools.BedTool(enhancerFilename).to_dataframe()
enhancerCols = enhancerDF.columns # Get the column names

# Find the name used for the chromosomes column and the enhancer group column
if("#chrom" in enhancerCols):
    chromName = "#chrom"
elif("chromosome" in enhancerCols):
    chromName = "chromosome"
else:
    chromName = "chrom"

if("deepTools_group" in enhancerCols):
    enhancerName = "deepTools_group"
else:
    enhancerName = "name"
# Rename columns
enhancerDF.rename(columns={chromName:"Chromosome", "start":"Start", "end":"End", enhancerName:"EnhancerClass"}, inplace=True)

# Drop columns except the ones we need
goodCols = ["Chromosome", "Start", "End", "EnhancerClass"]
for col in enhancerDF.columns:
    if(not (col in goodCols)):
        enhancerDF.drop(columns=[col], inplace=True)

# Merge columns if merge1 is set.        
if(args.merge1):
    for enhancer in args.merge1:
        print("Merging ", enhancer)
        enhancerDF.loc[enhancerDF["EnhancerClass"]==enhancer, "EnhancerClass"] = args.merge1Name
    


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


#
# Create a filename
#

# Create the basic filename base from the input files and add the sorting column.
if(args.filebase):
    outbase = "fisherTestByTFSite"+args.filebase
else:
    
    outbase  = "fisherTestByTFSite_"+re.sub(".bed", "", args.tfBed) + "_" + re.sub(".bed", "", args.enhancerBed)

outbase += "_sortOn"+args.sortBy


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
    
    outbase += "_cellLines"
    outbase += "".join(cellLines)
else:
    outbase +="_allLines"

# If we're only taking extreme values, add that to the filename
if(args.extremes):
    outbase += "_extreme"+str(args.extremes)
if(args.nonCRE):
    outbase += "_andCREless"
    

        
    



# Get the list of TF names and enhancer classes available. Allows flexibility
tfs = np.array(tfDF['TF'].unique(), dtype="<U100")
enhancerClasses = np.array(enhancerDF['EnhancerClass'].unique(), dtype="<U10")

if(args.dropEnhancers):
    enhancerClasses = np.setdiff1d(enhancerClasses, args.dropEnhancers)
    outbase += "_DropEnhancers"
    for enhancer in args.dropEnhancers:
        outbase += enhancer

if(args.nonCRE):
    enhancerClasses = np.append(enhancerClasses, ["NonCRE"])
nClasses = enhancerClasses.shape[0]



print("Outbase: " + outbase)

# Bookkeeping. Add the csv to the basename.
outfile = outbase + ".csv"

# We have three possible ways to order these classes. One is just to
# use whatever Python generated above. This is the default if neither
# sortBy nor plotOrder are keywords are set.  The second is to take an
# order from the command line. The third is, if the sortBy is set, to
# move that one to the front of the list and leave the rest as
# default.
print(enhancerClasses)
if(args.plotOrder):
    for (i, enhancer) in enumerate(args.plotOrder):
        # Tricky swap here, needs to record the original location of the swapped element
        oldEnhancerLoc = np.where(enhancerClasses==enhancer)[0][0]
        # No need to do anything if the class is already in the right spot.
        if(i != oldEnhancerLoc):
            enhancerClassesTmp = enhancerClasses[i] # Record what was in the location where moving the class to
            enhancerClasses[i] = enhancer           # Move the class to it's location
            enhancerClasses[oldEnhancerLoc] = enhancerClassesTmp # Put what was there into the classes old location

elif(any(np.isin(np.char.find(enhancerClasses, args.sortBy), 0))):
    # Re-sort the enhancer classes to put the sortBy column first
    
    # Slightly tricky swap here, gonna need a dummy variable for a few lines, plus a stored location for where the sortBy was
    oldSortByLoc = np.where(enhancerClasses==args.sortBy)[0][0]
    # No need to do anything if this class was already first
    if(oldSortByLoc != 0):
        enhancerClassesTmp = enhancerClasses[0]
        enhancerClasses[0] = args.sortBy
        enhancerClasses[oldSortByLoc] = enhancerClassesTmp



# Empty array to handle the items in the loop
enhancerClassArr = np.empty(tfs.shape[0]*(1+nClasses), dtype="<U100")
tfArr = np.empty(tfs.shape[0]*(1+nClasses), dtype="<U100")
inCountArr = np.zeros(tfs.shape[0]*(1+nClasses), dtype=int)
allCount = enhancerDF.shape[0]


df = pd.DataFrame(index=tfs, columns = enhancerClasses)


for tf in tfs:
    tfBED = pybedtools.BedTool.from_dataframe(tfDF[tfDF['TF']==tf])

    tfBindingsCount = 0
    for enhancerClass in enhancerClasses:
        if(enhancerClass == "NonCRE"):
            enhancerBED = pybedtools.BedTool.from_dataframe(enhancerDF)
            df.loc[tf, enhancerClass] = tfBED.count() - tfBED.intersect(enhancerBED, wa=True, u=True).count()
        else:
            enhancerBED = pybedtools.BedTool.from_dataframe(enhancerDF[enhancerDF['EnhancerClass']==enhancerClass])
            df.loc[tf, enhancerClass] = tfBED.intersect(enhancerBED, wa=True, u=True).count()
        

        

df.to_csv(outbase+"RawIntersects.csv")
    


oddsRatioTitles = [x + " Odds Ratio" for x in enhancerClasses]
log2OddsRatioTitles = ["Log2("+x+")" for x in oddsRatioTitles]
pvalueTitles = [x+ " p-value" for x in enhancerClasses]

resultDF = pd.DataFrame(index = tfs, columns = np.append(oddsRatioTitles, pvalueTitles))


allBindings = df.to_numpy().sum()



for tf in tfs:


    #
    # Arg, I changed the identity of tfDF here, so it means something else from here on out.
    #
    
    
    allBindingsOfTF = df.loc[tf, :].sum() # Sum across Row
    
    for className in enhancerClasses:
        allBindingsInEnhancer = df[className].sum()

        tfBindingsInEnhancer = df.loc[tf, className]
        notTFBindingsInEnhancer = allBindingsInEnhancer - tfBindingsInEnhancer
        tfBindingsInOtherClasses = allBindingsOfTF - tfBindingsInEnhancer
        notTFBindingsInOtherClasses = allBindings - allBindingsInEnhancer - tfBindingsInOtherClasses
       
        inputTable = [[tfBindingsInEnhancer, tfBindingsInOtherClasses],
                      [notTFBindingsInEnhancer, notTFBindingsInOtherClasses]]

            #print(tfs[ctr], enhancerClasses[i],  inputTable)
        #print("\t {enhancerClass}\n\t\t {oddsratio: 0.2f}, {pvalue:0.2E}".format(enhancerClass=enhancerClasses[i], oddsratio=oddsRatios[ctr, i], pvalue=pValues[ctr, i]))

        
        (oddsRatio, pValue) = stats.fisher_exact(inputTable)
        resultDF.loc[tf, className+" Odds Ratio"] = oddsRatio
        resultDF.loc[tf, "Log2("+className+" Odds Ratio)"] = np.log2(oddsRatio)
        resultDF.loc[tf, className+" p-value"] = pValue


resultDF = resultDF.astype(float)
        
# Do the stuff for the log(odds ratio) things.

if(args.noLogLabel):
    plotTitles = oddsRatioTitles
else:
    plotTitles = log2OddsRatioTitles

        
# Filter out low-count TFs. Do this after the above to avoid the mismatch between the oddsRatio array and the dataframe
#resultDF = resultDF[resultDF['All TF Peaks']>=args.threshold]


# If we're using k3, we need a less... specific... sort column
if(args.sortBy in enhancerClasses):
    resultDF.sort_values(args.sortBy+ " Odds Ratio", inplace=True)
else:
    resultDF.sort_values(logTitles[0], inplace=True)

                                             

# Save to csv file
resultDF.to_csv(outfile, index_label='TF')


# If we've set the extremes parameter, select the top and bottom N values
if(args.extremes):
    resultDF = resultDF.iloc[list(range(args.extremes)) + list(range(-args.extremes, 0))]


# Get the maximum absolute value of the reference ratios; we'll use that to set the color bar range.
maxRange = np.max(np.abs(resultDF[args.sortBy+ " Odds Ratio"]))


if(args.trimTFNames):
    tmpIndex = []
    for idx in resultDF.index:
        tmpIndex.append(re.sub(args.trimTFNames, "", idx))

    resultDF.index = tmpIndex

# Figure out how many TFs there are and if too many, don't label.
numTFs = resultDF.index.nunique()
if(numTFs > args.tfLimit):
    yLabelVal = False
else:
    yLabelVal = resultDF.index

# Plot the figure and save it
plt.figure(figsize=(6,14))
sns.heatmap(resultDF[plotTitles], vmin = -maxRange, vmax=maxRange,
            cmap='seismic', linewidths=0.0, yticklabels=yLabelVal, cbar_kws={'label': 'Log2(Odds Ratio)'})
plt.xticks(rotation=90) 
plt.tight_layout()
# writePDFPlot(outbase+"HeatMap.pdf", title="Fig 3c in paper")
plt.savefig(outbase+"HeatMap.pdf")

