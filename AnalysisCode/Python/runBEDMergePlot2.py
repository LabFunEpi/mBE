import numpy as np
import pandas as pd
import re
import glob
import pybedtools
import warnings

from atacIn import *
from atacOutput import *
from atacStats import *
from bubblePlot import *
from bedMerge import *




#
# Argument stuff
#

# Parse input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', help="Input file")
parser.add_argument('--outfile', help='Write out a pickle file')
parser.add_argument('--csvOut', help='Write out a csv file')
parser.add_argument('--jointTFs', help='Find the TFs in both E2 and Control', action="store_true")
parser.add_argument('--intersectE2Control', help='Intersect the E2 and control BEDs', action="store_true")
parser.add_argument('--intersectPeaks', help='Intersect the peaks matrix', action="store_true")
parser.add_argument('--buildStats', help='Build the stats DF.', action="store_true")
parser.add_argument('--leaveNegatives', help='Include the negative values', action="store_true")
parser.add_argument('--keptCells', help='Only include kept cells', action="store_true")
parser.add_argument('--natCommCells', help="Include kept cells from the NatComm set of June 2021.", action="store_true")
parser.add_argument('--plot', help='Make the bubble plot. Requires either an infile or --buildStats to be true. This defaults to true.', action="store_false")



args = parser.parse_args()


if(args.infile):
    df = pd.read_pickle(args.infile)
else:
    if((not args.buildStats) and (args.plot)):
        print("Plot needs either infile or --buildStats to be set.")
        exit()


# filebase = 'remap2020_MCF-7_all_macs2_hg19_liftedhg38'
filebase = 'remap2022_MCF-7_all_macs2_hg38_v1_0'
tfs=None
if(args.jointTFs):
    tfs=findOverlappingTFsRaw(filebase)
    
if(args.intersectE2Control):
    tfs = intersectE2ControlBEDs(filebase, tfs=tfs, sanityCheck=True)
    
if(args.intersectPeaks):
    tfs = intersectPeaksBEDs(filebase, tfs=tfs, baseDir="/")

    
if(args.buildStats):
    peaksWT   = readPeaksH5("/WT/filtered_peak_bc_matrix.h5", tag="WT")
    peaksKO   = readPeaksH5("/KO/filtered_peak_bc_matrix.h5", tag="KO")

    if(args.keptCells or args.natCommCells):
        if(args.keptCells):
            wtKept = pd.read_csv("/WT_cellskept_hg19.csv")
            koKept = pd.read_csv("/KO_cellskept_hg19.csv")
        else:
            wtKept = pd.read_csv("/WT_Cells.csv")
            koKept = pd.read_csv("/KO_Cells.csv")

        peaksWT = peaksWT.merge(wtKept, how="right", on="Barcode")
        peaksKO = peaksKO.merge(koKept, how="right", on="Barcode")

    
    df = buildStatsList(tfs=tfs, peaksWT=peaksWT, peaksKO=peaksKO, remapBase=filebase,
                        dirWT="/WT/",
                        dirKO="/KO/",
                        verbose=True)


    
    
    df["-Log10 p-value"] = -np.log10(df["p-value"])
    if(args.outfile):
        df.to_pickle(args.outfile)

    if(args.csvOut):
        df.to_csv(args.csvOut)
        
if(args.plot):

    if(args.leaveNegatives):
        vLim=None
    else:
        vLim = (0, 0.14)

    fig = bubblePlot(df, "Eff Size", "-Log10 p-value", "Condition", labelCol="TF",
               clusterNames = ["E2Only", "ControlOnly", "Both"], vertLabel = "Effect Size\nMCF7 mH2A2KO Over wt",
               sizeLabel="-Log10(p-value)", titles=["Estrogen-Specific", "Control-Specific", "Common"],
               colors = ["blue", "orange", "green"], labelFontSize=3, vLim=vLim, sizeFormat="0.0f",
               evenSpacing=True, largeLabels=['ESR1', 'GATA3', 'BRD4', 'FOXA1'], largeFontSize = 10, sizeVals=[10, 20, 30])
    fig.savefig('bubblePlot.pdf')
    # plt.show()


    
