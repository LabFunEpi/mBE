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


def findOverlappingTFsRaw(filebase):
    """Finds TFs in both _E2 and _control files with given
filebase. Writes out a BED file for every TF/(E2, Control)
combination. Returns the TF list, just in case."""
    
    #    bedE2 = pybedtools.BedTool(filebase+"_E2.bed").to_dataframe().rename(columns={"chrom":"Chromosome", "start":"Start", "end":"End"})
    #    bedControl = pybedtools.BedTool(filebase+"_control.bed").to_dataframe().rename(columns={"chrom":"Chromosome", "start":"Start", "end":"End"})
    
    rawRemap = pybedtools.BedTool(filebase+".bed").to_dataframe().rename(columns={"chrom":"Chromosome", "start":"Start", "end":"End"})
    rawRemap = rawRemap.merge(rawRemap['name'].str.extract('(^\w*)\.(\w*)\.([\w_\-]*)$'), how='left', left_index=True, right_index=True).rename(columns={0:'Source', 1:'TF', 2:"Line"})

    bedE2 = rawRemap[rawRemap['Line']== "MCF-7_E2"]
    bedControl = rawRemap[rawRemap['Line'].isin(["MCF-7", "MCF-7_ETOH"])]

    # Get the list of all of the E2 variants
    #for x in [y for y in rawRemap['Line'].unique() if str(y) != 'nan']:
    #    if(re.search('E2', x)):
    #        print(x)
    
    # Convert the "name" column (which is a mess of parts) to the
    # "Source", "TF", and "Cell Line" data we so crave.
    #bedE2=(bedE2.merge(bedE2['name'].str.extract('(^\w*).(\w*).([\w_]*)'), how='left', left_index=True, right_index=True)).rename(columns={0:'Source', 1:'TF', 2:"Line"})
    #bedControl=(bedControl.merge(bedControl['name'].str.extract('(^\w*).(\w*).([\w_]*)'), how='left', left_index=True, right_index=True)).rename(columns={0:'Source', 1:'TF', 2:"Line"})
    

    # Get a list of the transcription factors that appear in each file
    tf1 = bedE2['TF'].unique()
    tf2 = bedControl['TF'].unique()
    
    
    # Intersect those lists to find the overlap
    tfIntersect = np.intersect1d(tf1, tf2)
    
    # Also find the parts of the list that are unique to each
    tf1Only = np.setdiff1d(tf1, tf2)
    tf2Only = np.setdiff1d(tf2, tf1)
    

    # The lists had better not all be empty. That'd be Bad.
    if((tfIntersect is None) or (tf1Only is None) or (tf2Only is None)):
        print("Three lists are not all non-empty.")
        quit
        

    # For each of the TFs in the intersection, make new BED files of
    # just those peaks from the E2 and Control BEDs
    for tf in tfIntersect:
        dfE2 = bedE2[bedE2['TF']==tf].drop(columns=["Source", "TF", "Line"])
        dfControl = bedControl[bedControl['TF']==tf].drop(columns=["Source", "TF", "Line"])

        pybedtools.BedTool.from_dataframe(dfE2).saveas(filebase+"_E2_"+tf+".bed")
        pybedtools.BedTool.from_dataframe(dfControl).saveas(filebase+"_Control_"+tf+".bed")
        

    # We'll return that list of intersections for the next step, if this gets called. It's only polite.
    return(tfIntersect)




def findOverlappingTFs(filebase):
    """ Finds TFs in both _E2 and _control files with given filebase. Writes out a BED file for every TF/(E2, Control) combination. Returns the TF list, just in case."""
    
    bedE2 = pybedtools.BedTool(filebase+"_E2.bed").to_dataframe().rename(columns={"chrom":"Chromosome", "start":"Start", "end":"End"})
    bedControl = pybedtools.BedTool(filebase+"_control.bed").to_dataframe().rename(columns={"chrom":"Chromosome", "start":"Start", "end":"End"})
 
    
    # Convert the "name" column (which is a mess of parts) to the
    # "Source", "TF", and "Cell Line" data we so crave.
    bedE2=(bedE2.merge(bedE2['name'].str.extract('(^\w*).(\w*).([\w_]*)'), how='left', left_index=True, right_index=True)).rename(columns={0:'Source', 1:'TF', 2:"Line"})
    bedControl=(bedControl.merge(bedControl['name'].str.extract('(^\w*).(\w*).([\w_]*)'), how='left', left_index=True, right_index=True)).rename(columns={0:'Source', 1:'TF', 2:"Line"})
    

    # Get a list of the transcription factors that appear in each file
    tf1 = bedE2['TF'].unique()
    tf2 = bedControl['TF'].unique()
    
    
    # Intersect those lists to find the overlap
    tfIntersect = np.intersect1d(tf1, tf2)
    
    # Also find the parts of the list that are unique to each
    tf1Only = np.setdiff1d(tf1, tf2)
    tf2Only = np.setdiff1d(tf2, tf1)
    

    # The lists had better not all be empty. That'd be Bad.
    if((tfIntersect is None) or (tf1Only is None) or (tf2Only is None)):
        print("Three lists are not all non-empty.")
        quit
        

    # For each of the TFs in the intersection, make new BED files of
    # just those peaks from the E2 and Control BEDs
    for tf in tfIntersect:
        dfE2 = bedE2[bedE2['TF']==tf].drop(columns=["Source", "TF", "Line"])
        dfControl = bedControl[bedControl['TF']==tf].drop(columns=["Source", "TF", "Line"])

        pybedtools.BedTool.from_dataframe(dfE2).saveas(filebase+"_E2_"+tf+".bed")
        pybedtools.BedTool.from_dataframe(dfControl).saveas(filebase+"_Control_"+tf+".bed")
        

    # We'll return that list of intersections for the next step, if this gets called. It's only polite.
    return(tfIntersect)




def intersectE2ControlBEDs(filebase, tfs=None, sanityCheck=False):

    # If we aren't passed the TF list, generate it based on the
    # filenames in this directory.
    if(tfs is None):
        files = glob.glob(filebase+"_E2*.bed")
        tfs = []
        for file in files:
            match = re.search("E2_([\w\-_]+).bed", file)
            if(match):
                tfs.append(match.group(1))
            

    # Looping over TFs, open the BED files and intersect them, both
    # for overlap and for the parts that are unique to each
    # condition.x
    for tf in tfs:
        bedE2 = pybedtools.BedTool(filebase+"_E2_"+tf+".bed").sort().merge()
        bedControl = pybedtools.BedTool(filebase+"_Control_"+tf+".bed").sort().merge()

        bedBoth = bedE2.intersect(bedControl)
        bedE2Only = bedE2.intersect(bedControl, v=True)
        bedControlOnly = bedControl.intersect(bedE2, v=True)

        # This infor is sort of for a sanity check.
        if(sanityCheck):
            print("TF: {tf} \n\tE2 (unique): {e2u}\n\tControl (unique): {cou}\n\tBoth: {both}\n\tE2 (tot): {e2t}\n\tE2 (orig):{e2o}\n\tE2 percent diff: {epdiff}\n\tControl (tot): {cot}\n\tControl (orig) {coo}\n\tControl percent diff: {cpdiff}".format(tf=tf,e2u=len(bedE2Only),epdiff=(len(bedE2Only)+len(bedBoth) - len(bedE2))/(len(bedE2Only)+len(bedBoth) + len(bedE2))*200, cou=len(bedControlOnly), both=len(bedBoth), e2t=len(bedE2Only)+len(bedBoth),  e2o=len(bedE2), cot=len(bedControlOnly)+len(bedBoth), coo=len(bedControl), cpdiff=(len(bedControlOnly)+len(bedBoth) - len(bedControl))/(len(bedControlOnly)+len(bedBoth) + len(bedControl))*200))

        # Save all three intersections
        bedBoth.saveas(filebase+"_Both_"+tf+".bed")
        bedE2Only.saveas(filebase+"_E2Only_"+tf+".bed")
        bedControlOnly.saveas(filebase+"_ControlOnly_"+tf+".bed")

        
    # We'll return that list of intersections for the next step, if this gets called. It's only polite.
    return(tfs)
    




def intersectPeaksBEDs(filebase, tfs=None, baseDir="/", subdirWT="WT", subdirKO="KO"):

    # If we weren't passed the TF list, generate it ourselves from the filenames in this directory.
    if(tfs is None):
        files = glob.glob(filebase+"_E2*.bed")
        tfs = []
        for file in files:
            match = re.search("E2_([\w\-_]+).bed", file)
            if(match):
                tfs.append(match.group(1))

    # read in the matrix files for the two experiment conditions
    bedWT = pybedtools.BedTool(baseDir + "/" + subdirWT + "/filtered_peak_bc_matrix/peaks.bed")
    bedKO = pybedtools.BedTool(baseDir + "/" + subdirKO + "/filtered_peak_bc_matrix/peaks.bed")

    # For each TF, open the relevant BED files and intersect them with the WT and KO matrices
    for tf in tfs:
        bedE2Only = pybedtools.BedTool(filebase+"_E2Only_"+tf+".bed")
        bedControlOnly = pybedtools.BedTool(filebase+"_ControlOnly_"+tf+".bed")
        bedBoth = pybedtools.BedTool(filebase+"_Both_"+tf+".bed")

        # Note that the wa flag is on; this will produce BED files
        # with precisely the same peak locations as in the original
        # experimental matrix files. This is important for later.
        bedWT.intersect(bedE2Only, wa=True).saveas(baseDir + "/" + subdirWT + "/peaks_E2Only_"+tf+".bed")
        bedWT.intersect(bedControlOnly, wa=True).saveas(baseDir + "/" + subdirWT + "/peaks_ControlOnly_"+tf+".bed")
        bedWT.intersect(bedBoth, wa=True).saveas(baseDir + "/" + subdirWT + "/peaks_Both_"+tf+".bed")

        bedKO.intersect(bedE2Only, wa=True).saveas(baseDir + "/" + subdirKO + "/peaks_E2Only_"+tf+".bed")
        bedKO.intersect(bedControlOnly, wa=True).saveas(baseDir + "/" + subdirKO + "/peaks_ControlOnly_"+tf+".bed")
        bedKO.intersect(bedBoth, wa=True).saveas(baseDir + "/" + subdirKO + "/peaks_Both_"+tf+".bed")

        
    # We'll return that list of intersections for the next step, if this gets called. It's only polite.
    return(tfs)
    


        

def mergePeaks(peaksParentDir, tf, condition, tag=None, peaks=None):
    """Merge the peaks matrix with the tf BED file that overlaps the
peaks from previous step. tag set the tag column, peaks allows us to
pass in the peaks matrix data without reloading it over and over again
if this routine is being called in a loop. You're welcome."""

    if((peaks is None)):
        # Get the peaks matrix
        peaks  = readPeaks(peaksParentDir+"filtered_peak_bc_matrix/", tag=tag)


    # Get the overlapping peaks bed file, convert it to a
    # dataframe. The merge() method works here because the peaks have
    # the same locations in the matrix file as in the BED file we
    # produced. And it's hella faster than looping looking for partial
    # overlaps. Believe me.
    peaksOverlap =  pybedtools.BedTool(peaksParentDir+"peaks_"+condition+"_"+tf+".bed").to_dataframe()
    peaksMerged=peaks.merge(peaksOverlap, how="inner", left_on=['Chromosome', 'Start', 'End'], right_on=['chrom', 'start', 'end'])

    return(peaksMerged)


def cellCounts(peaks, tf, sum=False, groupOn="Cell"):
    """Count (or sum) the binding sites in a cell. Reaquires the peaks
dataframe. Keyword sum, when set to True, will use the sum of cutsites
rather than the number of peaks. groupOn (defaults to Cell) sets what
to use to group the cells. Consider "Barcode" instead, maybe."""

    if(sum):
        series = peaks.groupby([groupOn])["CutSites"].sum()
    else:
        series = peaks.groupby([groupOn])["CutSites"].count()


    df = series.to_frame().rename(columns={"CutSites":tf})

    df['Tag'] = peaks['Tag'].iloc[0]

    return(df)



def getStats(tf, condition, peaksWT=None, peaksKO=None, dirWT="/WT/", dirKO = "/KO/"):
    """ Computes the effect size, p-value, and differnetiability for a given TF and condition"""

    
    peaksWT = mergePeaks(dirWT,  tf, condition, peaks=peaksWT, tag="WT")
    dfWT = cellCounts(peaksWT, tf)

    peaksKO = mergePeaks(dirKO,  tf, condition, peaks=peaksKO)
    dfKO = cellCounts(peaksKO, tf)
    
    
    (dummy, pvalue) = compareMotifMediansBetweenSamples(dfWT, dfKO, tf)
    pvalue = pvalue[tf]
    diff = computePeakDiffStat(dfWT, dfKO, tf)
    (eff, description) = computeCohensEffectSize(dfWT, dfKO, tf)

    return((diff, eff, pvalue))
    

def buildStatsList(tfs=None, peaksWT=None, peaksKO=None, verbose=False, remapBase="Remap2020_MCF7_hg38", dirWT="/WT/", dirKO="/KO/"):
    """Builds a dataframe of TFs and conditions and the associated stats"""

    # If the TFs weren't passed, generate the list outselves from the
    # files in this directory.
    if(tfs is None):
        files = glob.glob(remapBase+"E2*.bed")
        tfs = []
        for file in files:
            match = re.search("E2_([\w\-_]+).bed", file)
            if(match):
                tfs.append(match.group(1))

    # The three conditions for the BED files
    conditions = ["E2Only", "ControlOnly", "Both"]

    # Empty lists that we'll fill in a moment.
    effsList = []
    diffsList = []
    pvaluesList = []
    conditionsList = []
    tfsList = []

    # For each condition and TF, work out the stats.
    for condition in conditions:
        if(verbose): print(condition)
        for tf in tfs:
            if(verbose): print("\t {tf}".format(tf=tf))
            # Get the effect size, p-value,  and differentiability.
            (eff, diff, pvalue) = getStats(tf, condition, peaksWT=peaksWT, peaksKO=peaksKO, dirWT=dirWT, dirKO=dirKO)

            # Put things into lists
            effsList.append(eff)
            diffsList.append(diff)
            pvaluesList.append(pvalue)
            conditionsList.append(condition)
            tfsList.append(tf)


    # Make a dataframe.
    df = pd.DataFrame({"TF":tfsList, "Condition":conditionsList, "Eff Size":effsList, "Differentiability":diffsList, "p-value":pvaluesList})

    # This doesn't work and I don't know why.
    #df.replace({"E2Only":"E2 Only", "ControlOnly":"Control Only"})

    return(df)
            


