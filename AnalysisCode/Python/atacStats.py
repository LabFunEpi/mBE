#
# Stats routines for ATAC data
#

import argparse
import scipy.stats as spstats
from matplotlib import gridspec
import plotly.graph_objects as go
import numpy as np
import pandas as pd



def calcPooledStdDev(df, tag1=None, tag2=None, ns=None):
    """Calculates the pooled standard deviations for two differently tagged groups in the same data frame"""
    
    if((not tag1) or (not tag2)):
       print("Need two tags to continue")
       return(-1)

    if(not ns):
       ns = [1, 1]
    
    entries = df.index.unique().to_numpy()
    stddevs = np.empty(entries.shape[0])
    
    for (i, entry) in enumerate(entries):
       stddevs[i] = np.sqrt((ns[0] * df[df['Tag']==tag1].loc[entry, "ColStdDev"]**2 + ns[0] * df[df['Tag']==tag2].loc[entry,"ColStdDev"]**2)/(ns[0]+ns[1]))
       

    return(pd.DataFrame({"Pooled Std Dev":stddevs}, index=entries))


def getStdDevsDF(df, tags=None, direction="column"):
    """Takes the dataframe of binding sites by cell and motif and creates a dataframe of standard deviations, either by row or column"""

    # Default action with tags: get all the unique items from the Tag column
    if(not tags):
        tags = df['Tag'].unique().tolist()
    


    # Check the direction of the standard deviations
    if (direction == 'row'):
        # This is probably not fully tested. But I'm also less sure of
        # the advantages of doing "stats by cell over all motifs", TBH
        nrows = df.shape[0] # Number of rows in the final df
        stddevs = pd.DataFrame([[tags[i]].std(axis=1).array,df['Tag']], columns = ['RowStdDev', tags], index=range(1,nrows+1))
    elif (direction == 'column'):
        nTags = len(tags)
        nrows = nTags * (df.shape[1]) # Number of rows in the final df; times n because of the tags
        tagsArray = np.repeat(tags, (df.shape[1]-2))
        stddevsList = []
        # Get the std devs 
        for tag in tags:
            stddevsList.append(df[df['Tag'] == tag].std(axis=0,numeric_only=True).array)


        stddevsList2 = []
        for i in range(len(stddevsList)):
            for j in range(len(stddevsList[0])):
                stddevsList2.append(stddevsList[i][j])

        stddevsArray = np.array(stddevsList2)

        finalArray=np.transpose(np.squeeze(np.stack((stddevsArray,tagsArray), axis=0)))

        stddevs = pd.DataFrame(finalArray, columns = ['ColStdDev', 'Tag'], index=np.tile(df.columns[1:-1], len(tags)))
        stddevs = stddevs.astype({"ColStdDev":float, "Tag":str})
        
    else:
        print(f"Direction passed ({direction}) is not valid. Use either 'row' or 'column'.")
        return(-1)


    return(stddevs)    



def getMediansDF(df, tags=None, direction='column'):
    """Takes the dataframe of binding sites by cell and motif and creates a dataframe of medians, either by row or column"""

    # Default action with tags: get all the unique items from the Tag column
    if(not tags):
        tags = df['Tag'].unique().tolist()
    


    # Check the direction of the averages
    if (direction == 'row'):
        # This is probably not fully tested. But I'm also less sure of
        # the advantages of doing "averages by cell over motifs", TBH
        nrows = df.shape[0] # Number of rows in the final df
        medians = pd.DataFrame([[tags[i]].median(axis=1).array,df['Tag']], columns = ['RowAve', tags], index=range(1,nrows+1))
    elif (direction == 'column'):
        nTags = len(tags)
        nrows = nTags * (df.shape[1]) # Number of rows in the final df; times n because of the tags
        tagsArray = np.repeat(tags, (df.shape[1]-2))
        mediansList = []
        # Get the medians 
        for tag in tags:
            mediansList.append(df[df['Tag'] == tag].median(axis=0,numeric_only=True).array)


        mediansList2 = []
        for i in range(len(mediansList)):
            for j in range(len(mediansList[0])):
                mediansList2.append(mediansList[i][j])

        mediansArray = np.array(mediansList2)

        finalArray=np.transpose(np.squeeze(np.stack((mediansArray,tagsArray), axis=0)))

        medians = pd.DataFrame(finalArray, columns = ['ColAve', 'Tag'], index=np.tile(df.columns[1:-1], len(tags)))
        medians = medians.astype({"ColAve":float, "Tag":str})
        
    else:
        print(f"Direction passed ({direction}) is not valid. Use either 'row' or 'column'.")
        return(-1)


    return(medians)
    




def getMeansDF(df, tags=None, direction='column'):
    """Takes the dataframe of binding sites by cell and motif and creates a dataframe of means, either by row or column"""

    # Default action with tags: get all the unique items from the Tag column
    if(not tags):
        tags = df['Tag'].unique().tolist()
    

    # Check the direction of the averages
    if (direction == 'row'):
        # This is probably not fully tested. But I'm also less sure of
        # the advantages of doing "averages by cell over motifs", TBH
        nrows = df.shape[0] # Number of rows in the final df
        means = pd.DataFrame([[tags[i]].mean(axis=1).array,df['Tag']], columns = ['RowAve', tags], index=range(1,nrows+1))
    elif (direction == 'column'):
        nTags = len(tags)
        nrows = nTags * (df.shape[1]) # Number of rows in the final df; times n because of the tags
        tagsArray = np.repeat(tags, (df.shape[1]-2))
        meansList = []
        # Get the means 
        for tag in tags:
            meansList.append(df[df['Tag'] == tag].mean(axis=0,numeric_only=True).array)

        meansList2 = []
        for i in range(len(meansList)):
            for j in range(len(meansList[0])):
                meansList2.append(meansList[i][j])

        meansArray = np.array(meansList2)

        finalArray=np.transpose(np.squeeze(np.stack((meansArray,tagsArray), axis=0)))

        means = pd.DataFrame(finalArray, columns = ['ColAve', 'Tag'], index=np.tile(df.columns[1:-1], len(tags)))
        means = means.astype({"ColAve":float, "Tag":str})
        
    else:
        print(f"Direction passed ({direction}) is not valid. Use either 'row' or 'column'.")
        return(-1)


    return(means)





def computeDifferentiability(df1, df2, motif, nBins=200, nSigma=5, maxX=None):
    """ Plots what amounts to 1-CDF of distributions. Basically, the fraction of cases more than a certain value."""

    data1 = df1[motif].to_numpy()
    data2 = df2[motif].to_numpy()

    # Drop NANS
    data1 = data1[~np.isnan(data1)]
    data2 = data2[~np.isnan(data2)]
    if(not maxX):
        maxX = max([np.mean(data1)+ nSigma * np.std(data1), np.mean(data2)+ nSigma * np.std(data2)])

    # maxX = np.max(np.concatenate((data1, data2), axis=None))
    xs = np.linspace(0, maxX, nBins)

    # Get the KDE functions
    kde1 = spstats.gaussian_kde(data1)
    kde2 = spstats.gaussian_kde(data2)
    

    # Generate the KDE values
    kdeVals1 = kde1(xs)
    kdeVals2 = kde2(xs)

    kdeVals1 /= np.sum(kdeVals1)
    kdeVals2 /= np.sum(kdeVals2)
    
    # Create empty arrays for the "effectiveness" metrics, aka, 1- CDF
    effs1 = np.empty(nBins)
    effs2 = np.empty(nBins)
    differentiability = np.empty(nBins)

    # The last element of these is, by defintion, the last KDE value
    effs1[-1] = kdeVals1[-1]
    effs2[-1] = kdeVals2[-1]

    for i in range(2, nBins+1):

        effs1[-i] = effs1[-i+1]+kdeVals1[-i]
        effs2[-i] = effs2[-i+1]+kdeVals2[-i]

    differentiability = effs2-effs1

    
    return(xs, effs1, effs2, differentiability)


def computePeakDiffStat(df1, df2, motif, nBins=200, nSigma=5, abs=False, both=False):
    """ Computes the peak of the difference in the differentiability/inverse-CDF curves of two different distributions of a given motif."""
    
    (xs, effs1, effs2, diffs) = computeDifferentiability(df1, df2, motif, nBins, nSigma)

    maxAbs = np.max(np.abs(diffs))
    maxLoc = np.argmax(np.abs(diffs))
    maxSigned = diffs[maxLoc]
    if(both):
        return((maxSigned, maxAbs))
    elif(abs):
        return(maxAbs)
    else:
        return(maxSigned)

    

def computeCohensEffectSize(df1, df2, motif):
    """ Computes Cohen's Effect Size statistics for two dataframes for a given motif. """

    data1 = df1[motif].to_numpy()
    data2 = df2[motif].to_numpy()

    # Drop NANS
    data1 = data1[~np.isnan(data1)]
    data2 = data2[~np.isnan(data2)]
    
    n1 = data1.shape[0]
    n2 = data2.shape[0]

    
    mean1 = np.mean(data1)
    mean2 = np.mean(data2)

    stddev1 = np.var(data1)
    stddev2 = np.var(data2)

    s = np.sqrt(((n1 - 1) * stddev1 + (n2-1) * stddev2)/(n1+n2-2))

    
    eff = abs(mean1-mean2)/s

    if(eff >= 0.8):
        description = "Large"
    elif(eff >=0.5):
        description = "Medium"
    elif(eff >= 0.2):
        description = "Small"
    else:
        description = "None"
    
    return((eff, description))

def compareMotifMediansBetweenSamples(df1, df2, motifs):
     """Compares the medians of the distribution of transcription factors
     in all the single-cell assays between two samples. Pass in
     dataframes with transcription factor counts by cell and motif and
     a list of motifs to check for."""

     if(type(motifs) != list):
         motifs=[motifs]

     nMotifs = len(motifs)

     # Create empty dictionaries to hold results
     stats = {}
     pvalues = {}

     
     for motif in motifs:
         type1 = df1[motif]
         type2 = df2[motif]

         (stat, pvalue) = spstats.mannwhitneyu(type1, type2, alternative='two-sided')

         stats.update({motif:stat})
         pvalues.update({motif:pvalue})


     return(stats, pvalues)

def motifsStatsDF(df1, df2, sortBy="Peak Differentiability Abs", nTrailingColumns = 1):
    motifs = df1.columns.to_list()[1:-(1+nTrailingColumns)]
    (stats, pvalues) = compareMotifMediansBetweenSamples(df1, df2, motifs)

    motifs.sort()
    
    pvals = np.empty(len(motifs))
    peakDiffs = np.empty(len(motifs))
    peakDiffsAbs = np.empty(len(motifs))
    effs = np.empty(len(motifs))
    descriptors = np.empty(len(motifs), dtype="U10")
    ratios = np.empty(len(motifs))
    
    
    for i, motifOfInterest in enumerate(motifs):
        (peakDiffs[i], peakDiffsAbs[i]) = computePeakDiffStat(df1,df2, motifOfInterest, both=True)
        pvals[i] = pvalues[motifOfInterest]
        (effs[i], descriptors[i])  = computeCohensEffectSize(df1, df2, motifOfInterest)
        if(pvals[i] != 0):
            ratios[i] = peakDiffs[i]/pvals[i]
        else:
            ratios[i] = np.NaN


    dfMotifs = pd.DataFrame({"p-value": pvals, "Peak Differentiability":peakDiffs, "Peak Differentiability Abs":peakDiffsAbs, "Eff Size":effs, "Eff Scale":descriptors}, index=motifs)
    dfMotifs=dfMotifs.sort_values(by=[sortBy], ascending=False)
    dfMotifs.index.name = "TF"
    
    return(dfMotifs)

