import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import seaborn as sns
import pandas as pd
import scipy.stats as spstats
from matplotlib import gridspec
import argparse




def bubblePlot(df, vertCol, sizeCol, clusterCol, labelCol=None, clusterNames=None, allowedLabels=None,
               colors = None, titles = None, alpha=0.5, vertLabel="", sizeLabel="", sizeVals = None,
               topTitle="",
               minSize=10, maxSize=500, vLim=None, sizeFormat=".0f", labelFontSize=None,
               largeLabels=None, largeFontSize=None, evenSpacing=False, linearSpacing=False, plotSize=None, 
               font={'family':'sans-serif','sans-serif':['Helvetica']}):
    """Makes a 'bubble plot' with rising arc points (it's cool) and points
clustered by some variable. Pass parameters are the dataframe with the
data (df), the title of the column used for the vertical (vertCol),
the title of the column used to set the bubble sizes (sizecol). If
labelCol is set, points will be labeled with the data in the column
with that title -- unless an associative array is also passed via
allowedLabels, in which case only labels in that list are
printed. (The associative array is structure with the keys being the
cluster names and the data being lists of the allowed labels.) If
clusterNames (an array or list) is passed, that will set the clusters
in clusterCol used to plot; if no clusterNames is passed, the program
will automatically use all of the unique items in clusterCol. Colors
sets the plotting colors (a list with the same number of items as
clusters to be plotted), titles sets the cluster titles for labeling
on the x-axis, alpha sets the opacity, vertLabel sets the label on the
y-axis, sizeLabel sets the label in the legend that gives the key to
the bubble sizes, sizeVals sets the values of the sizes in the legend
(if not set, quartiles are used), topTitle is the overall title of the
figure. minSize and maxSize are the min and max bubble sizes (in
plotting units, not data units), and vLim (a two-part tuple) sets the
vertical limits.

    """


    if(not (largeLabels is None)):
        if(largeFontSize is None):
            if(labelFontSize is None):
                largeFontSize = plt.rcParams['font.size'] * 2
            else:
                largeFontSize = labelFontSize * 2

    # Font
    #plt.rc('font', **font)
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['figure.figsize'] = 5, 6
            

    dfCopy = df
    if(not clusterNames):
        clusterNames = dfCopy[clusterCol].unique()
        
    nClusters = len(clusterNames)
    
    
    dfCopy['Horiz'] = np.zeros(dfCopy.shape[0])

    for (i, cluster) in enumerate(clusterNames):
        subDF = dfCopy[dfCopy[clusterCol]==cluster]
        if(evenSpacing):
            horizArray =subDF[vertCol].rank()
        elif(linearSpacing):
            horizArray = subDF[vertCol]
        else:
            horizArray = np.sqrt(subDF[vertCol]-np.min(subDF[vertCol]))

        width =  np.max(horizArray) - np.min(horizArray)
        if(width == 0):  width = 1.0
        horizArray = (horizArray - np.min(horizArray))/width + i

        for (index, row) in (dfCopy[dfCopy[clusterCol]==cluster]).iterrows():
            dfCopy.loc[index, 'Horiz'] = horizArray[index]


    # Plot
    fig = plt.figure(figsize=plotSize)
    ax=sns.scatterplot(x="Horiz", y=vertCol, size=sizeCol, data=dfCopy, sizes=(minSize, maxSize), 
                    hue=clusterCol, palette=colors, alpha=alpha, legend=None)

      
    # Plot limits. The maxRadius part makes sure that all of the bubbles are in-frame.
    oneUnit = (ax.transData.transform((1, 1))-ax.transData.transform((0,0))) * 72.0/fig.dpi
    maxRadiusX = np.sqrt(maxSize/np.pi)/oneUnit[0]
    maxRadiusY = np.sqrt(maxSize/np.pi)/oneUnit[1]
    if(vLim):
        plt.ylim(vLim)
    else:
        plt.ylim(np.min(dfCopy[vertCol])-maxRadiusY, np.max(dfCopy[vertCol])+maxRadiusY)

    plt.xlim(-maxRadiusX, nClusters+maxRadiusX)


    
    plt.xticks(np.arange(nClusters)+0.5, titles)
    
    plt.title(topTitle)
    plt.xlabel("")
    plt.ylabel(vertLabel)

    # Make the size legend
    labelPlts = []
    dummyLegends = []
    
    formatString = "{:"+sizeFormat+"}"
    if(sizeVals):
        for sz in sizeVals:
            labelPlts.append(formatString.format(sz))
            szPts = sz * maxSize/np.max(dfCopy[sizeCol])
            dummyLegends.append(plt.scatter([],[], s=szPts, edgecolors='black', facecolors='none'))
     
    else:
        quartilesPtSize = np.zeros(nClusters)
        quartilesColSize = np.zeros(nClusters)
        for i in range(3):
            quartilesPtSize[i] = (maxSize-minSize)/4.0 * (i+1) + minSize
            quartilesColSize[i] = (np.max(dfCopy[sizeCol])-np.min(dfCopy[sizeCol]))/4.0 * (i+1) + np.min(dfCopy[sizeCol])
            labelPlts.append(formatString.format(quartilesColSize[i]))
            dummyLegends.append(plt.scatter([],[], s=quartilesPtSize[i], edgecolors='black', facecolors='none'))
        
        
    leg = plt.legend(dummyLegends, labelPlts, ncol=3, frameon=True,  title=sizeLabel, scatterpoints = 1)

    if(labelCol):
        for cluster in clusterNames:
            subDF = dfCopy[dfCopy[clusterCol]==cluster]
            for (i, row) in subDF.iterrows():

                if((allowedLabels is None) or (row[labelCol] in allowedLabels[cluster])):
                    if((largeLabels is None) or (not (row[labelCol] in largeLabels))):
                        plt.annotate(row[labelCol], (row['Horiz'], row[vertCol]), xytext=(-10,0),
                                     textcoords='offset points', ha='right', fontsize=labelFontSize)
                    else:
                        plt.annotate(row[labelCol], (row['Horiz'], row[vertCol]), xytext=(-10,0),
                                     textcoords='offset points', ha='right', fontsize=largeFontSize, weight='bold')
    plt.tight_layout()
    return(fig)






                        
def bubblePlotFunc(vert, sz, labels=None, colors = ['blue', 'orange', 'green'], titles = ["Cluster 1", "Cluster 2", "Cluster 3"], labelList=None, scale=1, alpha=0.5, vertLabel="", sizeLabel="", topTitle="", font={'family':'sans-serif','sans-serif':['Helvetica']}):

    """Plotting routine that does a 'bubble plot'. Takes in arrays of the
vertical data and the bubble sizes. Arrays should be 2D, nxm where n
is the number of clusters of data nadm is the number of points in each
set. The bubbles are plotted with a rising arc, purely because it
looks cool. Labels for the points are passed via labels (also nxm). Colors set the colors of the columns. Titles gives the titles of the clusters."""

    # Read the inputs
    #
    
    
    wDim = vert.shape[0]
    horz = np.zeros((vert.shape[0], vert.shape[1]))

    # Font
    plt.rc('font', **font)

    for i in range(wDim):
        horz[i, :] = np.sqrt((vert[i,:]-np.min(vert[i,:])))
        
        width =  (np.max(horz[i,:]) + np.min(horz[i,:]))
        
        horz[i, :] = (horz[i]-np.min(horz[i,:]))/width + i
        
        
        
        ax=plt.scatter(horz[i,:] , vert[i,:], s=sz[i,:]*scale, c=colors[i], alpha=alpha)
        

#    if(args.vmax):
#        plt.ylim(np.min(vert), args.vmax)
#    else:
#        plt.ylim(np.min(vert), np.max(vert))

    plt.xticks([0.5, 1.5, 2.5], titles)

    plt.title(topTitle)
    plt.xlabel("")
    plt.ylabel(vertLabel)

    # Make the legend

    quartiles = np.zeros(wDim)
    labelPlts = []
    dummyLegends = []
    for i in range(3):
        quartiles[i] = np.quantile(sz, 0.25 * (i+1))
        labelPlts.append("{:.1E}".format(quartiles[i]))
        dummyLegends.append(plt.scatter([],[], s=quartiles[i]*scale, edgecolors='black', facecolors='none'))


    leg = plt.legend(dummyLegends, labelPlts, ncol=3, frameon=True,  title=sizeLabel, scatterpoints = 1)


    if(not (labels is None)):
        for i in range(wDim):
            for j, txt in enumerate(labels[i,:]):

                if((labelList is None) or (txt in labeList[i, :])):
                   plt.annotate(txt, (horz[i,j], vert[i,j]), xytext=(-20,0), textcoords='offset points', ha='right')
