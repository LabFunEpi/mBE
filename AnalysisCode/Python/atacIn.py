#
# Functions to read in ATAC-seq data
#
import os
import numpy as np
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sqlite3
import scipy.stats as spstats
import h5py
import pybedtools
import re




def readPeaksH5(infile, tag=None):
    """Read a Peaks file from a HDF5 file. Tag keyword is
optional, sets the tag column. This version (unlike the text files)
also sets the barcodes since those data are in the HDF5 file."""


    # Read the file
    hf = h5py.File(infile, 'r')

    # Read in the data we need
    shape = np.array(hf['matrix']['shape'])
    indices = np.array(hf['matrix'].get('indices'))
    indPtrs= np.array(hf['matrix'].get('indptr'))
    ndata=np.array(hf['matrix']['data'])

    # This assigns the data into the matrix
    dataArr = np.zeros((shape[1], shape[0]))
    for i in range(shape[1]):
        dataArr[i, indices[indPtrs[i]:indPtrs[i+1]]] = ndata[indPtrs[i]:indPtrs[i+1]].astype(int)

    # Get the indices of the locations of the non-zero elements of the above matrix
    (cells, peaks) = np.nonzero(dataArr)

    # Get the locations.
    chrom = []
    start=[]
    end = []
    for feat in np.array(hf['matrix']['features']['id']):
        chrloc = re.split(":|-", feat.decode())
        chrom.append(chrloc[0])
        start.append(chrloc[1])
        end.append(chrloc[2])

    chromosomeArr = np.array(chrom)
    startArr = np.array(start, dtype=np.int)
    endArr = np.array(end, dtype=np.int)

    chromosomes = chromosomeArr[peaks]
    startPts = startArr[peaks]
    endPts = endArr[peaks]
    

    # Get the barcodes
    bcodes = []
    for bcode in np.array(hf['matrix']['barcodes']):
        bcodes.append(bcode.decode())

    barcodeArr = np.array(bcodes)

    barcodes = barcodeArr[cells]

    
    cutSites = dataArr[cells, peaks].astype(int)

    
    # This is more like the right thing
    df = pd.DataFrame({"Cell":cells, "Chromosome":chromosomes, "Start":startPts, "End":endPts, "Barcode":barcodes, "CutSites":cutSites})

    if(tag is not None):
        df['Tag']=tag
    
    
    df.sort_values(by=["Cell", "Chromosome", "Start"], inplace=True)



    
    return(df)


def readPeakMatrixH5(infile, tag=None):
    """Read a Peaks Matrix from a HDF5 file. Tag keyword is
optional, sets the tag column. This version (unlike the text files)
also sets the barcodes since those data are in the HDF5 file."""


    # Read the file
    hf = h5py.File(infile, 'r')

    # Read in the data we need
    shape = np.array(hf['matrix']['shape'])
    indices = np.array(hf['matrix'].get('indices'))
    indPtrs= np.array(hf['matrix'].get('indptr'))
    ndata=np.array(hf['matrix']['data'])

    # This assigns the data into the matrix
    dataArr = np.zeros((shape[1], shape[0]))
    for i in range(shape[1]):
        dataArr[i, indices[indPtrs[i]:indPtrs[i+1]]] = ndata[indPtrs[i]:indPtrs[i+1]].astype(int)
        
    featID =  cleanH5ByteArray(hf['matrix']['features']['name'])

    # Create the dataframe
    df = pd.DataFrame(dataArr,  columns=featID,  dtype="int")

    df['Cell']=np.arange(1, df.shape[0]+1)
    
    # If the tag keyword is set, add a Tag column
    if(tag):
        df['Tag']=tag

    # Add barcodes
    barcodes = cleanH5ByteArray(hf['matrix']['barcodes'])
    df['Barcode']=barcodes


    
    return(df)



def readTFMatrixH5(infile, tag=None):
    """Read a Transcription Factor Matrix from a HDF5 file. Tag keyword is
optional, sets the tag column. This version (unlike the text files)
also sets the barcodes since those data are in the HDF5 file."""


    # Read the file
    hf = h5py.File(infile, 'r')

    # Feature Type
    if((cleanH5ByteArray(hf['matrix']['features']['feature_type']))[0] !="Motifs"):
        print("Feature type isn't motifs. Exiting.")
        return(-1)
    
    # Read in the data we need
    shape = np.array(hf['matrix']['shape'])
    indices = np.array(hf['matrix'].get('indices'))
    indPtrs= np.array(hf['matrix'].get('indptr'))
    ndata=np.array(hf['matrix']['data'])

    # This assigns the data into the matrix
    dataArr = np.zeros((shape[1], shape[0]))
    for i in range(shape[1]):
        dataArr[i, indices[indPtrs[i]:indPtrs[i+1]]] = ndata[indPtrs[i]:indPtrs[i+1]].astype(int)
        

    featID =  cleanH5ByteArray(hf['matrix']['features']['name'])

    # Create the dataframe and label the index "Cell"
    df = pd.DataFrame(dataArr, columns=featID, dtype="int")

    df['Cell']=np.arange(1, df.shape[0]+1)
    
    # If the tag keyword is set, add a Tag column
    if(tag):
        df['Tag']=tag

    # Add barcodes
    barcodes = cleanH5ByteArray(hf['matrix']['barcodes'])


    df['Barcode']=barcodes

    return(df)


def readTFAuxH5(infile):
    """ Reads the auxiliary information about the motifs from an HDF5 file """

    
    # Read the file
    hf = h5py.File(infile, 'r')

    # Motif name
    tfname = cleanH5ByteArray(hf['matrix']['features']['name'])

    # Full name
    tfFullName = cleanH5ByteArray(hf['matrix']['features']['id'])

    # Feature Type
    featureType =  cleanH5ByteArray(hf['matrix']['features']['feature_type'])
    
    # Derivation
    derivation = cleanH5ByteArray(hf['matrix']['features']['derivation'])
    
    # Genome
    genome = cleanH5ByteArray(hf['matrix']['features']['genome'])

    return(pd.DataFrame({"Motif":tfname, "Full Motif Name":tfFullName, "Derivation":derivation, "Genome":genome, "Type":featureType}))



    
def cleanH5ByteArray(h5Array, datatype="U30"):
    """Takes in a reference from an HDF5 file's byte array that we want
    to clean up into a character array. Optional argument is datatype,
    defaults to U30"""

    bytArr = np.array(h5Array)
    newArr = np.empty(bytArr.shape[0], dtype=datatype)
    for (i, elem) in enumerate(bytArr):
        newArr[i] = elem.decode()

    return(newArr)


def integrateOneBed(df, bed, columName="State"):

    """ Integrate four dataframes by looking for where peaks overlap with the active/inactive/mBE regions.
    Inputs are the peaks df, the active, inactive, and mBE dataframes.
    We'll actually use a SQL query to do this since it's WAY (~25 times) faster than doing it as a loop."""

    #Make the db in memory
    conn = sqlite3.connect(':memory:')
    #write the tables
    df.to_sql('dfDB', conn, index=False)
    bed.to_sql('bedDB', conn, index=False)



    # The query. We'll do three left-joins and use the "between" syntax (as well as matching the chromosome, of course) to do this
    qry = '''
    SELECT  
        dfDB.*,
        bedDB.'''+columName+'''
    FROM
        dfDB LEFT JOIN bedDB ON
                bedDB.Chromosome=dfDB.Chromosome AND
                    (bedDB.Start BETWEEN dfDB.Start AND dfDB.End OR
                    dfDB.Start BETWEEN bedDB.Start AND bedDB.End OR
                    bedDB.End BETWEEN dfDB.Start AND dfDB.End)
    ORDER BY
        dfDB.Chromosome, dfDB.Start, dfDB.End
    '''

    # Execute the query and return the df it produces
    df = pd.read_sql_query(qry, conn)
    df.drop_duplicates(inplace=True)
    return(df)



def integrateActiveInactivemBE(df, active, inactive, mBE):
    """ Integrate four dataframes by looking for where peaks overlap with the active/inactive/mBE regions.
    Inputs are the peaks df, the active, inactive, and mBE dataframes.
    We'll actually use a SQL query to do this since it's WAY (~25 times) faster than doing it as a loop."""

    #Make the db in memory
    conn = sqlite3.connect(':memory:')
    #write the tables
    df.to_sql('dfDB', conn, index=False)
    active.to_sql('activeDB', conn, index=False)
    inactive.to_sql('inactiveDB', conn, index=False)
    mBE.to_sql('mBEDB', conn, index=False)



    # The query. We'll do three left-joins and use the "between" syntax (as well as matching the chromosome, of course) to do this
    qry = '''
    SELECT  
        dfDB.*,
        activeDB.Active,
        inactiveDB.Inactive,
        mBEDB.mBE
    FROM
        dfDB LEFT JOIN activeDB ON
                activeDB.Chromosome=dfDB.Chromosome AND
                    (activeDB.Start BETWEEN dfDB.Start AND dfDB.End OR
                    activeDB.End BETWEEN dfDB.Start AND dfDB.End OR
                    dfDB.Start BETWEEN activeDB.Start AND activeDB.End)
            LEFT JOIN inactiveDB ON
                inactiveDB.Chromosome=dfDB.Chromosome AND
                    (inactiveDB.Start BETWEEN dfDB.Start AND dfDB.End OR
                    inactiveDB.End BETWEEN dfDB.Start AND dfDB.End OR
                    dfDB.Start BETWEEN inactiveDB.Start AND inactiveDB.End)
            LEFT JOIN mBEDB ON
                mBEDB.Chromosome=dfDB.Chromosome AND
                    (mBEDB.Start BETWEEN dfDB.Start AND dfDB.End OR
                    mBEDB.End BETWEEN dfDB.Start AND dfDB.End OR
                    dfDB.Start BETWEEN mBEDB.Start AND mBEDB.End)
    ORDER BY
        dfDB.Chromosome, dfDB.Start, dfDB.End
    '''

    # Execute the query and return the df it produces
    df = pd.read_sql_query(qry, conn)
    df.drop_duplicates(inplace=True)
    return(df)


def loadBED4(filename, state=None, stateName=None, column4Name="Column4"):
    """Read a slightly more complex, four-column BED file."""

    # Read the file itself
    data = np.loadtxt(filename, skiprows=1, dtype={'names': ('Chromosome', 'Start', 'End', 'Column4'),'formats': ('U30', 'int', 'int', 'U60')}, delimiter='\t')

    # Number of lines in the file. Important for later.
    nLines = data.shape[0]

    # Create empty arrays to hold the data read in above.
    chromosomes= np.empty(nLines, dtype="U5")
    startPts = np.zeros(nLines, dtype=np.int)
    endPts =np.zeros(nLines, dtype=np.int)
    column4s = np.empty(nLines, dtype="U60")
    
    # For each of the tuples returned by loadtxt (which returns an
    # array of tuples), unpack it and send each variable to the
    # appropriate array
    for ctr in range(nLines):
        (chromosomes[ctr], startPts[ctr], endPts[ctr], column4s[ctr]) = data[ctr]

    df = pd.DataFrame({"Chromosome":chromosomes, "Start":startPts, "End":endPts, column4Name:column4s})
    if(state):
        if(not stateName): stateName = "State"
        df[stateName] = state

    # Create and return the dataframe
    return(df)




def loadBED(filename, state="Unknown", useBedTools=False):
    """Read a simple, three-column BED file."""

    if(useBedTools):
        df = pybedtools.BedTool(filename).to_dataframe().rename(columns={"chrom":"Chromosome", "start":"Start", "end":"End"})
        df[state]=state
        return(df)
        
    else:
        # Read the file itself
        data = np.loadtxt(filename, skiprows=1, dtype={'names': ('Chromosome', 'Start', 'End'),'formats': ('U30', 'int', 'int')}, delimiter='\t')

        # Number of lines in the file. Important for later.
        nLines = data.shape[0]

        # Create empty arrays to hold the data read in above.
        chromosomes= np.empty(nLines, dtype="U5")
        startPts = np.zeros(nLines, dtype=np.int)
        endPts =np.zeros(nLines, dtype=np.int)

        # For each of the tuples returned by loadtxt (which returns an
        # array of tuples), unpack it and send each variable to the
        # appropriate array
        for ctr in range(nLines):
            (chromosomes[ctr], startPts[ctr], endPts[ctr]) = data[ctr]

        # Create and return the dataframe
        return(pd.DataFrame({"Chromosome":chromosomes, "Start":startPts, "End":endPts, state:np.repeat(state, nLines)}))




def loadActiveInactive(filename, state="Unknown"):
    """Read the active/inactive/mBE groupings"""

    # Read the file itself
    data = np.loadtxt(filename, skiprows=1, dtype={'names': ('Chromosome', 'Start', 'End', 'Name'),'formats': ('U30', 'int', 'int', 'U30')}, delimiter='\t')

    # Number of lines in the file. Important for later.
    nLines = data.shape[0]

    # Create empty arrays to hold the data read in above.
    chromosomes= np.empty(nLines, dtype="U5")
    startPts = np.zeros(nLines, dtype=np.int)
    endPts =np.zeros(nLines, dtype=np.int)

    # For each of the tuples returned by loadtxt (which returns an
    # array of tuples), unpack it and send each variable to the
    # appropriate array
    for ctr in range(nLines):
        (chromosomes[ctr], startPts[ctr], endPts[ctr], junk) = data[ctr]

    # Create and return the dataframe
    return(pd.DataFrame({"Chromosome":chromosomes, "Start":startPts, "End":endPts, state:np.repeat(state, nLines)}))
    


def loadClusterData(filename):
    """Read the clustering data. Returns a dictionary of barcode/cluster pairs."""

    # Read the file
    data = np.loadtxt(filename, skiprows=1, dtype={'names': ('Barcode', 'Cluster'),'formats': ('U30', 'int')}, delimiter=',')

    # Get the number of lines in the files
    nLines = data.shape[0]

    # Create empty lists to hold the file data
    barcodes = []
    clusters = []

    # Unpack the each tuple in the array of tuples returned by the loadtxt function and save the results to the appropriate list
    for ctr in range(nLines):
        # Convert the list of tuples into lists
        (bc, cl) = data[ctr]
        barcodes.append(bc)
        clusters.append(cl)

    # Return a dictionary of the barcodes and the clusters. (Former is the key, the later is the entry.)
    return(dict(zip(barcodes, clusters)))





def combineCluster(df, clustersDict):
    """Combine a dataframe with barcodes with the cluster. df has to be a
    dataframe with a "barcode" column, clusters is a Dict output from
    reading cluster data"""

    # Add a column of cluster data to the dataframe using the cluster dictionary and the barcode column
    df['Cluster']= df['Barcode'].map(clustersDict)
        
    return(df)
    



def loadMotifsBED(filename, compressRows = True):
    """ Read in a motif .bed file and turn it into a dataframe"""

    # Actually read in the file. Sadly, into an array of tuples.
    motifsList = np.loadtxt(filename,  dtype={'names': ('chroms', 'starts', 'fins', 'Motifs'),'formats': ('U5', 'int', 'int', 'U30')}, delimiter='\t')
    nLines = motifsList.shape[0] # Number of lines read in


    
    # Create empty arrays for holding the data before turing the whole thing into a dataframe
    chromosomes= np.empty(nLines, dtype="U5")
    startPts = np.zeros(nLines, dtype=np.int)
    endPts =np.zeros(nLines, dtype=np.int)
    motifs= np.empty(nLines, dtype = "U30")
    
    for ctr in range(nLines):
        (chromosomes[ctr], startPts[ctr], endPts[ctr], motifs[ctr]) = motifsList[ctr]

    # Make the dataframe
    df = pd.DataFrame({"Chromosome":chromosomes, "Start":startPts, "End":endPts, "Motif":motifs})
    print(len(df))
    # Now to remove redundancies
    df.drop_duplicates(inplace=True)
    print(len(df))
    

    if(compressRows):
        # We're going to compress the motifs from one peak into one line to make the dataframe more compact
        
        # Arrays to hold things
        chromosomes2 = []
        startPts2 = []
        endPts2 = []
        motifs2 = []
        
        
        for chr in df.Chromosome.unique():
            dfChr = df[df['Chromosome'] == chr]
            for strt in dfChr.Start.unique():
                dfStrt =  dfChr[dfChr['Start'] == strt]
                for end in dfStrt.End.unique():
                    mtfs = ",".join(dfStrt[dfStrt['End']==end].Motif.values.tolist())
                    chromosomes2.append(chr)
                    startPts2.append(strt)
                    endPts2.append(end)
                    motifs2.append(mtfs)

        df2 =pd.DataFrame({"Chromosome":chromosomes2, "Start":startPts2, "End":endPts2, "Motif":motifs2})
        print(len(df2))
        return(df2)
        
    else:
        # Not compressing the motif data into a single row for each location, so just return what we already have
        return(df)

        

    
def readBarcodes(infile, reverse=False, cellCharacter=False):
    """Read the barcodes and return a hash. If a directory is passed, the
routine guesses the filename in that directory. By default, the key is
the barcode and the data is the cell number. If reverse is set to
True, that's swapped."""

    # If the passed "infile" is a direcotry, add "barcodes.tsv" to the end to guess the filename.
    if(os.path.isdir(infile)):
        infile = os.path.join(infile, "barcodes.tsv")

    # Read the barcodes file.
    barcodesA = np.loadtxt(infile, dtype="U30")

    # Cells are from 1 to number of lines in barcodes file. They map
    # in the obvious way to the file lines.
    cell = np.arange(1, barcodesA.shape[0]+1)

    # If the cellCharacter flag is set, make the cell array a character array. Otherwise, it's ints
    if(cellCharacter):
        cell = np.char.mod("%d", cell)

    # Return a dictionary of cell number and barcodes.
    # If reverse is true, we'll return a dictionary with cell as the key. Otherwise, the bacode will be.
    if(reverse):
        return(dict(zip(cell, barcodesA)))
    else:
        return(dict(zip(barcodesA, cell)))


    
    
    
def readPeaks(peaksFile, tag='', barcodesFile=None, matrixFile=None):
    """ Parse a peaks.bed file into a dataframe"""

    # If the peaksFile passed is really a directory, then let's guess the filenames
    if(os.path.isdir(peaksFile)):
        # Guess barcode file
        barcodesFile = os.path.join(peaksFile, "barcodes.tsv")

        # Also guess the matrix files
        matrixFile= os.path.join(peaksFile, "matrix.mtx")
        # And add the peaks file to the directory path
        peaksFile = os.path.join(peaksFile, "peaks.bed")
        

    # Check to see if files exist
    if(not (barcodesFile and os.path.exists(barcodesFile))):
        barcodesFile = None
    if(not (matrixFile and os.path.exists(matrixFile))):
        matrixFile = None
    if(not (peaksFile and os.path.exists(peaksFile))):
        peaksFile = None

    # If any of the necessary files aren't there, quite. Nicely.
    if(not (peaksFile and barcodesFile and matrixFile)):
        print("Either not all input files were provided or some of the files didn't exist.")
        return(None)

    # Read the matrix file and save the output
    data = np.loadtxt(matrixFile, skiprows=3, dtype='int', delimiter=' ') # THis just limits to the loading time:
    lineInPeaks = data[:,0]
    lineInBarcodes = data[:,1]
    cutSites = data[:,2]
    nLines = cutSites.shape[0]

    # Read the peaks datafile
    peaksA = np.loadtxt(peaksFile,dtype={'names': ('chroms', 'strts', 'fins'),'formats': ('U5', 'int', 'int')} , delimiter='\t')


    # Read the barcodes file
    barcodesA = np.loadtxt(barcodesFile, dtype="U30")

    # Get the barcodes as they correspond to the matrix file.
    barcodes = np.empty(nLines, dtype = "U20")
    chromosomes=np.empty(nLines, dtype="U5")
    startPts = np.zeros(nLines, dtype=np.int)
    endPts =np.zeros(nLines, dtype=np.int)

    # Unpack the array of tuples of peaks data
    for ctr in range(nLines):
        (chromosomes[ctr], startPts[ctr], endPts[ctr]) = peaksA[lineInPeaks[ctr]-1]
        barcodes[ctr] = barcodesA[lineInBarcodes[ctr]-1]
        
    # Assemble it all into a dataframe with cells, chromsome, start, end, barcode, and number of cut sites
    df = pd.DataFrame({"Cell":lineInBarcodes-1, "Chromosome":chromosomes, "Start":startPts, "End":endPts, "Barcode":barcodes, "CutSites":cutSites})

    df.sort_values(by=["Cell", "Chromosome", "Start"], inplace=True)

    # If the tag parameter was passed, add a column of that tag in
    # every entry. (Useful in you combine dataframes laters.)
    if (tag != ''):
        df['Tag'] = np.repeat(tag, nLines)

    # And return it.
    return(df)



def readTFMatrix(matrixFile, tag='', motifsFile=None, fullMotifNames=False):
    """Read a matrix file of motifs and cells and turn it into a dataframe that can be more easily read and manipulated"""

    # if the matrixFile is a directory, we can probably work out what the correct filenames are.
    if(os.path.isdir(matrixFile)):
        motifsFile = os.path.join(matrixFile, "motifs.tsv")
        matrixFile = os.path.join(matrixFile, "matrix.mtx")

    if(not os.path.exists(motifsFile)):
        motifsFile=None
    if(not os.path.exists(matrixFile)):
        matrixFile = None

    if(not (matrixFile and motifsFile)):
        print("Motifs and matrix file (or a path to a directory containing them) don't both exist or weren't specified. Exiting.")
        return(None)

    
    # Read the data in
    data = np.loadtxt(matrixFile, skiprows=3, dtype='int', delimiter=' ')#, max_rows=(579-1)*10-1) # THis just limits to the loading time:

    motifs = data[:, 0]
    cells = data[:, 1]
    sites = data[:, 2]
    maxMotifs = max(motifs)
    maxCells = max(cells)

    sitesArr = np.full((maxCells, maxMotifs+1), 0)


    for i in range(sites.shape[0]):
        sitesArr[cells[i]-1,motifs[i]] = sites[i]


    # If we have motif data file information, read the file and use the names
    if(motifsFile):
        motifs = np.array2string(motifs) # Stupid, but need to make it strings now or else the assigment blows up

        motifData = np.loadtxt(motifsFile, dtype='S25')
        if(fullMotifNames):
            motifs=motifData[:,0] # Full names
        else:
            motifs = motifData[:,1] # Common names, default action
        
    else:
        # If we have no motif data file, just call 'em "Motif1", "Motif2", etc. It's something, at least.
        motifs = np.char.mod('Motif%d', range(1,maxMotifs+1))

        
    df = pd.DataFrame(sitesArr, columns=np.append(['Cell'], motifs))

    df['Cell']= np.char.mod("%d", np.arange(1,maxCells+1, dtype=np.int))

    
    if (tag != ''):
        df['Tag'] = np.repeat(tag, maxCells)

        
    return(df)


def mergeDFs(df1, df2):
    """Merge two dataframes. Keeps the 'tag' column intact"""

    # Pretty straight forward, to be honest.
    return(pd.concat([df1, df2]))



    

def splitDFbyTag(df, tags=None):
    """Takes a dataframe and splits the columns out into separate columns by the tags. Rows are matched by index"""

    
    # Default action with tags: get all the unique items from the Tag column
    if(not tags):
        tags = df['Tag'].unique().tolist()

    # Get the number of tags
    nTags = len(tags)

    # Get the indices
    indices = df.index.unique().tolist()
    nIndices = len(indices)

    # Get the column names, except "Tag". We don't need that anymore.
    cols = df.columns
    cols = cols[cols != 'Tag']
    nCols = cols.shape[0]


    # Temp array to hold the variables before they're popped into a dataframe
    newArray = np.full((nIndices, nTags*nCols), np.nan)
    newCols = np.empty((nTags*nCols), dtype='U25')#np.full((nTags*nCols), " ")

    # Loop over all the tags and pull up a sub-dataframe of just matches to that tag
    for (j, tag) in enumerate(tags):

        dfTemp = df[df['Tag']==tag]

        # Now loop over all the columns, in part to create the column list
        for (k, col) in enumerate(cols.tolist()):

            newCols[j+k] = col + " " + tag

            # And finally, loop over the indices. Populate the temporary array with the entries in order to eventually fill the new DF
            for (i, idx) in enumerate(indices):

                newArray[i, j+k] = dfTemp.loc[idx, col]


    # Create the dataframe from the column headers, data, and indices
    dfNew = pd.DataFrame(newArray, columns = newCols, index=indices)

    return(dfNew)

        


