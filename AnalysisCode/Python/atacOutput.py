import argparse
import scipy.stats as spstats
from matplotlib import gridspec
import plotly.graph_objects as go

import atacIn
import atacStats

def printCSV(df1, df2, filename, sortBy=None):


    dfMotifs=motifsStatsDF(df1, df2)

    if(sortBy):
        dfMotifs=dfMotifs.sort_values(by=[sortBy], ascending=False)
    else:
        dfMotifs.sort_index(inplace=True)
        
    with open(filename, 'w') as f:
        f.write(dfMotifs.to_string( float_format='%.3G'))



def writeBED(df, filename, columns=3):

    df.to_csv(filename, sep="\t", index=False, columns=columns)
