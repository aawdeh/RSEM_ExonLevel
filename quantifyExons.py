# ------------------------------------
# python modules
# ------------------------------------
from __future__ import print_function
import os
import sys
import math
import numpy as np
import pandas as pd
import os.path
from os import path
import shutil
import sys

def weightedCount(fileName, outputDir):
    """
    """
    pairname = fileName.split("/")
    pairname = pairname[len(pairname) - 1].replace(".Exon.txt","")
    chunksize = 10 ** 6
    list_dfs = []

    for chunk in pd.read_csv(fileName, sep="\t", header=None, chunksize=chunksize):
        chunk.columns = ['transcript_id', 'start','end','exon_id', 'transcript_id_file', 'start_file', 'end_file','posterior_prob']
        chunk['posterior_prob'] = chunk['posterior_prob'].str.replace("ZW:f:","").astype(float)
        chunk['new_start'] = np.select([chunk.start_file < chunk.start], [chunk.start], default=chunk.start_file)
        chunk['new_end'] = np.select([chunk.end_file > chunk.end], [chunk.end], default=chunk.end_file)
        chunk['new_len'] = chunk.new_end - chunk.new_start
        chunk['length'] = chunk.end_file - chunk.start_file
        chunk['new_prob'] = chunk.posterior_prob * (chunk.new_len / chunk.length)
        chunk =  chunk.groupby(['transcript_id','exon_id'])['new_prob'].sum().to_frame()
        list_dfs.append(chunk)

    df = pd.concat(list_dfs)
    df =  df.groupby(['transcript_id','exon_id'])['new_prob'].sum().to_frame()
    print(df)
    df.to_csv(outputDir + pairname + ".txt", sep="\t")

if __name__== "__main__":
    fileName = sys.argv[1]
    outputDir = "/scratch/aaseel/Aim2/Data/RNAseq/Exon_ReadCounts/"
    print(fileName)
    weightedCount(fileName, outputDir)
