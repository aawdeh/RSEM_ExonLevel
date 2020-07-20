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

####################################################################################################################################
# Prepare exon interval file "exonFile"
# Using the gtf file and the transcript rsem bam file generated by RSEM for each RNA experiment, we will generate the weighted
# exon counts -- "expected counts per exon"
####################################################################################################################################
def prepareGTF():
    """
        Using the gtf file, we create a file of exons with their corresponding gene ids, transcript ids and exon ids.
        The chr, start and end of each exon corresponds to the genomic mapping.
        However, in this case, we are interested in the transcriptome mapping. -- So we need to convert genomic mapping to transcriptome mapping.

        Takes forever to run -- need to make more efficient. But only ran it once to create the exon intervals so not top priority.

    """
    df = pd.read_csv("/project/6019283/aaseel/Aim2/Data/Metadata/intervalsExon_fullDetails.txt", sep="\t")
    df['diff']=df.end-df.start
    groups = df.sort_values(['exon_number']).groupby(["transcript_id"])
    df_mod = pd.DataFrame()

    for method, group in groups:
        print(method)
        group.reset_index(inplace=True)
        del group['index']

        start = []
        end = []
        for index,row in group.iterrows():
            if index == 0:
                start = [1]
                end = [group['diff'][0] + 1]
            else:
                start.append(end[index-1])
                end.append(start[index] + row['diff'])

        group.insert(2, "start_transcript", start, True)
        group.insert(4, "end_transcript", end, True)

        df_mod = pd.concat([df_mod, group], axis=0)

    df_mod.to_csv("/project/6019283/aaseel/Aim2/Data/Metadata/intervalsExon_fullDetails_new.txt", sep="\t", index=False)  
