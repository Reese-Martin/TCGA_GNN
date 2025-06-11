# the gdc-client downloads the requested data in a series of separate folders. This script takes the files that were
# generated and then puts them together in a single csv

import pandas as pd
import numpy as np
import os
files = os.listdir("star_counts")
files = [i for i in files if '.' not in i]
metadata = pd.read_csv('Pt_metadata.csv')
# trim metadata file to the files that were successfully downloaded (ideally no files errored, but this is a good catch)
trim_metadata = metadata.loc[[True if i in files else False for i in metadata['file_id']], :]
# names of the folders storing data
ids = trim_metadata['file_id'].tolist()
# names of the individual tsv files
names = trim_metadata['file_name'].tolist()

# tpm_df = pd.DataFrame(index=ids, columns=['pct_Assigned', 'pct_noFeat', 'pct_unmapped'])
tpm_dfs = []
fpkm_dfs = []

for (fid, name) in zip(ids, names):
    tmp = pd.read_csv('star_counts/'+fid+'/'+name, sep='\t', comment='#')
    tmpReads = tmp.loc[4:]
    tmpQC = pd.DataFrame(data=tmp.loc[0:3, ['unstranded', 'stranded_first', 'stranded_second']].values,
                         index=tmp.loc[0:3]['gene_id'].tolist(),
                         columns=['unstranded', 'stranded_first', 'stranded_second'])
    N_unmapped = tmpQC.loc['N_unmapped', 'unstranded']
    N_ambig = tmpQC.loc['N_ambiguous', 'unstranded']
    N_noFeat = tmpQC.loc['N_noFeature', 'unstranded']

    assigned_reads = sum(tmpReads['unstranded'])
    total_reads = assigned_reads + N_unmapped
    total_mapped = assigned_reads + N_ambig + N_noFeat
    # for first pass QC I will focus on the percent assigned, percent unmapped, and the percent with no feature
    # if  pct_Assigned < .8, pct_unmapped < 15%, or the pct_noFeat < 20% the sample is rejected. These are relatively
    # stringent requirements.
    pct_Assigned = assigned_reads / total_mapped
    pct_unmapped = N_unmapped / total_reads
    pct_noFeat = N_noFeat / total_mapped

    if (pct_Assigned >= .8) and (pct_unmapped <= .15) and (pct_noFeat <= .2):
        tpm_tmp = pd.DataFrame(columns=[trim_metadata[trim_metadata['file_id'] == fid]['sample_barcode'].values],
                           index=tmpReads['gene_id'].values, data=tmpReads['tpm_unstranded'].values)
        fpkm_tmp = pd.DataFrame(columns=[trim_metadata[trim_metadata['file_id'] == fid]['sample_barcode'].values],
                               index=tmpReads['gene_id'].values, data=tmpReads['tpm_unstranded'].values)
        tpm_dfs.append(tpm_tmp)
        fpkm_dfs.append(fpkm_tmp)
    else:
        trim_metadata = trim_metadata[trim_metadata['file_id'] != fid]

tpmDF = pd.concat(tpm_dfs, axis=1)