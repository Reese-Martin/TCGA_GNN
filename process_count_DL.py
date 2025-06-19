# the gdc-client downloads the requested data in a series of separate folders. This script takes the files that were
# generated and then puts them together in a single csv

import pandas as pd
import os


metadata = pd.read_csv('Pt_metadata.csv')
# Get valid sample folders
files = [f for f in os.listdir("star_counts") if os.path.isdir(os.path.join("star_counts", f))]
valid_file_ids = set(files)

# trim metadata file to the files that were successfully downloaded (ideally no files errored, but this is a good catch)
trim_metadata = metadata[metadata['file_id'].isin(valid_file_ids)].copy()

# Precompute file_id → filename and barcode
file_id_to_filename = dict(zip(trim_metadata['file_id'], trim_metadata['file_name']))
file_id_to_barcode = dict(zip(trim_metadata['file_id'], trim_metadata['sample_barcode']))


# tpm_df = pd.DataFrame(index=ids, columns=['pct_Assigned', 'pct_noFeat', 'pct_unmapped'])
raw_dfs = []
tpm_dfs = []
fpkm_dfs = []

for fid in trim_metadata['file_id']:
    path = os.path.join("star_counts", fid, file_id_to_filename[fid])

    # Identify QC rows and gene expression rows
    df = pd.read_csv(path, sep='\t', comment='#')
    qc_rows = df[df['gene_id'].str.startswith('N_')]
    qc = qc_rows.set_index('gene_id')['unstranded'].astype(float)
    gene_rows = df[~df['gene_id'].str.startswith('N_')]

    N_unmapped = qc['N_unmapped']
    N_ambig = qc['N_ambiguous']
    N_noFeat = qc['N_noFeature']

    assigned_reads = gene_rows['unstranded'].sum()
    total_reads = assigned_reads + N_unmapped
    total_mapped = assigned_reads + N_ambig + N_noFeat

    # for first pass QC I will focus on the percent assigned, percent unmapped, and the percent with no feature
    # if  pct_Assigned < .8, pct_unmapped < 15%, or the pct_noFeat < 20% the sample is rejected. These are relatively
    # stringent requirements.
    pct_Assigned = assigned_reads / total_mapped
    pct_unmapped = N_unmapped / total_reads
    pct_noFeat = N_noFeat / total_mapped

    # Apply QC filters
    if (pct_Assigned >= .8) and (pct_unmapped <= .15) and (pct_noFeat <= .2):
        sample_ID = [file_id_to_barcode[fid]]
        ind_Genes = gene_rows['gene_id']
        raw_dfs.append(pd.DataFrame(columns=sample_ID, index=ind_Genes, data=gene_rows['unstranded'].values))
        tpm_dfs.append(pd.DataFrame(columns=sample_ID, index=ind_Genes, data=gene_rows['tpm_unstranded'].values))
        fpkm_dfs.append(pd.DataFrame(columns=sample_ID, index=ind_Genes, data=gene_rows['tpm_unstranded'].values))
    else:
        # remove bad sample from metadata
        trim_metadata = trim_metadata[trim_metadata['file_id'] != fid]

# Concatenate files that passed QC into gene X sample dataframes
rawDF = pd.concat(raw_dfs, axis=1)
tpmDF = pd.concat(tpm_dfs, axis=1)
fpkmDF = pd.concat(fpkm_dfs, axis=1)

# save data as raw counts, tpm normalized, fpkm normalized, and the meta data files
rawDF.to_csv('star_counts/raw_counts_test.csv')
tpmDF.to_csv('star_counts/tpm_counts_test.csv')
fpkmDF.to_csv('star_counts/fpkm_counts_test.csv')
trim_metadata.to_csv('star_counts/trim_metadata_test.csv', index=False)
