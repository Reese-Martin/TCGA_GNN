# this script will perform differential expression analysis on the data from process_count_DL.py
# Note: rather than use the R implementation of DESeq2, I am using the python implementation to keep
# this project 'monolingual'. the paper on pyDESeq2 can be found here:https://doi.org/10.1093/bioinformatics/btad547

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pandas as pd
import numpy as np


def DEG_analysis(loc, rawFile, metadata, factor, group1, group2, skip_zr=True,analysis='normalize'):
    # load sample data and metadata dataframes
    samples = pd.read_csv(loc + rawFile, index_col=0)
    full_metadata = pd.read_csv(loc + metadata)

    # Reshape data: pyDESeq2 wants data to have sampleIDs x features (genes) shape, so transpose samples
    # then trim metadata to a subset for DE with the sampleIDs as index
    samples = samples.T
    trimMD = pd.DataFrame({factor: full_metadata[factor].values}, index=full_metadata['sample_barcode'])

    # drop samples with NaNs in the factor column
    samples = samples[trimMD[factor].notna()]
    trimMD = trimMD[trimMD[factor].notna()]

    # trim to just samples comparing the groups of interest (i.e. if there are 3 groups in the specified factor,
    # focus the analysis on two of them)
    y = [True if (i == group1) or (i == group2) else False for i in trimMD[factor]]
    samples = samples.loc[y, :]
    trimMD = trimMD.loc[y, :]

    # as a quick point of QC prior to running the DESeq algorithm, optionally remove genes where >50% of samples have
    # no expression (count = 0 in raw counts).
    if not skip_zr:
        NumSamples = samples.shape[0]
        samples = samples[[i for i in samples.columns if np.count_nonzero(samples[i].values) > NumSamples*.5]]
    # for my samples (I will include the manifest I used on github to make this a little more reproducible in the future)
    # I dropped about 30k genes in this way

    # single factor analysis
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(counts=samples, metadata=trimMD, design_factors=factor, refit_cooks=True, inference=inference)
    dds.deseq2()
    print(dds)
    if analysis == 'normalize':
        NormDF = pd.DataFrame(index=samples.index, columns=samples.columns, data=dds.layers['normed_counts'])
        print("NormDF size: ", NormDF.shape)
        NormDF.to_csv(loc + 'DESeq2_Norm.csv')
        trimMD.to_csv(loc + 'DESeq2_MetaData.csv')

    elif analysis == "differential expression":
        # now we can see things like dispersion and LogFoldChange from the .varm section of the AnnData object
        dds.varm['LFC'].head()
        dds.varm['dispersions'].head()

        #Statistical Analysis
        # Now that dispersions and LFCs have been calculated, stats can be calculated
        ds = DeseqStats(dds, contrast=["sample-type", "Solid Tissue Normal", "Primary Tumor"], inference=inference)
        ds.summary()
        ds.results_df.head()
        # Note: the L2FC values here follow X vs Y specified in dds.varm['LFC'], so negative l2FC means that
        # the gene is down regulated in the X case compared to the Y (here that would mean Solid Tissue Normal
        # has lower expression than Primary Tumor)

        # Now we have generated a DESeq dataframe with normalized counts, and a stats results dataframe with
        # significance values. Save both for future use, and in next steps we will can filter data using DES
        # significance values.
        ds.results_df.to_csv(loc + 'DESeq2_stats.csv')
        NormDF = pd.DataFrame(index=samples.index, columns=samples.columns, data=dds.layers['normed_counts'])
        NormDF.to_csv(loc + 'DESeq2_Norm.csv')
        trimMD.to_csv(loc + 'DESeq2_MetaData.csv')


direct = 'STN_star_counts/'
rFile = 'raw_counts.csv'
mdFile = 'trim_metadata.csv'
# note, if you are strictly interested in the normalization, factor and groups do not matter
factor = 'gender'
group1 = 'female'
group2 = 'male'

DEG_analysis(direct, rFile, mdFile, factor, group1, group2, skip_zr=True, analysis='normalize')

