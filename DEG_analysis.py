# this script will perform differential expression analysis on the data from process_count_DL.py
# Note: rather than use the R implementation of DESeq2, I am using the python implementation to keep
# this project 'monolingual'. the paper on pyDESeq2 can be found here:https://doi.org/10.1093/bioinformatics/btad547

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pandas as pd
import numpy as np

# load sample data and metadata dataframes
samples = pd.read_csv('star_counts/raw_counts.csv', index_col=0)
full_metadata = pd.read_csv('star_counts/trim_metadata.csv')

# Reshape data: pyDESeq2 wants data to have sampleIDs x features (genes) shape, so transpose samples
# then trim metadata to a subset for DE with the sampleIDs as index
samples = samples.T
trimMD = pd.DataFrame({"sample_type": full_metadata['sample_type'].values}, index=full_metadata['sample_barcode'])

# now I am interested in comparing Primary Tumor and Solid Tissue Normal RNAseq outputs using DESeq (this is not the only
# comparison that can be conducted, using this framework any categorical data you can retrieve for your samples could be
# analyzed in this style of differential expression)
NormTumor = [True if (i == "Primary Tumor") or (i == "Solid Tissue Normal") else False for i in trimMD.sample_type.values]
samples = samples.loc[NormTumor, :]
trimMD = trimMD.loc[NormTumor, :]

# finally, as a quick point of QC prior to running the DESeq algorithm, I am removing genes where >50% of samples have
# no expression (count = 0 in raw counts). There are smarter ways to do this (like preserving presence/abscence genes)
# but those are a little more in depth than necessary for now
NumSamples = samples.shape[0]
samples = samples[[i for i in samples.columns if np.count_nonzero(samples[i].values) > NumSamples*.5]]
# for my samples (I will include the manifest I used on github to make this a little more reproducible in the future)
# I dropped about 30k genes in this way

# single facor analysis
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(counts=samples, metadata=trimMD, design_factors="sample_type", refit_cooks=True, inference=inference)
dds.deseq2()
print(dds)

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
ds.results_df.to_csv('star_counts/DESeq2_stats.csv')
NormDF = pd.DataFrame(index=samples.index, columns=samples.columns, data=dds.layers['normed_counts'])
NormDF.to_csv('star_counts/DESeq2_Norm.csv')
trimMD.to_csv('star_counts/DESeq2_MetaData.csv')