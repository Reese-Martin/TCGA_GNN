# this script will use the DESeq2 outputs to generate putatitve gene interaction networks from the data

import pandas as pd
import numpy as np
import networkx as nx

### Load data, DESeq stats, and sample metadata
normDF = pd.read_csv('star_counts/DESeq2_Norm.csv', index_col=0)
# update gene names removing the .XX transcript version to make searching more straightforward.
normDFcols = [i.split('.')[0] for i in normDF.columns]
normDF.columns = normDFcols

statsDF = pd.read_csv('star_counts/DESeq2_stats.csv', index_col=0)
statsDFind = [i.split('.')[0] for i in statsDF.index]
statsDF.index = statsDFind

MD = pd.read_csv('star_counts/DESeq2_metadata.csv', index_col=0)

### filter full dataset to consider just the sig DE genes (using adjusted significance)
filtDF = normDF[[i for i in statsDF.index if statsDF.loc[i, 'padj'] < .05]]

### Identify subgraphs of importance in the network
# This step could be done in any of a number of ways: modules from WGCNA, Louvain community detection, label propagation
# by finding highly central genes, etc. For simplicity (and speed) I will proceed with label propagation.
# if this ends up failing downstream I will likely revisit this step

# first generate a gene x gene correlation matrix and convert to absolute correlation (networkX label propagation does
# not consider directed edges)
corr_Mat = np.abs(np.corrcoef(filtDF.values, rowvar=False))
# then convert to an adjacent, using the .7 threshold to identify strongly correlated genes as being connected
corr_Mat[corr_Mat > .7] = 1
corr_Mat[corr_Mat < .7] = 0

# convert adjacency matrix to networkx graph
G = nx.convert_matrix.from_numpy_array(corr_Mat)

