# this script will use the DESeq2 outputs to generate putatitve gene interaction networks from the data

import pandas as pd
import numpy as np
import networkx as nx
import pickle

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
# This step could be done in any of a number of ways: modules from WGCNA, Louvain/leiden community detection,
# label propagation by finding highly central genes, k-cliques, etc.
# For computational simplicity, I will proceed using the basic process of identifying subgraphs in an adjacency matrix
# based on the correlation matrix. From there collections of genes can be evaluated for their relationship to the
# classification task (i.e. are all genes upregulated/downregulated, most significant/proportion of top 100
# most significant genes, etc.). From there the analysis can hone in on just a few graphs for the classification task.

# first generate a gene x gene correlation matrix, convert to absolute correlation
corr_Mat = np.abs(np.corrcoef(filtDF.values, rowvar=False))

# then convert to an adjacency matrix, using the .7 threshold to identify strongly correlated genes as being connected
corr_Mat[corr_Mat > .7] = 1
corr_Mat[corr_Mat < .7] = 0

# convert adjacency matrix to networkx graph
G = nx.convert_matrix.from_numpy_array(corr_Mat)

connected_components = nx.connected_components(G)
subgraphs = (G.subgraph(c) for c in connected_components)

# discard the subgraphs with fewer than 5 members (this will be the majority of the subgraphs identified.)
# I elected to keep the graphs with at least 5 connected members be
subgraphs = [i for i in subgraphs if len(i.nodes) > 4]
NodeCounts = [len(i) for i in subgraphs]
print(f"{len(subgraphs)} subgraphs remaining after filtering")
print(f"the largest graph is {max(NodeCounts)} genes and the median graph size is {np.median(NodeCounts)}")

# Now we can filter our subgraphs based on the direction of their expression, their significance in DESeq2 (Note:
# DESeq2 significance != biological significance)

graphDF = pd.DataFrame(columns=['size', 'Prop.Upreg', 'Prop.Top1000LFC', 'Prop.Top1000padj'], index=['g_'+str(i) for i in range(1,len(subgraphs)+1)])
absLFC = np.abs(statsDF['log2FoldChange'].values)
absLFC.sort()
LFCThresh = absLFC[-1000]

padj = statsDF['padj'].values
padj.sort()
padjThresh = padj[1000]

for gNum in range(len(subgraphs)):
    ind = 'g_' + str(gNum+1)
    graph = subgraphs[gNum]
    g_genes = [normDF.columns.tolist()[i] for i in graph.nodes()]
    graphDF.loc[ind, 'size'] = len(g_genes)
    graphDF.loc[ind, 'Prop.Upreg'] = sum(statsDF.loc[g_genes]["log2FoldChange"].values > 0)/len(g_genes)
    graphDF.loc[ind, 'Prop.Top1000LFC'] = (sum(np.abs(statsDF.loc[g_genes]["log2FoldChange"].values) > LFCThresh)/len(g_genes))
    graphDF.loc[ind, 'Prop.Top1000padj'] = (sum(statsDF.loc[g_genes]["padj"].values < padjThresh)/len(g_genes))

# explanations for the metrics I chose to track
# size: important for limiting the graphs we take into the next step, as larger graphs will be more computationally
# intensive
# Prop.UpReg: here I am looking for genes that did the same thing ( all upregulated or all downregulated) indicating
# regions of global up or down regulation
# Prop.Top1000LFC: how many of the genes where in the top 1000 based on LFC, my thinking here is that genes that undergo
# greater swings in expression will be easier for the models to get a good sense for
# Prop.Top1000padj: a bit of a tie-breaker column, not intended to be the ultimate deciding factor, but could inform
# the choice between close graphs.

# selecting graphs: I am planning to select one based on percent upregulated, one based on percent downregulated (0.0 in
# Prop.UpReg), and then the graph with the greatest percent in top 1000 LFC

graphDF = graphDF.sort_values(by='Prop.Upreg', ascending=False)
graphUR = graphDF.sort_values(by='Prop.Upreg', ascending=False).index[0]
graphDR = graphDF.sort_values(by='Prop.Upreg', ascending=True).index[0]
graphLFC = graphDF.sort_values(by='Prop.Top1000LFC', ascending=False).index[0]

# save the graphDF as a csv, then create a graph dictionary to save the subgraphs with the g_ label
graphDF.to_csv('graphs/graphDataFrame.csv', index=True)

graphDict = {}
for gID in [graphUR, graphDR, graphLFC]:
    graphDict[gID] = subgraphs[int(gID.split('_')[1])-1]

with open("graphs/graph_dictionary.pkl", "wb") as f:
    # Use pickle.dump() to write the dictionary to the file
    pickle.dump(graphDict, f)

# now the next steps will be to construct a GNN and construct graphs for each sample
