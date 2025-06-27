# Overview
This repo contains the work I am doing to analyze TCGA data through the lens of network biology. I am trying to write this code so that it canbe flexibly applied to other TCGA data sources, and starting at the DESeq2 level the scripts
should be applicable to any RNAseq data.

Note: This is what I would call a learning repo, it is a personal improvement project so that I can learn new skills (particularly interested in learning to wrangle TCGA data, as well as the application of Graph Neural Networks for disease 
classification and biomarker discovery). In the future I may beautify the code here and share it/any results I find interesting, but for now it is a side project that I have a lot of interest in, rather than my day to day life.

# Goals
## 1. Access TCGA data in a programmatic way: 
   Done with Generate_Manifest.py and the gdc-client (using the client because downloading hundreds of files will get you rate limited/ timed out if you jsut try to use API calls
   a. Note: by changing values in the `Params` variable you can change which data you access as well as the project it orignites from

## 2. Retrieve sample metadata for downloaded samples: 
   Done with PT_metadata.py.

## 3. Trim metadata for successful downloads and generate unified count .csv (for raw, tpm, fpkm) from downloaded patient data: 
   done using process_count_DL.py

## 4. Differential Expression analysis: 
   Done using DEG_analysis.py. Becuase this is a learning opportunity, I elected to use the pyDESeq2 package rather than R and DESeq2. This was driven by curiosity to see if python had caught up to R
   in this realm and the results seem promising. DE genes make sense with the conditions analyzed, so that is promising. I would say this is effectively a drop in python native replacement for DESeq2 (which python has sorely needed)

## 5. Network Analysis: 
   In progress with Network_Analysis.py. So far, using gene correlation analysis on the significantly differentially expressed genes I have generated a GxG adjacency matrix. The next step is to identify communities
   of genes, but this is unfortunately a very computationally intensive task (..and my laptop is not up to it).
   
   Plan to overcome the harsh realities of computation: From the large ~22kx22k adjacency matrix, identify subgraphs. This matrix is not fully connected, so communities of genes will by necessity be disjoint (cannot span subgraphs)
   From here, there are two options: 1) community detection per subgraph and serialize the detection process for all graphs. 2) identify subgraphs that are associated with the classification task (i.e. all genes show same direction of
   DE, all genes are top 100 by abs(LFC), etc) and then community detection on the top subgraphs. My current plan is to proceed with 2 using the task as a heuristic to limit the graph space, but if nothing pops out then I will be 
   stuck with 1.

## 6. Graph classification: 
   Not started.
   This section of the project with combine the graphs from 5 (nodes are genes, edges are interactions between genes) with expression values for each node coming from a sample's noramlized RNAseq expression count and feed these graphs
   into a pytorch geometric GCN. I have worked with GCNs before in Julia using the Flux and GeometricFlux packages, so I am looking forward to this chance to expand my ML experience to the geometric_torch implementation.
   a. I am also itnerested in providing aditional node information to the GCN (like number of string interactions, previous association with disease state, etc) to see how that effects performance.

## 8. Comparison to traditional methods: 
   Not started.
   Here I plan to compare the GCN performance to other ML tools (regression on top 100 expressed genes, NN on same, RandomForest on same) to see if adding graph information improves results
   I will also compare with things like PCA plots and more traditional expression based analysis.
   
