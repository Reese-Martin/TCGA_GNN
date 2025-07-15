# This script loads the networks generated in Network_Analysis.py, generates a per sample
# graph where the edges are shared across all samples, but the nodes (genes) carry that samples specific expression
# level. These graphs are then fed into a GNN constructed using torch_geometric
from torch_geometric.utils import from_networkx
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
import torch
import pickle as pkl
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.utils.class_weight import compute_class_weight
import random
# from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay


class GCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, num_classes):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.classifier = torch.nn.Linear(hidden_channels, num_classes)

    def forward(self, x, edge_index, batch):
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()
        x = global_mean_pool(x, batch)  # Aggregate node features per sample
        return self.classifier(x)

    def copy(self):
        new_model = GCN(
            in_channels=self.conv1.in_channels,
            hidden_channels=self.conv1.out_channels,
            num_classes=self.classifier.out_features
        )
        new_model.load_state_dict(self.state_dict())
        return new_model


def prepare_graphs(graph, samples, geneDF, metadataDF):  # this assumes you want a single graph template for all samples
    # convert networkx to torch_geometric
    graph_template = from_networkx(graph)
    # Save edge_index (used for all samples)
    edge_template = graph_template.edge_index

    labels = metadataDF['sample_type'].values
    le = LabelEncoder()
    metadataDF['int_label'] = le.fit_transform(labels)
    classCount = len(le.classes_)

    sample_graphs = []
    for sample in samples:
        # [num_nodes, 1], in the future node features could be added as additional columns
        node_feats = torch.tensor(geneDF.loc[sample].values, dtype=torch.float).unsqueeze(1)

        graph = Data(x=node_feats, edge_index=edge_template,
                     y=torch.from_numpy(np.array(metadataDF.loc[sample, 'int_label'])).long())
        sample_graphs.append(graph)

    return sample_graphs, classCount


def model_training(model, loader, optimizer, loss_fn, epochs):
    loss_track = np.zeros((epochs, len(loader)))
    batchAcc = np.zeros((epochs, len(loader)))
    for epoch in range(epochs):
        BtchCnt = 0
        for batch in loader:
            out = model(batch.x, batch.edge_index, batch.batch)
            loss = loss_fn(out, batch.y)
            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
            # measure accuracy
            pred = out.argmax(dim=1)
            correct = (pred == batch.y).sum().item()
            accuracy = correct / batch.y.size(0)
            batchAcc[epoch, BtchCnt] = accuracy
            loss_track[epoch, BtchCnt] = loss.item()
            BtchCnt = BtchCnt + 1
        print(f'Epoch {epoch + 1}/{epochs}, Loss: {np.mean(loss_track[epoch, :])}')
        print(f'Accuracy: {np.mean(batchAcc[epoch,:]):.2f}')

    return loss_track, batchAcc


# load networks and expression data
with open('graphs/graph_dictionary.pkl', 'rb') as file:
    graph_dict = pkl.load(file)

graphDF = pd.read_csv('graphs/graphDataFrame.csv', index_col=0)
metadataDF = pd.read_csv('star_counts/DESeq2_MetaData.csv', index_col=0)
normExpDF = pd.read_csv('star_counts/DESeq2_Norm.csv', index_col=0)


# training for upregulated network
upregGraph = graph_dict[list(graph_dict.keys())[0]]
upregGraphGenes = normExpDF.iloc[:, list(upregGraph.nodes)]
samples = normExpDF.index.values

# for each sample build a graph with gene expression data as node data
upreg_sample_graphs, classNum = prepare_graphs(upregGraph, samples, upregGraphGenes, metadataDF)

# Get the unique classes
classes = np.unique(metadataDF['int_label'])
weights = compute_class_weight(class_weight='balanced', classes=classes, y=metadataDF['int_label'])

# Convert to tensor
class_weights = torch.tensor(weights, dtype=torch.float)

# Prepare loader, (the loader groups graphs together to speed up the training process)
upreg_loader = DataLoader(upreg_sample_graphs, batch_size=16, shuffle=True)

# Model definition and training parameters
upreg_model = GCN(in_channels=1, hidden_channels=64, num_classes=classNum)
optimizer = torch.optim.Adam(upreg_model.parameters(), lr=1e-3)
loss_fn = torch.nn.CrossEntropyLoss()
epochs = 1000

upRegLoss, upRegAcc = model_training(upreg_model, upreg_loader, optimizer, loss_fn, epochs)

# note: model loss is low, and accuracy is high, but this could be an artifact of having the dataset be severely biased
# The normal tissue is underrepresented by ~9:1 so the model seems to be always predicting that the graph is
# tumor and accuracy is high because of that (note that accuracy sits at .9 for several hundred generations)

# let's see what happens when we use a weighted loss function
# Get the unique classes
classes = np.unique(metadataDF['int_label'])
weights = compute_class_weight(class_weight='balanced', classes=classes, y=metadataDF['int_label'])
# Convert to tensor
class_weights = torch.tensor(weights, dtype=torch.float)

wloss_fn = torch.nn.CrossEntropyLoss(weight=class_weights)
upreg_w_model = GCN(in_channels=1, hidden_channels=64, num_classes=classNum)
epochs = 1000
upRegLoss_w, upRegAcc_w = model_training(upreg_w_model, upreg_loader, optimizer, wloss_fn, epochs)

# It seems that even class weighting isn't enough to push the models towards correctly classifying the minority
# class, although it did serve the purpose of making the loss spectacularly high. I will evaluate the other 2 networks
# saved, but may need to go back to the graph construction drawing board or find more normal samples
# (I would like to avoid other class balancing efforts like SMOTE for now)

# training for downregulated networks
dregGraph = graph_dict[list(graph_dict.keys())[1]]
dregGraphGenes = normExpDF.iloc[:, list(dregGraph.nodes)]

# for each sample build a graph with gene expression data as node data
dreg_sample_graphs, classNum = prepare_graphs(dregGraph, samples, dregGraphGenes, metadataDF)
dreg_loader = DataLoader(dreg_sample_graphs, batch_size=16, shuffle=True)

dreg_w_model = GCN(in_channels=1, hidden_channels=64, num_classes=classNum)
epochs = 1000
dRegLoss_w, dRegAcc_w = model_training(dreg_w_model, dreg_loader, optimizer, wloss_fn, epochs)

# training for LFC networks
LFCregGraph = graph_dict[list(graph_dict.keys())[2]]
LFCregGraphGenes = normExpDF.iloc[:, list(dregGraph.nodes)]

# for each sample build a graph with gene expression data as node data
LFCreg_sample_graphs, classNum = prepare_graphs(LFCregGraph, samples, LFCregGraphGenes, metadataDF)
LFCreg_loader = DataLoader(LFCreg_sample_graphs, batch_size=16, shuffle=True)

LFCreg_w_model = GCN(in_channels=1, hidden_channels=64, num_classes=classNum)
epochs = 1000
LFCRegLoss_w, LFCRegAcc_w = model_training(LFCreg_w_model, LFCreg_loader, optimizer, wloss_fn, epochs)

# seems that there were no magic bullets in the other graphs (very similar accuracies and log losses.)
