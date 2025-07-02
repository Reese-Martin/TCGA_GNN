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

# load networks and expression data
with open('graphs/graph_dictionary.pkl', 'rb') as file:
    graph_dict = pkl.load(file)

graphDF = pd.read_csv('graphs/graphDataFrame.csv', index_col=0)
metadataDF = pd.read_csv('star_counts/DESeq2_MetaData.csv', index_col=0)
normExpDF = pd.read_csv('star_counts/DESeq2_Norm.csv', index_col=0)


upregGraph = graph_dict[list(graph_dict.keys())[0]]
upregGraphGenes = normExpDF.iloc[:, list(upregGraph.nodes)]

# for now this step will be done with just the upreg graphs in mind, but when I get the classification task working
# for those networks I will functionalize this section to make it applicable to the other two graph types we have saved
# so far

# convert networkx to torch_geometric
graph_template = from_networkx(upregGraph)
# Save edge_index (used for all samples)
edge_template = graph_template.edge_index

labels = metadataDF['sample_type'].values
le = LabelEncoder()
metadataDF['int_label'] = le.fit_transform(labels)
classCount = len(le.classes_)

samples = normExpDF.index.values
sample_graphs = []
for sample in samples:
    # [num_nodes, 1], in the future node features could be added as additional columns
    node_feats = torch.tensor(upregGraphGenes.loc[sample].values, dtype=torch.float).unsqueeze(1)

    graph = Data(x=node_feats, edge_index=edge_template,
                 y=torch.from_numpy(np.array(metadataDF.loc[sample, 'int_label'])).long())
    sample_graphs.append(graph)


# Prepare loader, the loader groups graphs together to speed up the training process
loader = DataLoader(sample_graphs, batch_size=16, shuffle=True)

# Training loop (basic)
model = GCN(in_channels=1, hidden_channels=64, num_classes=classCount)
test = model.copy()  # tmp bug fixing copied the model to make sure the training was actually taking place
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
loss_fn = torch.nn.CrossEntropyLoss()

epochs = 1000
loss_track = np.zeros((epochs, len(loader)))
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
        print(f'Accuracy: {accuracy:.2f}')
        loss_track[epoch, BtchCnt] = loss.item()
        BtchCnt = BtchCnt + 1
    print(f'Epoch {epoch+1}/{epochs}, Loss: {np.mean(loss_track[epoch, :])}')

# note: model loss is low, and accuracy is high, but this is an artifact of having the dataset be severely biased
# The normal tissue is underrepresented by ~9:1 so the model is just always predicting that the graph is
# tumor and accuracy is high because of that. We will fix that.
