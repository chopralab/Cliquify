# import torch
# import torch.nn as nn

# x = torch.tensor([[1.0, -1.0],
#                   [0.0,  1.0],
#                   [0.0,  0.0]])

# in_features = x.shape[1]  # = 2
# out_features = 2

# m = nn.Linear(in_features, out_features)
# print(m.weight)
# print(m.bias)

# y = m(x)

# print(y)

import networkx as nx

G = nx.Graph()
G.add_edge(*[1, 2])
G.add_edge(*[4, 2])
G.add_edge(*[7, 2])

print(G.edges)

