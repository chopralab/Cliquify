import sys
sys.path.append('../')
sys.path.append('../jtnn')
from rdkit import Chem
import concurrent.futures
import networkx as nx
import multiprocessing
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


from jtnn import chemutils
from jtnn.mol_tree import MolTree

with open("../../zinc/all.txt") as f:
    smiles_list = f.readlines()


def tree_similarity(idx):
    chosen_smiles = smiles_list[idx]
    tree = MolTree(chosen_smiles)

    chosen_graph = nx.Graph()
    # chosen_graph.add_edges_from(tree.edges)
    chosen_graph.add_nodes_from(tree.edges)

    ged_value = []
    for i, smiles in enumerate(smiles_list[:1000]):
        temp_tree = MolTree(smiles)
        temp_graph = nx.Graph()
        # temp_graph.add_edges_from(temp_tree.edges)
        temp_graph.add_nodes_from(temp_tree.edges)


        ged_value.append(nx.graph_edit_distance(chosen_graph, temp_graph))


    return ged_value

with open("ged_mean_values.txt", "a") as myfile:
    ged_mean_list = []
    for i, smiles in enumerate(smiles_list[:50]):
        ged = tree_similarity(i)
        ged_mean = np.mean(np.array(ged))
        print(ged_mean)
        ged_mean_list.append(ged_mean)
        myfile.writelines("{},{}\n".format(i, ged_mean))


plt.bar(x=[i for i, smiles in enumerate(smiles_list[:50])], height=ged_mean_list)
# plt.xticks(rotation = 90)
plt.xlabel('smiles idx')
plt.ylabel('ged score')
plt.savefig("tree_ged_score.png")

# with concurrent.futures.ThreadPoolExecutor() as executor:
#     executor.map(tree_similarity, [i for i in range(len(smiles_list))])