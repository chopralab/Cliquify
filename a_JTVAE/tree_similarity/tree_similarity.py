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


def tree_similarity2(idx):
    chosen_smiles = smiles_list[idx]
    tree = MolTree(chosen_smiles)

    chosen_graph = nx.Graph()
    chosen_graph.add_edges_from(tree.edges)
    # chosen_graph.add_nodes_from(tree.edges)

    ged_value = []
    for i, smiles in enumerate(smiles_list[:5000]):
        temp_tree = MolTree(smiles)
        temp_graph = nx.Graph()
        temp_graph.add_edges_from(temp_tree.edges)
        # temp_graph.add_nodes_from(temp_tree.edges)


        # ged_value.append(nx.graph_edit_distance(chosen_graph, temp_graph))
        val = nx.graph_edit_distance(chosen_graph, temp_graph, timeout=2)
        # val = val if val else 0
        if val: ged_value.append(val)
        if i % 200 == 0: print(i)
    
    # ged_mean = np.mean(np.array(ged_value))
    ged_mean = sum(ged_value) / len(ged_value)

    # with open("ged_mean_values2.txt", "a") as myfile:
    #     myfile.writelines("{},{}\n".format(idx, ged_mean))

    return idx, ged_mean

# with open("ged_mean_values.txt", "w") as myfile:
#     ged_mean_list = []
#     for i, smiles in enumerate(smiles_list[:50]):
#         ged = tree_similarity(i)
#         ged_mean = np.mean(np.array(ged))
#         print(ged_mean)
#         ged_mean_list.append(ged_mean)


ged_mean_list = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    val = [i for i in range(len(smiles_list[:50]))]
    ged_mean_list = executor.map(tree_similarity2, val)


ged_mean_list = list(ged_mean_list)
print(ged_mean_list)

with open("ged_mean_values2.txt", "w") as myfile:
    ged_mean_list_new = []
    for i, ged_mean in ged_mean_list:
        myfile.writelines("{},{}\n".format(i, ged_mean))
        ged_mean_list_new.append(ged_mean)


with open("ged_mean_values2.txt", "a") as myfile:
    myfile.writelines("--------------------\n")


plt.bar(x=[i for i, smiles in enumerate(smiles_list[:50])], height=ged_mean_list_new)
# plt.xticks(rotation = 90)
plt.xlabel('smiles idx')
plt.ylabel('ged score')
plt.savefig("tree_ged_score.png")

