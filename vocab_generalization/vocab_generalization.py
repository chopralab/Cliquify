import sys
sys.path.append('../')
from rdkit import Chem
import networkx as nx
import numpy as np
import networkx.algorithms.isomorphism as iso

from enumerate_random import MolTreeNode, dfs_random_assemble, remove_edges_reset_idx, set_atommap, set_atommap_graph
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol

with open("vocab2.txt") as f:
    smiles_list = f.readlines()

fragments = ["C.CC", "CCC"]

np.random.seed(42)

MAX_COUNT = 11

def random_sample_tree(root, depth):
    global count
    global edges_list
    global nodes_dict

    if count > MAX_COUNT:
        return
        
    # neighbors_count = np.random.randint(0, 5-depth) # max of 4 neigbors
    neighbors_count = np.random.randint(1,3) # max of 4 neigbors


    for i in range(neighbors_count):
        if count > MAX_COUNT: return

        count += i + 1
        topology = np.random.rand()

        if topology > 0.5:
            idx = np.random.randint(0, len(fragments))
            neigh = MolTreeNode(Chem.MolFromSmiles(fragments[idx]), smiles=fragments[idx])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh

            random_sample_tree(neigh, depth+1)

        else:
            neigh = MolTreeNode(Chem.MolFromSmiles(fragments[1]), smiles=fragments[1])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh



gen_smiles_list = []
for _ in range(1000):

    # if _ != 1:
    #     continue

    #--------------------------------------------------
    count = 0
    root = MolTreeNode(Chem.MolFromSmiles(fragments[1]), smiles=fragments[1])
    root.idx = 0
    nodes_dict = {root.idx: root}
    edges_list = []

    random_sample_tree(root, 0)

    for x,y in edges_list:
        nodes_dict[x].add_neighbor(nodes_dict[y])
        nodes_dict[y].add_neighbor(nodes_dict[x])

    # print(edges_list)    
    
    list_of_nodes = list(nodes_dict.values())

    for i,node in enumerate(list_of_nodes):
        node.nid = i + 1
        if len(node.neighbors) > 1: #Leaf node mol is not marked
            set_atommap(node.mol, node.nid)
            set_atommap(node.tri_mol, node.nid)
            set_atommap_graph(node.graph, node.nid)
        node.is_leaf = (len(node.neighbors) == 1)

    #---------------------------------------------------

    root_idx = 0
    root = list_of_nodes[root_idx]
    cur_graph = root.graph.copy()
    global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
    global_amap[root_idx+1] = {node:node for node in cur_graph.nodes()}

    dfs_random_assemble(cur_graph, global_amap, [], root, None, print_out=False)

    # draw_mol(cur_graph, 9999, folder="../vocab_generalization/subgraph")
    final_graph = remove_edges_reset_idx(cur_graph)
    # draw_mol(final_graph, 10_000, folder="../vocab_generalization/subgraph")

    # remove unconnected nodes remnants 
    final_graph.remove_nodes_from(list(nx.isolates(final_graph)))

    cur_mol = nx_to_mol(mol_to_nx(nx_to_mol(final_graph)))
    dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)
    dec_mol  = Chem.MolFromSmiles(dec_smiles)

    if cur_mol and dec_mol:
        print(dec_smiles)

        # # if dec_smiles == "CCCC1CCCC2CCC2C1":
        # if dec_smiles == "C1CC2(C1)CCC2":
        #     draw_mol(cur_graph, 7778, folder="../vocab_generalization/subgraph")
        #     print(_)
        #     raise

        if dec_smiles not in gen_smiles_list:
            gen_smiles_list.append(dec_smiles)

    else:
        print(None)


with open("gen_mol_count_test_nei_enclosed(no_bridge).txt", "a") as myfile:
    for gen_smiles in gen_smiles_list:
        myfile.writelines("{}\n".format(gen_smiles))
    # raise

