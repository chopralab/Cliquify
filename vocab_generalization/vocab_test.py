import sys
sys.path.append('../')
from rdkit import Chem
import networkx as nx
import numpy as np
import networkx.algorithms.isomorphism as iso
import random
import pickle
import json
import matplotlib.pyplot as plt
import pickle


from enumerate_random import MolTreeNode, dfs_random_assemble, remove_edges_reset_idx, set_atommap, set_atommap_graph
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol, data_to_mol
from search_tree_vocab import get_common_vocabs, tree_vocab_byAtom, tree_vocab_byBond, tree_vocab_byBondType, vocab_usage


KEY = "ghost"
folder="../vocab_generalization/tree"
idx = 1

with open("{}/show_{}_tree.pkl".format(folder, idx), 'rb') as handle:
    nodes_dict = pickle.load(handle)

    list_of_nodes = list(nodes_dict.values())
    for i,node in enumerate(list_of_nodes):
        node.nid = i + 1
        if len(node.neighbors) > 1: #Leaf node mol is not marked
            set_atommap(node.tri_mol, node.nid)
            set_atommap_graph(node.graph, node.nid)
        node.is_leaf = (len(node.neighbors) == 1)

    # fix boolprop not saved in pickle
    for node in nodes_dict.values():
        node.tri_mol = nx_to_mol(node.graph)

    root_idx = 0
    root = list_of_nodes[root_idx]
    cur_graph = root.graph.copy()
    global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
    global_amap[root_idx+1] = {node:node for node in cur_graph.nodes()}


    dfs_random_assemble(cur_graph, global_amap, [], root, None, print_out=False)

    draw_mol(cur_graph, idx, folder="../vocab_generalization/subgraph")