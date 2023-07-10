import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
from multiprocessing import Pool
import pickle
import os

from enumerate3 import dfs_assemble, node_labelling, remove_edges_reset_idx, reconstruction_evaluation, MolTreeNode
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, draw_mol, data_to_mol, node_exact, edge_exact
import matplotlib.pyplot as plt

import timeit
with open("../zdata/zinc/all.txt") as f:
    smiles_list = f.readlines()


def vocab(idx):
    return data_to_mol("../vocab_generalization/graph_vocab/{}".format(idx))

def load_graph_list(file_path="/storage/fong/Cliquify/vocab_generalization/graph_vocab"):
    total_vocab_graphs=[]
    for filename in os.listdir(file_path):
        vocab_file_path = file_path + "/{}".format(filename)
        G = pickle.load(open(vocab_file_path, "rb"))
        total_vocab_graphs.append(G)
    return total_vocab_graphs

def match_graph(graph, total_vocab_graph):
    for g_idx, G in enumerate(total_vocab_graph):
        if nx.is_isomorphic(G, graph, node_match=node_exact, edge_match=edge_exact):
            return g_idx

def recover_graph(triangulated_graph, mol_node_u, mol_node_v):
    u_clique = mol_node_u.clique
    v_clique = mol_node_v.clique

    atoms_subgraph = list(set(u_clique + v_clique))
    assemble_label = triangulated_graph.subgraph(atoms_subgraph)

    return assemble_label

def save_tree(triangulated_graph, nodes_dict, edges_list, idx=7000, folder="../vocab_generalization/trees2"):
    plt.clf()
    graph = nx.Graph()
    saved_graph = nx.Graph()
    for u, v in edges_list:
        u_gid, v_gid = nodes_dict[u].gid, nodes_dict[v].gid
        graph.add_node(u, smiles=Chem.MolToSmiles(vocab(u_gid)))
        graph.add_node(v, smiles=Chem.MolToSmiles(vocab(v_gid)))
        graph.add_edge(u, v)

        saved_graph.add_node(u, mol_node=nodes_dict[u])
        saved_graph.add_node(v, mol_node=nodes_dict[v])
        assemble_node = recover_graph(triangulated_graph, nodes_dict[u], nodes_dict[v])
        saved_graph.add_edge(u, v, assemble=assemble_node)

    nx.write_gpickle(saved_graph, "{}/show{}_tree.pkl".format(folder, idx))

    # with open("{}/show{}_tree.pkl".format(folder, idx), 'wb') as handle:
    #     pickle.dump(nodes_dict, handle)

    pos = nx.spring_layout(graph,k=0.15,iterations=100)
    nx.draw(graph, pos)
    node_labels = nx.get_node_attributes(graph, "smiles")
    node_labels = {k : "     ({})".format(v) for k, v in node_labels.items()}
    nx.draw_networkx_labels(graph, pos, node_labels)
    # nx.draw(graph, pos, node_color="yellow", with_labels=True)
    plt.savefig("{}/show{}_tree.png".format(folder, idx))
    plt.clf()

unusable_ = []
usable = []
def prepare_tree(idx):
    chosen_smiles = smiles_list[idx]

    total_vocab_graphs = load_graph_list()
    mol = Chem.MolFromSmiles(chosen_smiles)
    
    cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)

    nodes_dict = {}
    list_of_nodes = []

    for i, clique in enumerate(cliques):
        if isinstance(clique, tuple):
            c = list(clique)
            m = MolTreeNode(mol, c)
            if not match_graph(m.graph, total_vocab_graphs):
                print(idx)
                print(i, clique)
                raise
            m.gid = match_graph(m.graph, total_vocab_graphs)
            nodes_dict[i] = m
        list_of_nodes.append(m)

    for x,y in molTreeEdges:
        list_of_nodes[x].add_neighbor(list_of_nodes[y])
        list_of_nodes[y].add_neighbor(list_of_nodes[x])

    
    save_tree(triangulated_graph, nodes_dict, molTreeEdges, idx=idx, folder="../space_time_benchmarking/trees")
    
    # try:
    #     usable.append(idx)
    #     print(usable)
    #     print(unusable_)
    #     print()
    # except:
    #     unusable_.append(idx)
    #     print(usable)
    #     print(unusable_)
    #     print()
    # # raise


    # root_idx = 0
    # root = list_of_nodes[root_idx]
    # cur_graph = root.graph.copy()
    # global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
    # global_amap[root_idx+1] = {node:node for node in cur_graph.nodes()}
    
    # dfs_assemble(cur_graph, global_amap, [], root, None, print_out=False)

    
    return 



for idx, smiles in enumerate(smiles_list):
    prepare_tree(idx)

print(unusable_)