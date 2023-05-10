import sys
sys.path.append('../')
from rdkit import Chem
import networkx as nx
import numpy as np
import networkx.algorithms.isomorphism as iso
import itertools
import os

from enumerate_all import MolTreeNode, dfs_all_assemble, remove_edges_reset_idx, set_atommap, set_atommap_graph
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol, data_to_mol



def get_trees(nodes):
    # generate all pairwise combinations of nodes
    edges =  [a for a in itertools.product(range(nodes), range(nodes))]

    # use sets to lose..
    # ..symmetric edges: (0,1), (1,0) => keep only (0,1) 
    edges = list(set([tuple(set(e)) for e in edges]))
    # ..and self-loops: (0,0)
    edges = [e for e in edges if len(e)>1]

    trees = []
    # generate all graphs that have nodes-1 edges
    for o in itertools.combinations(edges, nodes-1):
        #make sure that all nodes are in the edgelist:
        flattened = [item for sublist in o for item in sublist]
        
        if len(set(flattened)) == nodes:
            G = nx.Graph()
            G.add_edges_from(o)
            
            nei_larger_2 = [1 if len(nei) > 3 else 0 for node, nei in G.adjacency()] 
            if sum(nei_larger_2): continue
                
            # make sure all nodes are connected
            if len(list(nx.connected_components(G)))==1:
                trees.append(G)

    return trees

def get_combinations_vocab(nums, r):
    combos = list(itertools.combinations(nums, r))
    return combos

def get_combinations_attachment(nums):
    combos = list(itertools.product(*(range(num) for num in nums)))
    return combos

def prepare_mol_tree(tree, vocab_list, attachment_idxs=None):
    nodes_dict = {}
    for idx, vocab_idx in enumerate(vocab_list): 
        mol = data_to_mol("../vocab_generalization/graph_vocab/{}".format(vocab_idx))
        node = MolTreeNode(mol)
        node.idx = idx
        nodes_dict[idx] = node
    
    for x, y in tree.edges:
        nodes_dict[x].add_neighbor(nodes_dict[y])
        nodes_dict[y].add_neighbor(nodes_dict[x])
    
    list_of_nodes = list(nodes_dict.values())

    for i, node in enumerate(list_of_nodes):
        node.nid = i + 1
        if len(node.neighbors) > 1: #Leaf node mol is not marked
            set_atommap(node.mol, node.nid)
            set_atommap(node.tri_mol, node.nid)
            set_atommap_graph(node.graph, node.nid)
        if attachment_idxs:
            node.cand_idx = attachment_idxs[i]
        node.is_leaf = (len(node.neighbors) == 1)

    #-----------------------------------------------
    root_idx = 0
    root = list_of_nodes[root_idx]
    cur_graph = root.graph.copy()
    global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
    global_amap[root_idx+1] = {node:node for node in cur_graph.nodes()}

    return cur_graph, global_amap, root

def main():

    n = 5
    trees = get_trees(n)

    VOCAB_SIZE = len(os.listdir("../vocab_generalization/graph_vocab"))

    vocab_combination = get_combinations_vocab(range(VOCAB_SIZE), n)

    input = itertools.product(trees, vocab_combination)


    for tree, vocab_list in input:
        
        cur_graph, global_amap, root = prepare_mol_tree(tree, vocab_list)

        # -----------------------FIND ALL POSSIBLE ATTACHMENT-------------------------#

        possible_attachment = {}
        dfs_all_assemble(cur_graph, global_amap, [], root, None, possible_attachment)

        print(possible_attachment)

        attachment = [possible_attachment[i+1] for i in range(len(possible_attachment))]
        attachment_idxs = get_combinations_attachment(attachment)        


        #------------------ENUMERATE ALL POSSIBLE ATTACHMENTS-------------------------#

        for attachment_idx in attachment_idxs:
            prepare_mol_tree(tree, vocab_list, attachment_idx)

        final_graph = remove_edges_reset_idx(cur_graph)

        # if it forms valid graph
        if len(list(nx.connected_components(final_graph))) == 1:
            cur_mol = nx_to_mol(mol_to_nx(nx_to_mol(final_graph)))
            dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)
            dec_mol  = Chem.MolFromSmiles(dec_smiles)

            if cur_mol and dec_mol:
                return dec_smiles

    
    nums = [5, 2, 4, 2]
    combos = get_combinations_attachment(nums)
    print(combos)

main()