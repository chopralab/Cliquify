import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
import multiprocessing
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import pickle
import json


from enumerate3 import get_fragments, get_fragments2, get_mol2, get_smarts_fragments, get_mol2, get_triangulated_graph
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, node_equal_iso, ring_edge_equal_iso, mol_to_data



with open("../zdata/zinc/all.txt") as f:
    smiles_list = f.readlines()


total_vocab_smiles = set()
total_vocab_smiles_smarts = set()
total_vocab_graph = []
graph_count = defaultdict(int)
cond_probability_one = defaultdict(int)
cond_probability_two = defaultdict(int)

total_one = defaultdict(int)
total_two = defaultdict(int)

def vocab_count(idx):
    chosen_smiles = smiles_list[idx]
    mol = Chem.MolFromSmiles(chosen_smiles)

    # remove stereochemistry
    gold_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
    mol = Chem.MolFromSmiles(gold_smiles)

    cliques, edges, _ = tree_decomp(mol)

    tree = nx.Graph()
    tree.add_edges_from(edges)
    cliqId_to_graphId = {}
    
    for i, clique in enumerate(cliques):
        if isinstance(clique, tuple):
            fragment_smiles = get_fragments(mol, clique)
            canon_smiles = Chem.CanonSmiles(fragment_smiles)
            fragment_smarts = get_smarts_fragments(mol, clique)
            _, fragment_graph = get_triangulated_graph(get_mol2(mol, clique))

            duplicate=False
            for g_idx, G in enumerate(total_vocab_graph):
                if nx.is_isomorphic(G, fragment_graph, node_match=node_equal_iso, edge_match=ring_edge_equal_iso):
                    duplicate = True
                    graph_count[g_idx] += 1
                    cliqId_to_graphId[i] = g_idx
                    break
            if duplicate: continue

            total_vocab_smiles.add(fragment_smiles)
            total_vocab_smiles_smarts.add(fragment_smiles + "|" + canon_smiles + "|" + fragment_smarts)
            total_vocab_graph.append(fragment_graph)
            
            graph_count[len(total_vocab_graph)] += 1
            cliqId_to_graphId[i] = len(total_vocab_graph)
    
    # for conditional sampling
    for st, ed in edges:
        u, v = cliqId_to_graphId[st], cliqId_to_graphId[ed]

        key = "{}|{}".format(u, v)
        cond_probability_one[key] += 1
        total_one[str(v)] += 1

        key = "{}|{}".format(v, u)
        cond_probability_one[key] += 1
        total_one[str(u)] += 1

        for nei in tree.neighbors(st):
            if nei == ed: continue
            nei = cliqId_to_graphId[nei]
            key = "{}|{},{}".format(u, v, nei)
            cond_probability_two[key] += 1
            total_two["{},{}".format(v, nei)] += 1

        for nei in tree.neighbors(ed):
            if nei == st: continue
            nei = cliqId_to_graphId[nei]
            key = "{}|{},{}".format(v, u, nei)
            cond_probability_two[key] += 1
            total_two["{},{}".format(u, nei)] += 1
    
    # printing
    if idx % 1000 == 0:
        print("smiles", len(total_vocab_smiles))
        print("smiles_smarts", len(total_vocab_smiles_smarts))
        print("graph", len(total_vocab_graph))
        print()
        # x = dict(sorted(graph_count.items(), key=lambda x:x[1], reverse=True))
        # pprint.pprint(graph_count)
        # print(x)


if __name__ == '__main__':

    for i in range(len(smiles_list)):
        vocab_count(i)

    print()
    total_vocab_smiles = list(total_vocab_smiles)
    print("total_vocab_smiles", len(total_vocab_smiles))
    print("total_vocab_graph", len(total_vocab_graph))

    # save graph vocabs
    for i, vocab_G in enumerate(total_vocab_graph):
        pickle.dump(vocab_G, open("../vocab_generalization/graph_vocab/{}".format(i), "wb"))

    with open('../vocab_generalization/cond_proba/cond_probability_one.json', 'w', encoding ='utf8') as json_file:
        json.dump(cond_probability_one, json_file) 
    with open('../vocab_generalization/cond_proba/cond_probability_two.json', 'w', encoding ='utf8') as json_file:
        json.dump(cond_probability_two, json_file) 
    with open('../vocab_generalization/cond_proba/total_one.json', 'w', encoding ='utf8') as json_file:
        json.dump(total_one, json_file) 
    with open('../vocab_generalization/cond_proba/total_two.json', 'w', encoding ='utf8') as json_file:
        json.dump(total_two, json_file) 