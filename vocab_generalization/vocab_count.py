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


from enumerate3 import get_fragments, get_mol2, get_smarts_fragments, get_mol2, get_triangulated_graph
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, node_equal_iso, ring_edge_equal_iso, mol_to_data



with open("../zinc/all.txt") as f:
    smiles_list = f.readlines()


total_vocab_smiles = set()
total_vocab_smiles_smarts = set()
# total_vocab_inchi = set()
# total_vocab_graph = set()
total_vocab_graph = []
def vocab_count(idx):
    chosen_smiles = smiles_list[idx]
    mol = Chem.MolFromSmiles(chosen_smiles)

    # remove stereochemistry
    gold_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
    mol = Chem.MolFromSmiles(gold_smiles)

    cliques, _, _ = tree_decomp(mol)

    for i, clique in enumerate(cliques):
        if isinstance(clique, tuple):
            # if len(clique) > 4: continue
            fragment_smiles = get_fragments(mol, clique)
            canon_smiles = Chem.CanonSmiles(fragment_smiles)
            fragment_smarts = get_smarts_fragments(mol, clique)
            tri_mol, fragment_graph = get_triangulated_graph(get_mol2(mol, clique))
            # inchiKey = Chem.MolToInchiKey(get_mol2(mol, clique))            
            duplicate=False
            for G in total_vocab_graph:
                if nx.is_isomorphic(G, fragment_graph, \
                                    node_match=node_equal_iso, \
                                    edge_match=ring_edge_equal_iso):
                    duplicate = True
                    break
            if duplicate: continue

            # print(canon_smiles == fragment_smiles)

            total_vocab_smiles.add(fragment_smiles)
            total_vocab_smiles_smarts.add(fragment_smiles + "|" + canon_smiles + "|" + fragment_smarts)
            # total_vocab_inchi.add(inchiKey)
            total_vocab_graph.append(fragment_graph)

    if idx % 1000 == 0:
        print("smiles", len(total_vocab_smiles))
        print("smiles_smarts", len(total_vocab_smiles_smarts))
        # print("Inchi Key", len(total_vocab_inchi))
        print("graph", len(total_vocab_graph))
        print()

for i in range(len(smiles_list)):
    vocab_count(i)

print()
total_vocab_smiles = list(total_vocab_smiles)
print("total_vocab_smiles", len(total_vocab_smiles))
print("total_vocab_graph", len(total_vocab_graph))

# save graph vocabs
for i, vocab_G in enumerate(total_vocab_graph):
    pickle.dump(vocab_G, open("../vocab_generalization/graph_vocab/{}".format(i), "wb"))

