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


from enumerate3 import get_fragments, get_fragments2, get_smarts_fragments
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso


with open("../zinc/all.txt") as f:
    smiles_list = f.readlines()


total_vocab = set()
def vocab_count(idx):
    chosen_smiles = smiles_list[idx]
    mol = Chem.MolFromSmiles(chosen_smiles)
    cliques, _, _ = tree_decomp(mol)

    for i, clique in enumerate(cliques):
        if isinstance(clique, tuple):
            if len(clique) > 4: continue
            fragment_smiles = get_fragments(mol, clique)
            # fragment_smiles = get_smarts_fragments(mol, clique)

            total_vocab.add(fragment_smiles)

    if idx % 2000 == 0: 
        print("OK", len(total_vocab))

# with concurrent.futures.ThreadPoolExecutor() as executor:
#     executor.map(vocab_count, [i for i in range(len(smiles_list))])

for i in range(len(smiles_list)):
    vocab_count(i)

total_vocab = list(total_vocab)
print()
print("total_vocab", len(total_vocab))
print("unkekulizable_set", unkekulizable_set)

with open("vocab_smarts.txt", "a") as myfile:
    for vocab in total_vocab:       
        myfile.writelines("{}\n".format(vocab))
