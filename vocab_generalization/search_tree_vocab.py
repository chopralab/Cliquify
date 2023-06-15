import sys
sys.path.append('../')
from rdkit import Chem
from rdkit import RDLogger
from collections import defaultdict
import pickle
import os

import networkx as nx
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from enumerate_random import get_fragments, get_fragments2
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol, data_to_mol

VOCAB_SIZE = len(os.listdir("../vocab_generalization/graph_vocab2"))

def get_common_vocabs():
    with open("../vocab_generalization/vocab_usage.txt") as f:
        vocab_counts = f.readlines()
        return [int(vocab_count.split(": ")[0]) for vocab_count in vocab_counts]

#----------------------------------------------------------------------
symbol_set = set()
formal_charge_set = set()
num_exp_Hs_set = set()
is_arom_set = set()

for i in get_common_vocabs():
    # mol = data_to_mol("../vocab_generalization/graph_vocab/{}".format(i))
    G = pickle.load(open("../vocab_generalization/graph_vocab2/{}".format(i), "rb"))

    symbol_set |= set(nx.get_node_attributes(G, "symbol").values())
    formal_charge_set |= set(nx.get_node_attributes(G, "formal_charge").values())
    num_exp_Hs_set |= set(nx.get_node_attributes(G, "num_explicit_hs").values())
    is_arom_set |= set(nx.get_node_attributes(G, "is_aromatic").values())

tree_vocab_byAtom = {}
for sym in symbol_set:
    tree_vocab_byAtom[sym] = {}
    for charge in formal_charge_set:
        tree_vocab_byAtom[sym][charge] = {}
        for num in num_exp_Hs_set:
            tree_vocab_byAtom[sym][charge][num] = {}
            for arom in is_arom_set:
                tree_vocab_byAtom[sym][charge][num][arom] = set()

for i in get_common_vocabs():
    G = pickle.load(open("../vocab_generalization/graph_vocab2/{}".format(i), "rb"))
    mol = nx_to_mol(G)
    for idx, data in G.nodes.data():
        sym, chrg, hs, arom, _ = data.values()
        tree_vocab_byAtom[sym][chrg][hs][arom].add(i)
    
#----------------------------------------------------------------------
tree_vocab_byBond = {
    0: set(),
    1: set(),
    3: {"2T": set(), "2F": set()},
    4: set()
}

for i in get_common_vocabs():
    G = pickle.load(open("../vocab_generalization/graph_vocab2/{}".format(i), "rb"))
    mol = data_to_mol("../vocab_generalization/graph_vocab2/{}".format(i))

    num_of_bonds = len(mol.GetBonds())
    if num_of_bonds < 2:
        tree_vocab_byBond[num_of_bonds].add(i)
    elif num_of_bonds == 3:
        ghost_list = list(nx.get_edge_attributes(G, "ghost").values())
        ghost_list.remove(False)
        ghost_bonds_id = sum(ghost_list)
        bondGhost = "2T" if ghost_bonds_id == 2 else "2F"
        tree_vocab_byBond[num_of_bonds][bondGhost].add(i)
    else:
        tree_vocab_byBond[4].add(i)


tree_vocab_byBondType = defaultdict(set)
"""{
1: set(), 1.5: set(), 2: set(), 3: set()
}
"""

for i in get_common_vocabs():
    mol = data_to_mol("../vocab_generalization/graph_vocab2/{}".format(i))
    
    for bond in mol.GetBonds():
        key = bond.GetBondTypeAsDouble()
        tree_vocab_byBondType[key].add(i)
    else:
        tree_vocab_byBondType[0.0].add(i)

# print(tree_vocab_byBondType.keys())

with open("../vocab_generalization/vocab_usage.txt") as f:
    vocab_counts = f.readlines()
    vocab_usage = {int(vocab_count.split(": ")[0]) : int(vocab_count.split(": ")[1].strip()) 
                   for vocab_count in vocab_counts}
    

tree_vocab_bySMILES = defaultdict(set)

for i in get_common_vocabs():
    mol = data_to_mol("../vocab_generalization/graph_vocab2/{}".format(i))
    
    for bond in mol.GetBonds():
        u, v = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        smiles = get_fragments(mol, [u, v])
        tree_vocab_bySMILES[smiles].add(i)
