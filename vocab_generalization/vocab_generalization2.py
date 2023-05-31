import sys
sys.path.append('../')
from rdkit import Chem
import networkx as nx
import numpy as np
import networkx.algorithms.isomorphism as iso
import random
import pickle

from enumerate_random import MolTreeNode, dfs_random_assemble, remove_edges_reset_idx, set_atommap, set_atommap_graph
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol, data_to_mol
from search_tree_vocab import get_common_vocabs, tree_vocab_byAtom, tree_vocab_byBond, vocab_usage


def vocab(idx):
    return data_to_mol("../vocab_generalization/graph_vocab2/{}".format(idx))

def graph(idx):
    return pickle.load(open("../vocab_generalization/graph_vocab2/{}".format(idx), "rb"))

def get_start_vocabs():
    suitable_vocabs = set()
    suitable_vocabs |= tree_vocab_byBond[1]
    suitable_vocabs |= tree_vocab_byBond[3]["2F"]
    suitable_vocabs |= tree_vocab_byBond[4]
    return list(suitable_vocabs)

def get_vocab_via_conn(prev, current, suitable_vocabs):
    
    num_current =current.graph.number_of_nodes()
    if prev:
        num_prev = prev.graph.number_of_nodes()

        if num_prev == 3 and num_current == 0: # complete spiro
            return list(tree_vocab_byBond[3]["2F"] & suitable_vocabs)
        elif num_prev == 1 and num_current == 0: # complete single bonds -|-
            return list(tree_vocab_byBond[1] & suitable_vocabs)
    
    all_cliq_types = list(tree_vocab_byBond.keys())
    random_rank = []
    if num_current == 0: random_rank = [1,3]; np.random.shuffle(random_rank)
    elif num_current == 1: random_rank = all_cliq_types; np.random.shuffle(all_cliq_types) 
    elif num_current == 3: random_rank = all_cliq_types; np.random.shuffle(all_cliq_types)
    elif num_current == 4: random_rank = [1, 3]; np.random.shuffle(random_rank)
    
    for clq_size in random_rank:
        if clq_size == 3: filter_vocab = tree_vocab_byBond[clq_size]["2F"] | tree_vocab_byBond[clq_size]["2T"]
        else: filter_vocab = tree_vocab_byBond[clq_size]
        
        if filter_vocab & suitable_vocabs:
            return list(filter_vocab & suitable_vocabs)
        
    return list(suitable_vocabs)


def get_suitable_vocab(prev, root, is_leaf=False):
    suitable_vocabs = set()
    if not is_leaf:
        for idx, data in root.graph.nodes.data():
            sym, chrg, hs, arom, _ = data.values()
            suitable_vocabs |= tree_vocab_byAtom[sym][chrg][hs][arom]
    else:
        suitable_vocabs |= tree_vocab_byBond[1]
        suitable_vocabs |= tree_vocab_byBond[3]["2F"]

    return get_vocab_via_conn(prev, root, suitable_vocabs)

def get_proba(suitable_vocabs):
    total = sum([vocab_usage[vocab] for vocab in suitable_vocabs])
    return [vocab_usage[vocab]/total for vocab in suitable_vocabs]
    

#----------------------------------------------------#

common_vocabs = get_common_vocabs()

np.random.seed(42)

MAX_ATOM_COUNT = 30

def random_sample_tree(prev, root, depth):
    global count
    global atom_count
    global edges_list
    global nodes_dict

    if atom_count > MAX_ATOM_COUNT: return
            

    if root.graph.number_of_nodes() == 1:
        neighbors_count = np.random.randint(1, root.tri_mol.GetAtoms()[0].GetTotalNumHs()) # only spiro is possible // [tri/prev] - [atom/root] - [tri/neigh]
    else:
        neighbors_count = np.random.randint(1, root.graph.number_of_nodes()) # max of 4 neigbors

    atom_count += root.graph.number_of_nodes() - random.choice([1, 2])

    for i in range(neighbors_count):
        if atom_count > MAX_ATOM_COUNT: return

        count += 1
        topology = np.random.rand()

        if topology > 0.5:
            suitable_vocabs = get_suitable_vocab(prev, root, is_leaf=False)
            probs = get_proba(suitable_vocabs)
            idx = np.random.choice(suitable_vocabs, p=probs)

            neigh = MolTreeNode(vocab(idx))
            atom_count += neigh.graph.number_of_nodes() - random.choice([1, 2])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh

            random_sample_tree(root, neigh, depth+1)

        else:
            suitable_vocabs = get_suitable_vocab(prev, root, is_leaf=True)
            probs = get_proba(suitable_vocabs)
            idx = np.random.choice(suitable_vocabs, p=probs)

            neigh = MolTreeNode(vocab(idx))
            atom_count += neigh.graph.number_of_nodes() - np.random.choice([1, 2])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh


gen_smiles_list = []
for _ in range(1000):

    #--------------------------------------------------
    count = 0
    atom_count = 0

    start_idx = np.random.choice(get_start_vocabs())
    root = MolTreeNode(vocab(start_idx))
    root.idx = 0
    nodes_dict = {root.idx: root}
    edges_list = []

    random_sample_tree(None, root, 0)

    for x,y in edges_list:
        nodes_dict[x].add_neighbor(nodes_dict[y])
        nodes_dict[y].add_neighbor(nodes_dict[x])
    
    list_of_nodes = list(nodes_dict.values())

    for i,node in enumerate(list_of_nodes):
        node.nid = i + 1
        if len(node.neighbors) > 1: #Leaf node mol is not marked
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

    if cur_mol and dec_mol and dec_smiles:
        print(dec_smiles)

        if dec_smiles not in gen_smiles_list:
            gen_smiles_list.append(dec_smiles)

    # else:
        # print(None)


with open("gen_mol_count_rand42_modified.txt", "a") as myfile:
    for gen_smiles in gen_smiles_list:
        myfile.writelines("{}\n".format(gen_smiles))
    # raise
