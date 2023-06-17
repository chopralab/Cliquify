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


from enumerate_random import MolTreeNode, dfs_random_assemble, remove_edges_reset_idx, set_atommap, set_atommap_graph, get_fragments, get_fragments2 
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol, data_to_mol
from search_tree_vocab import get_common_vocabs, tree_vocab_byAtom, tree_vocab_byBond, tree_vocab_byBondType, tree_vocab_bySMILES, vocab_usage

def load_cond_proba():
    global cond_probability_one
    global cond_probability_two
    global total_one
    global total_two
    
    with open('../vocab_generalization/cond_proba/cond_probability_one.json', 'r') as json_file:
        cond_probability_one = json.load(json_file) 
    with open('../vocab_generalization/cond_proba/cond_probability_two.json', 'r') as json_file:
        cond_probability_two = json.load(json_file) 
    with open('../vocab_generalization/cond_proba/total_one.json', 'r') as json_file:
        total_one = json.load(json_file) 
    with open('../vocab_generalization/cond_proba/total_two.json', 'r') as json_file:
        total_two = json.load(json_file)

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
            return list((tree_vocab_byBond[3]["2F"] | tree_vocab_byBond[3]["2T"]) & suitable_vocabs)
        elif num_prev == 1 and num_current == 0: # complete single bonds -|-
            return list(tree_vocab_byBond[1] & suitable_vocabs)
    
    all_cliq_types = list(tree_vocab_byBond.keys())
    random_rank = []
    if num_current == 0: random_rank = [1,3]; np.random.shuffle(random_rank)
    elif num_current == 1: random_rank = [0, 1, 3]; np.random.shuffle(random_rank) 
    elif num_current == 3: random_rank = [0, 1, 3]; np.random.shuffle(random_rank)
    
    # elif num_current == 1: random_rank = all_cliq_types; np.random.shuffle(all_cliq_types) 
    # elif num_current == 3: random_rank = all_cliq_types; np.random.shuffle(all_cliq_types)
    # elif num_current == 4: random_rank = [1, 3]; np.random.shuffle(random_rank)
    
    for clq_size in random_rank:
        if clq_size == 3: filter_vocab = tree_vocab_byBond[clq_size]["2F"] | tree_vocab_byBond[clq_size]["2T"]
        else: filter_vocab = tree_vocab_byBond[clq_size]
        
        if filter_vocab & suitable_vocabs:
            # print(len(filter_vocab & suitable_vocabs), len(suitable_vocabs), len(filter_vocab & suitable_vocabs))
            return list(filter_vocab & suitable_vocabs)
        
    return list(suitable_vocabs)


def get_suitable_vocab(prev, root, is_leaf=False, filter_edge=False):
    suitable_vocabs = set()

    if root.graph.number_of_nodes() <= 2: # bond or atom
        for idx, data in root.graph.nodes.data():
            sym, chrg, hs, arom, _ = data.values()
            suitable_vocabs |= tree_vocab_byAtom[sym][chrg][hs][arom]
    else: # tri_clique or bridge or full compound
        for bond in root.tri_mol.GetBonds():
            u, v = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            smiles = get_fragments(root.tri_mol, [u, v])

            sym, chrg, hs, arom, _ = root.graph.nodes[u].values()
            vocab1 = tree_vocab_byAtom[sym][chrg][hs][arom]
            sym, chrg, hs, arom, _ = root.graph.nodes[v].values()
            vocab2 = tree_vocab_byAtom[sym][chrg][hs][arom]
            
            suitable_vocabs |= (vocab1 & vocab2 & tree_vocab_bySMILES[smiles])

    if is_leaf:
        leaf_filter = tree_vocab_byBond[1] | tree_vocab_byBond[3]["2F"]
        if suitable_vocabs & leaf_filter:
            suitable_vocabs = (suitable_vocabs & leaf_filter)

    if filter_edge:
        vocab_filter = set()
        for atom in root.tri_mol.GetAtoms():
            num_electrons = atom.GetAtomicNum() % 8 - 2 # valence electrons
            # print(num_electrons, atom.GetSymbol(), atom.GetTotalNumHs())
            for bond in atom.GetBonds(): 
                num_electrons -= bond.GetBondTypeAsDouble()

            if num_electrons % 2 == 0: # even
                vocab_filter |= tree_vocab_byBondType[2.0]
                vocab_filter |= tree_vocab_byBondType[1.0]
            elif num_electrons % 2 == 1: # odd
                if num_electrons >= 3: vocab_filter |= tree_vocab_byBondType[3.0]
                vocab_filter |= tree_vocab_byBondType[1.0]
            elif num_electrons % 2 == 0.5: # decimal (aromatic)
                vocab_filter |= tree_vocab_byBondType[1.5]

        if suitable_vocabs & vocab_filter and len(vocab_filter) < len(total_one.keys()):
            suitable_vocabs = suitable_vocabs & vocab_filter
        
    # # if nothing left, use this
    # if not suitable_vocabs:
    #     for idx, data in root.graph.nodes.data():
    #         sym, chrg, hs, arom, _ = data.values()
    #         suitable_vocabs |= tree_vocab_byAtom[sym][chrg][hs][arom]

    # remove all bridge[OPTIONAL]
    suitable_vocabs = suitable_vocabs - tree_vocab_byBond[4]

    return get_vocab_via_conn(prev, root, suitable_vocabs)

def get_proba(prev, current, suitable_vocabs):
    
    laplace = 5
    two_key = "{},{}".format(current.gid, prev.gid) if prev else ""
    one_key = "{}".format(current.gid)

    probs = []
    for vocab in suitable_vocabs:
        cond_two_key = "{}|{}".format(vocab, two_key)
        cond_one_key = "{}|{}".format(vocab, one_key)
        if cond_probability_two.get(cond_two_key) and total_two.get(two_key):
            probs.append((cond_probability_two.get(cond_two_key)*1.3+laplace)/(total_two.get(two_key)+laplace))
        elif cond_probability_one.get(cond_one_key) and total_one.get(one_key):
            probs.append((cond_probability_one.get(cond_one_key)+laplace)/(total_one.get(one_key)+laplace))
        else:
            probs.append((laplace)/(total_one.get(one_key)+laplace))

    # normalize
    total = sum([prob for prob in probs])
    return [prob/total for prob in probs]
    

#----------------------------------------------------#


MAX_ATOM_COUNT = 30
def random_sample_tree(prev, root, depth):
    global count
    global atom_count
    global edges_list
    global nodes_dict

    if atom_count > MAX_ATOM_COUNT: return
            

    # if root.graph.number_of_nodes() == 1:
    #     neighbors_count = np.random.randint(1, root.tri_mol.GetAtoms()[0].GetTotalNumHs()) # only spiro is possible // [tri/prev] - [atom/root] - [tri/neigh]
    # else:
    #     neighbors_count = np.random.randint(1, root.graph.number_of_nodes()) # max of 4 neigbors

    neighbors_count = np.random.randint(1, 3)

    atom_count += root.graph.number_of_nodes() - random.choice([1, 2])

    for i in range(neighbors_count):
        # if atom_count > MAX_ATOM_COUNT: return

        count += 1
        topology = np.random.rand()

        # if topology > 0.05: 
        if atom_count < MAX_ATOM_COUNT:
            suitable_vocabs = get_suitable_vocab(prev, root, is_leaf=False)
            if not suitable_vocabs: return # if unable to generate suitable terminate
            probs = get_proba(prev, root, suitable_vocabs)
            idx = np.random.choice(suitable_vocabs, p=probs)

            neigh, neigh.gid = MolTreeNode(vocab(idx)), idx
            atom_count += neigh.graph.number_of_nodes() - random.choice([1, 2])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh

            random_sample_tree(root, neigh, depth+1)

        else:
            suitable_vocabs = get_suitable_vocab(prev, root, is_leaf=True)
            if not suitable_vocabs: return # if unable to generate suitable terminate
            probs = get_proba(prev, root, suitable_vocabs)
            idx = np.random.choice(suitable_vocabs, p=probs)

            neigh, neigh.gid = MolTreeNode(vocab(idx)), idx
            atom_count += neigh.graph.number_of_nodes() - np.random.choice([1, 2])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh


def save_tree(nodes_dict, edges_list, idx=7000, folder="../vocab_generalization/tree"):
    plt.clf()
    graph = nx.Graph()
    for u, v in edges_list:
        u_gid, v_gid = nodes_dict[u].gid, nodes_dict[v].gid
        graph.add_node(u, smiles=Chem.MolToSmiles(vocab(u_gid)))
        graph.add_node(v, smiles=Chem.MolToSmiles(vocab(v_gid)))
        graph.add_edge(u, v)
    
    with open("{}/show{}_tree.pkl".format(folder, idx), 'wb') as handle:
        pickle.dump(nodes_dict, handle)

    pos = nx.spring_layout(graph)
    nx.draw(graph, pos)
    node_labels = nx.get_node_attributes(graph, "smiles")
    node_labels = {k : "     ({})".format(v) for k, v in node_labels.items()}
    nx.draw_networkx_labels(graph, pos, node_labels)
    # nx.draw(graph, pos, node_color="yellow", with_labels=True)
    plt.savefig("{}/show{}_tree.png".format(folder, idx))
    plt.clf()


if __name__ == '__main__':
    np.random.seed(42)

    load_cond_proba()

    gen_smiles_list = []
    for sample_idx in range(1000):

        #--------------------------------------------------
        count = 0
        atom_count = 0

        probs = [total_one[str(vocab)]/sum(total_one.values()) for vocab in get_start_vocabs()]
        probs = [prob/sum(probs) for prob in probs]
        start_idx = np.random.choice(get_start_vocabs())
        root, root.gid = MolTreeNode(vocab(start_idx)), start_idx
        root.idx = 0
        nodes_dict = {root.idx: root}
        edges_list = []

        random_sample_tree(None, root, 0)

        for x,y in edges_list:
            nodes_dict[x].add_neighbor(nodes_dict[y])
            nodes_dict[y].add_neighbor(nodes_dict[x])
        
        
        #----
        save_tree(nodes_dict, edges_list, idx=sample_idx)
        #----

        list_of_nodes = list(nodes_dict.values())
        for i,node in enumerate(list_of_nodes):
            node.nid = i + 1
            if len(node.neighbors) > 1: #Leaf node mol is not marked
                set_atommap(node.tri_mol, node.nid)
                set_atommap_graph(node.graph, node.nid)
            node.is_leaf = (len(node.neighbors) == 1)

        #---------------------------------------------------

        possible_smiles = ""
        possible_graph = None
        for start_idx in range(len(list_of_nodes)):

            root_idx = start_idx
            root = list_of_nodes[root_idx]

            cur_graph = root.graph.copy()
            global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
            global_amap[root_idx+1] = {node:node for node in cur_graph.nodes()}

            dfs_random_assemble(cur_graph, global_amap, [], root, None, print_out=False)

            final_graph = remove_edges_reset_idx(cur_graph)

            # remove unconnected nodes remnants 
            final_graph.remove_nodes_from(list(nx.isolates(final_graph)))

            #----------------------------------------------#
            cur_mol = nx_to_mol(mol_to_nx(nx_to_mol(final_graph)))
            dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)
            dec_mol  = Chem.MolFromSmiles(dec_smiles)

            if cur_mol and dec_mol and dec_smiles:
                if len(dec_smiles) > len(possible_smiles):
                    possible_smiles = dec_smiles
                    possible_graph = cur_graph


        # if dec_smiles not in gen_smiles_list:
        if possible_smiles and possible_graph and ("." not in possible_smiles) and (possible_smiles not in gen_smiles_list):

            draw_mol(possible_graph, sample_idx, ["symbol", "bond_type", "color"], folder="tree", label="dec_graph")
            print(possible_smiles, sample_idx)
            # print("num_of_nodes", len(list_of_nodes))      
            gen_smiles_list.append(possible_smiles)


    with open("gen_mol_count_rand42_new.txt", "w") as myfile:
        for gen_smiles in gen_smiles_list:
            myfile.writelines("{}\n".format(gen_smiles))
        # raise
